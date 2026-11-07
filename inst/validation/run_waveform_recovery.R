source(file.path("inst", "validation", "waveform_recovery.R"), local = TRUE)

wv_execute <- function(data, settings, sources = "oracle", keep_example = FALSE) {
  fit <- if ("learned" %in% sources) wv_hrf_fit(data) else NULL
  metrics <- leaks <- diagnostics <- examples <- list()
  for (source in sources) {
    kernels <- wv_kernels(data,source,fit)
    cache <- new.env(parent=emptyenv())
    for (i in seq_len(nrow(settings))) {
      setting <- settings[i,,drop=FALSE]; runs <- list()
      for (r in seq_along(data$bold)) {
        k <- kernels[[r]]
        answer <- wv_reconstruct(data$bold[[r]], k$h,k$rho,data$masks[[r]],setting,data$groups,cache)
        counter <- wv_reconstruct(data$without_previous[[r]], k$h,k$rho,data$masks[[r]],setting,data$groups,cache,
          scale_reference=data$bold[[r]])
        score <- wv_metrics(answer,counter,data,r)
        runs[[r]] <- score$metrics
        leaks[[length(leaks)+1L]] <- data.frame(source=source,id=setting$id,run=r,score$leakage)
        diagnostic <- rbind(data.frame(fit="observed",answer$diagnostic),data.frame(fit="counterfactual",counter$diagnostic))
        diagnostics[[length(diagnostics)+1L]] <- data.frame(source=source,id=setting$id,run=r,diagnostic)
        if (keep_example && r==1L) examples[[paste(source,setting$id,sep="|")]] <- answer$activity
      }
      values <- rowMeans(do.call(cbind,runs))
      metrics[[length(metrics)+1L]] <- data.frame(source=source,id=setting$id,method=setting$method,
        lambda=setting$lambda,spatial=setting$spatial,as.list(values))
    }
  }
  list(metrics=do.call(rbind,metrics),leakage=do.call(rbind,leaks),diagnostics=do.call(rbind,diagnostics),
    example=if(keep_example) list(data=data,activity=examples) else NULL,
    hrf=if(is.null(fit)) NULL else lapply(seq_along(data$bold),function(r) data.frame(run=r,
      voxel=seq_len(data$config$voxels),kernel=fit$selections[[r]]$kernel_index,
      status=fit$selections[[r]]$status,rho=fit$noise[[r]]$rho,
      kernel_nmse=colSums((fit$kernels[,fit$selections[[r]]$kernel_index,drop=FALSE]-data$kernels)^2)/colSums(data$kernels^2))))
}

# Select a global setting per method using mean continuous NMSE over independent
# development datasets. Evaluation truth never enters selection. Nonconverged
# candidates are ineligible. The existing lambda=.01 diff2 is retained separately.
wv_select <- function(metrics) {
  if (any(metrics$source!="oracle")) stop("Selection requires oracle-HRF development only")
  grid <- wv_grid(); scores <- vapply(grid$id,function(id) mean(metrics$neural_nmse[metrics$id==id]),numeric(1))
  eligible <- vapply(grid$id,function(id) all(metrics$convergence[metrics$id==id]==1) && any(metrics$id==id),logical(1))
  selected <- lapply(unique(grid$method),function(method) {
    i <- which(grid$method==method & eligible & is.finite(scores))
    if (!length(i)) stop("No converged development candidate for ",method)
    grid[i[which.min(scores[i])],,drop=FALSE]
  })
  out <- do.call(rbind,selected)
  baseline <- grid[grid$method=="diff2" & grid$lambda==.01,,drop=FALSE]
  baseline$id <- "diff2_fixed/0.010"; baseline$method <- "diff2_fixed"
  # The solver sees diff2 for both fixed and development-selected references.
  rbind(out,baseline)
}

run_waveform_recovery <- function(out_dir, phase=c("development","evaluation"), repetitions=NULL,
    scenario_names=NULL, selection_path=NULL, cache_dir="/tmp/waveform-study/cpp-cache") {
  phase <- match.arg(phase)
  if (is.null(repetitions)) repetitions <- if (phase=="development") 4L else 20L
  stopifnot(repetitions>=1L,repetitions==as.integer(repetitions),repetitions<10000L)
  if (dir.exists(out_dir) && length(list.files(out_dir,all.files=TRUE,no..=TRUE))) stop("Use an empty output directory")
  if (is.null(scenario_names)) scenario_names <- if(phase=="development") c("brief","blocks") else names(wv_scenarios())
  scenarios <- wv_scenarios()[scenario_names]
  if (!length(scenarios) || any(vapply(scenarios,is.null,logical(1)))) stop("Unknown or empty scenarios")
  if (phase=="development" && !all(scenario_names %in% c("brief","blocks"))) stop("Held-out family used in development")
  if (phase=="evaluation") {
    if (is.null(selection_path)) stop("Supply frozen development selection")
    selection <- readRDS(selection_path)
    if (!identical(selection$source_md5,tools::md5sum(wv_sources()))) stop("Development sources changed")
    settings <- selection$settings
  } else settings <- wv_grid()
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  protocol <- list(phase=phase,repetitions=repetitions,scenarios=scenarios,settings=settings,
    source_md5=tools::md5sum(wv_sources()),
    selection_md5=if(!is.null(selection_path)) tools::md5sum(selection_path) else NULL,
    seed_base=if(phase=="development") 24000000L else 26000000L,
    seed_formula="base + original_scenario_index * 10000 + replicate",
    sources=if(phase=="development") "oracle" else c("oracle","learned"),
    primary="oracle HRF/rho; trial waveform NMSE and paired previous-source excess leakage",
    selection="mean continuous neural NMSE over brief/block development only; all-voxel convergence required",
    offset="evaluation-only single constant alignment per voxel over usable interior; no lag/scale correction",
    solver=list(admm_rho=.01,max_iterations=12000L,absolute_tolerance=1e-6,relative_tolerance=1e-5,ridge=1e-8),
    threads=Sys.getenv(c("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS")),rng_kind=RNGkind())
  saveRDS(protocol,file.path(out_dir,"settings.rds"))
  wv_load_solver(cache_dir)
  metrics <- leakage <- diagnostics <- hrfs <- list(); start <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) for (replicate in seq_len(repetitions)) {
    seed <- protocol$seed_base + match(scenario,names(wv_scenarios()))*10000L + replicate
    message(scenario," ",replicate,"/",repetitions)
    data <- wv_generate(scenarios[[scenario]],seed)
    execute_settings <- settings
    execute_settings$method[execute_settings$method=="diff2_fixed"] <- "diff2"
    result <- wv_execute(data,execute_settings,protocol$sources,keep_example=replicate==1L)
    result$metrics$method <- settings$method[match(result$metrics$id,settings$id)]
    prefix <- data.frame(scenario=scenario,replicate=replicate,seed=seed)
    metrics[[length(metrics)+1L]] <- cbind(prefix,result$metrics)
    leakage[[length(leakage)+1L]] <- cbind(prefix,result$leakage)
    diagnostics[[length(diagnostics)+1L]] <- cbind(prefix,result$diagnostics)
    if(!is.null(result$hrf)) hrfs[[length(hrfs)+1L]] <- cbind(prefix,do.call(rbind,result$hrf))
    if(replicate==1L) saveRDS(result$example,file.path(out_dir,paste0("example_",scenario,".rds")),compress="xz")
    saveRDS(do.call(rbind,metrics),file.path(out_dir,"replicate_metrics.rds"))
    saveRDS(do.call(rbind,leakage),file.path(out_dir,"leakage.rds"))
    saveRDS(do.call(rbind,diagnostics),file.path(out_dir,"solver_diagnostics.rds"),compress="xz")
    if(length(hrfs)) saveRDS(do.call(rbind,hrfs),file.path(out_dir,"hrf_diagnostics.rds"))
    message("elapsed ",round(proc.time()[["elapsed"]]-start,1)," s")
  }
  if(phase=="development") saveRDS(list(settings=wv_select(do.call(rbind,metrics)),
    source_md5=protocol$source_md5,development_seed_base=protocol$seed_base,
    development_metrics_md5=tools::md5sum(file.path(out_dir,"replicate_metrics.rds"))),file.path(out_dir,"selection.rds"))
  writeLines(c(capture.output(sessionInfo()),paste("Elapsed seconds",proc.time()[["elapsed"]]-start)),file.path(out_dir,"environment.txt"))
  invisible(do.call(rbind,metrics))
}

if (sys.nframe()==0L) {
  .libPaths(c(Sys.getenv("FMRI_VALIDATION_LIBRARY",unset=.libPaths()[1L]),.libPaths()))
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)<2L) stop("Supply output directory and phase")
  run_waveform_recovery(args[1L],args[2L],if(length(args)>2L) as.integer(args[3L]) else NULL,
    if(length(args)>3L && nzchar(args[4L])) strsplit(args[4L],",",fixed=TRUE)[[1L]] else NULL,
    if(length(args)>4L) args[5L] else NULL)
}
