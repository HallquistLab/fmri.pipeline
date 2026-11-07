# Resume one scenario from a STOPPED evaluation shard into a new directory.
# Seeds, settings, source versions and completed dataset records remain fixed.
source(file.path("inst","validation","run_waveform_recovery.R"),local=TRUE)

resume_waveform_recovery<-function(input_dir,scenario,out_dir,cache_dir="/tmp/waveform-study/cpp-cache",workers=1L) {
  stopifnot(workers>=1L,workers==as.integer(workers))
  protocol<-readRDS(file.path(input_dir,"settings.rds"))
  if(protocol$phase!="evaluation" || !scenario%in%names(protocol$scenarios))stop("Unknown evaluation scenario")
  if(!identical(protocol$source_md5,tools::md5sum(wv_sources())))stop("Study sources changed")
  if(dir.exists(out_dir) && length(list.files(out_dir,all.files=TRUE,no..=TRUE)))stop("Use an empty output directory")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  protocol$scenarios<-protocol$scenarios[scenario]
  saveRDS(protocol,file.path(out_dir,"settings.rds"))
  read_subset<-function(name) {
    path<-file.path(input_dir,paste0(name,".rds"))
    if(!file.exists(path))return(NULL)
    x<-readRDS(path);x[x$scenario==scenario,,drop=FALSE]
  }
  metrics<-read_subset("replicate_metrics");leakage<-read_subset("leakage")
  diagnostics<-read_subset("solver_diagnostics");hrfs<-read_subset("hrf_diagnostics")
  complete<-if(is.null(metrics))integer()else sort(unique(metrics$replicate))
  expected<-expand.grid(source=protocol$sources,id=protocol$settings$id,stringsAsFactors=FALSE)
  key<-function(x)paste(x$source,x$id)
  for(rep in complete) {
    x<-metrics[metrics$replicate==rep,,drop=FALSE]
    if(anyDuplicated(key(x)) || !setequal(key(x),key(expected)))stop("Incomplete saved dataset")
    if(any(x$seed!=protocol$seed_base+match(scenario,names(wv_scenarios()))*10000L+rep))stop("Incorrect saved seed")
  }
  example_path<-file.path(input_dir,paste0("example_",scenario,".rds"))
  if(1L%in%complete) {
    if(!file.exists(example_path))stop("Missing saved example")
    file.copy(example_path,file.path(out_dir,basename(example_path)))
  }
  wv_load_solver(cache_dir);started<-proc.time()[["elapsed"]]
  missing<-setdiff(seq_len(protocol$repetitions),complete)
  message(scenario,": retaining ",length(complete)," completed datasets; fitting ",length(missing))
  append<-function(old,new)if(is.null(old)||!nrow(old))new else rbind(old,new)
  flush<-function() {
    saveRDS(metrics,file.path(out_dir,"replicate_metrics.rds"))
    saveRDS(leakage,file.path(out_dir,"leakage.rds"))
    saveRDS(diagnostics,file.path(out_dir,"solver_diagnostics.rds"),compress="xz")
    saveRDS(hrfs,file.path(out_dir,"hrf_diagnostics.rds"))
  }
  if(length(complete))flush()
  calculate<-function(rep) {
    seed<-protocol$seed_base+match(scenario,names(wv_scenarios()))*10000L+rep
    message(scenario," ",rep,"/",protocol$repetitions)
    data<-wv_generate(protocol$scenarios[[scenario]],seed)
    settings<-protocol$settings;settings$method[settings$method=="diff2_fixed"]<-"diff2"
    answer<-wv_execute(data,settings,protocol$sources,keep_example=rep==1L)
    answer$metrics$method<-protocol$settings$method[match(answer$metrics$id,protocol$settings$id)]
    prefix<-data.frame(scenario=scenario,replicate=rep,seed=seed)
    list(replicate=rep,metrics=cbind(prefix,answer$metrics),leakage=cbind(prefix,answer$leakage),
      diagnostics=cbind(prefix,answer$diagnostics),hrf=cbind(prefix,do.call(rbind,answer$hrf)),example=answer$example)
  }
  retain<-function(answer) {
    if(inherits(answer,"try-error"))stop("A resumed dataset failed; inspect worker messages")
    metrics<<-append(metrics,answer$metrics);leakage<<-append(leakage,answer$leakage)
    diagnostics<<-append(diagnostics,answer$diagnostics);hrfs<<-append(hrfs,answer$hrf)
    if(answer$replicate==1L)saveRDS(answer$example,file.path(out_dir,paste0("example_",scenario,".rds")),compress="xz")
    flush();message("elapsed ",round(proc.time()[["elapsed"]]-started,1)," s")
  }
  if(workers==1L)for(rep in missing)retain(calculate(rep))else {
    # Each generator sets its explicit dataset seed. No fork RNG streams,
    # truth-based selection, or shared native plans cross worker boundaries.
    answers<-parallel::mclapply(missing,calculate,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
    for(answer in answers)retain(answer)
  }
  original_environment<-file.path(input_dir,"environment.txt")
  writeLines(c(capture.output(sessionInfo()),paste("Resumed stopped shard",input_dir),
    paste("Retained replicates",paste(complete,collapse=",")),
    paste("Resume source MD5",unname(tools::md5sum("inst/validation/resume_waveform_recovery.R"))),
    paste("Resume workers",workers),
    paste("Additional elapsed seconds",proc.time()[["elapsed"]]-started),
    if(file.exists(original_environment))readLines(original_environment)else "Original shard stopped before final environment write"),
    file.path(out_dir,"environment.txt"))
  invisible(metrics)
}

if(sys.nframe()==0L) {
  .libPaths(c(Sys.getenv("FMRI_VALIDATION_LIBRARY",unset=.libPaths()[1L]),.libPaths()))
  args<-commandArgs(trailingOnly=TRUE)
  if(!length(args)%in%c(3L,4L))stop("Supply stopped input shard, scenario, new output directory, optional workers")
  resume_waveform_recovery(args[1L],args[2L],args[3L],workers=if(length(args)>3L)as.integer(args[4L])else 1L)
}
