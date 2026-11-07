source(file.path("inst","validation","run_waveform_recovery.R"),local=TRUE)

wv_key<-function(x,fields)do.call(paste,c(x[fields],sep="|"))
wv_collect<-function(input_dir) {
  paths<-if(file.exists(file.path(input_dir,"settings.rds")))file.path(input_dir,"settings.rds")else
    list.files(input_dir,pattern="^settings[.]rds$",recursive=TRUE,full.names=TRUE)
  if(!length(paths))stop("No frozen settings")
  protocols<-lapply(paths,readRDS);metrics<-leakage<-diagnostics<-hrfs<-list();seen<-character()
  invariant<-setdiff(names(protocols[[1]]),"scenarios")
  for(i in seq_along(paths)) {
    p<-protocols[[i]];dir<-dirname(paths[i]);scenarios<-names(p$scenarios)
    if(!identical(p$source_md5,tools::md5sum(wv_sources())))stop("Study sources changed")
    if(!identical(p[invariant],protocols[[1]][invariant]))stop("Incompatible protocols")
    if(any(scenarios%in%seen))stop("Overlapping scenario shards")
    seen<-c(seen,scenarios)
    expected<-expand.grid(scenario=scenarios,replicate=seq_len(p$repetitions),source=p$sources,id=p$settings$id,
      stringsAsFactors=FALSE)
    d<-readRDS(file.path(dir,"replicate_metrics.rds"));fields<-c("scenario","replicate","source","id")
    if(anyDuplicated(wv_key(d,fields)) || !setequal(wv_key(d,fields),wv_key(expected,fields)))stop("Incomplete or duplicate metrics")
    verify_seed<-function(x) {
      seed<-p$seed_base+match(x$scenario,names(wv_scenarios()))*10000L+x$replicate
      if(anyNA(seed)||any(x$seed!=seed))stop("Incorrect dataset seeds")
    }
    verify_seed(d)
    s<-readRDS(file.path(dir,"solver_diagnostics.rds"));verify_seed(s)
    e<-do.call(rbind,lapply(scenarios,function(scenario) {
      config<-p$scenarios[[scenario]]
      expand.grid(scenario=scenario,replicate=seq_len(p$repetitions),source=p$sources,id=p$settings$id,
        run=seq_len(config$runs),voxel=seq_len(config$voxels),fit=c("observed","counterfactual"),stringsAsFactors=FALSE)
    }))
    fields<-c("scenario","replicate","source","id","run","voxel","fit")
    if(anyDuplicated(wv_key(s,fields)) || !setequal(wv_key(s,fields),wv_key(e,fields)))stop("Incomplete solver diagnostics")
    l<-readRDS(file.path(dir,"leakage.rds"));verify_seed(l)
    leakage_expected<-list()
    for(scenario in scenarios) {
      config<-p$scenarios[[scenario]];pairs<-length(seq(36,config$n-50,by=40))
      le<-expand.grid(scenario=scenario,replicate=seq_len(p$repetitions),source=p$sources,id=p$settings$id,
        run=seq_len(config$runs),pair=seq_len(pairs),stringsAsFactors=FALSE)
      gaps<-rep(c(3,6,10),length.out=pairs)
      le$gap<-gaps[(le$pair+le$run-2L)%%pairs+1L];le$null<-le$pair%%2L==1L
      leakage_expected[[length(leakage_expected)+1L]]<-le
    }
    fields<-c("scenario","replicate","source","id","run","gap","null")
    observed_counts<-table(wv_key(l,fields));expected_counts<-table(wv_key(do.call(rbind,leakage_expected),fields))
    if(!identical(names(observed_counts),names(expected_counts)) ||
        !identical(as.vector(observed_counts),as.vector(expected_counts)))stop("Incomplete or duplicate leakage diagnostics")
    metrics[[i]]<-d;leakage[[i]]<-l;diagnostics[[i]]<-s
    if("learned"%in%p$sources) {
      hrfs[[i]]<-readRDS(file.path(dir,"hrf_diagnostics.rds"));verify_seed(hrfs[[i]])
      fields<-c("scenario","replicate","run","voxel")
      expected_hrf<-unique(e[fields])
      if(anyDuplicated(wv_key(hrfs[[i]],fields)) || !setequal(wv_key(hrfs[[i]],fields),wv_key(expected_hrf,fields)))
        stop("Incomplete HRF diagnostics")
    }
  }
  list(paths=paths,protocols=protocols,metrics=do.call(rbind,metrics),leakage=do.call(rbind,leakage),
    diagnostics=do.call(rbind,diagnostics),hrf=do.call(rbind,hrfs))
}

wv_summarize<-function(data,fields,measures) {
  groups<-split(seq_len(nrow(data)),wv_key(data,fields))
  do.call(rbind,lapply(groups,function(rows)do.call(rbind,lapply(measures,function(metric) {
    x<-data[[metric]][rows];x<-x[is.finite(x)];n<-length(x);se<-if(n>1L)stats::sd(x)/sqrt(n)else NA_real_
    half<-if(n>1L)stats::qt(.975,n-1L)*se else NA_real_
    data.frame(data[rows[1L],fields,drop=FALSE],metric=metric,n=n,mean=if(n)mean(x)else NA_real_,
      mc_se=se,lower=if(n)mean(x)-half else NA_real_,upper=if(n)mean(x)+half else NA_real_)
  }))))
}

wv_paired<-function(data,baseline="diff2_fixed",converged_only=FALSE) {
  metrics<-setdiff(names(data),c("scenario","replicate","seed","source","id","method","lambda","spatial"))
  base<-data[data$method==baseline,];other<-data[data$method!=baseline,]
  joined<-merge(other,base,by=c("scenario","replicate","seed","source"),suffixes=c("","_base"))
  if(converged_only)joined<-joined[joined$convergence==1 & joined$convergence_base==1,,drop=FALSE]
  for(metric in metrics)joined[[metric]]<-joined[[metric]]-joined[[paste0(metric,"_base")]]
  wv_summarize(joined,c("scenario","source","method"),metrics)
}

freeze_waveform_development<-function(input_dir,out_dir,selection_path) {
  d<-wv_collect(input_dir)
  if(d$protocols[[1]]$phase!="development" || !setequal(unique(d$metrics$scenario),c("brief","blocks")))
    stop("Selection requires complete brief/block development")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  selection<-list(settings=wv_select(d$metrics),source_md5=tools::md5sum(wv_sources()),
    development_seed_base=d$protocols[[1]]$seed_base,
    development_metrics_md5=tools::md5sum(file.path(dirname(d$paths),"replicate_metrics.rds")))
  saveRDS(selection,selection_path);saveRDS(selection,file.path(out_dir,"selection.rds"))
  utils::write.csv(selection$settings,file.path(out_dir,"selection.csv"),row.names=FALSE)
  utils::write.csv(d$metrics,file.path(out_dir,"development_metrics.csv"),row.names=FALSE)
  utils::write.csv(wv_summarize(d$metrics,c("id","method","lambda","spatial"),c("neural_nmse","convergence","elapsed")),
    file.path(out_dir,"development_summary.csv"),row.names=FALSE)
  utils::write.csv(stats::aggregate(converged~id,d$diagnostics,function(x)c(n=length(x),failures=sum(!x))),
    file.path(out_dir,"development_convergence.csv"),row.names=FALSE)
  saveRDS(d$protocols,file.path(out_dir,"development_settings.rds"))
  writeLines(unlist(lapply(d$paths,function(p)c(paste("Shard",dirname(p)),
    readLines(file.path(dirname(p),"environment.txt"))))),file.path(out_dir,"development_environment.txt"))
  invisible(selection)
}

report_waveform_recovery<-function(input_dir,out_dir) {
  d<-wv_collect(input_dir)
  if(d$protocols[[1]]$phase!="evaluation" || !setequal(unique(d$metrics$scenario),names(wv_scenarios())))
    stop("Final report requires all evaluation scenarios")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  measures<-setdiff(names(d$metrics),c("scenario","replicate","seed","source","id","method","lambda","spatial"))
  summary<-wv_summarize(d$metrics,c("scenario","source","id","method","lambda","spatial"),measures)
  paired<-wv_paired(d$metrics)
  utils::write.csv(summary,file.path(out_dir,"summary.csv"),row.names=FALSE)
  utils::write.csv(paired,file.path(out_dir,"paired_vs_fixed.csv"),row.names=FALSE)
  utils::write.csv(wv_paired(d$metrics,converged_only=TRUE),file.path(out_dir,"paired_converged_vs_fixed.csv"),row.names=FALSE)
  utils::write.csv(wv_paired(d$metrics,baseline="diff2"),file.path(out_dir,"paired_vs_tuned.csv"),row.names=FALSE)
  utils::write.csv(wv_paired(d$metrics,baseline="diff2",converged_only=TRUE),
    file.path(out_dir,"paired_converged_vs_tuned.csv"),row.names=FALSE)
  primary<-reshape(summary[c("scenario","source","method","lambda","metric","mean")],
    idvar=c("scenario","source","method","lambda"),timevar="metric",direction="wide")
  names(primary)<-sub("^mean[.]","",names(primary))
  utils::write.csv(primary,file.path(out_dir,"primary.csv"),row.names=FALSE)
  saveRDS(d$metrics,file.path(out_dir,"replicate_metrics.rds"),compress="xz")
  saveRDS(d$leakage,file.path(out_dir,"leakage.rds"),compress="xz")
  saveRDS(d$protocols,file.path(out_dir,"settings.rds"))
  # Aggregate trials and runs WITHIN each independent dataset before uncertainty.
  fields<-c("scenario","replicate","seed","source","id","gap","null")
  leakage<-stats::aggregate(d$leakage[c("excess_leakage","previous_activity","recovered_previous_activity","target_error")],
    d$leakage[fields],mean)
  utils::write.csv(wv_summarize(leakage,c("scenario","source","id","gap","null"),
    c("excess_leakage","previous_activity","recovered_previous_activity","target_error")),
    file.path(out_dir,"leakage_by_spacing.csv"),row.names=FALSE)
  solver<-d$diagnostics
  utils::write.csv(stats::aggregate(cbind(converged,iterations)~scenario+source+id,solver,
    function(x)c(mean=mean(x),max=max(x))),file.path(out_dir,"convergence.csv"),row.names=FALSE)
  groups<-split(seq_len(nrow(solver)),wv_key(solver,c("scenario","source","id")))
  failures<-do.call(rbind,lapply(groups,function(rows) {
    x<-solver[rows,]
    data.frame(x[1,c("scenario","source","id")],records=nrow(x),failed=sum(!x$converged),
      maximum_iterations=max(x$iterations),maximum_primal_ratio=max(x$primal/pmax(x$primal_limit,1e-30)),
      maximum_dual_ratio=max(x$dual/pmax(x$dual_limit,1e-30)))
  }))
  utils::write.csv(failures,file.path(out_dir,"solver_failures.csv"),row.names=FALSE)
  saveRDS(solver,file.path(out_dir,"solver_diagnostics.rds"),compress="xz")
  saveRDS(d$hrf,file.path(out_dir,"hrf_diagnostics.rds"),compress="xz")
  configs<-unlist(lapply(d$protocols,function(p)p$scenarios),recursive=FALSE)
  true_rho<-vapply(configs,function(x)x$noise_rho,numeric(1))
  hd<-d$hrf;hd$rho_error<-abs(hd$rho-true_rho[hd$scenario]);hd$fallback<-as.numeric(grepl("^fallback",hd$status))
  hd<-stats::aggregate(hd[c("kernel_nmse","rho","rho_error","fallback")],hd[c("scenario","replicate","seed")],mean)
  utils::write.csv(wv_summarize(hd,"scenario",c("kernel_nmse","rho","rho_error","fallback")),
    file.path(out_dir,"hrf_noise_summary.csv"),row.names=FALSE)
  methods<-c("diff2_fixed","diff2","tv","mixture","tv_spatial")
  labels<-c("L2 .01","L2 tuned","TV","Transient + sustained","TV + spatial")
  colors<-c("#666666","#b35806","#2166ac","#1b7837","#762a83")
  grDevices::png(file.path(out_dir,"waveform_recovery.png"),width=2200,height=1500,res=170)
  graphics::par(mfrow=c(2,3),mar=c(6,4,2.5,1))
  for(metric in c("trial_deviation_r","null_target_excess"))for(scenario in c("brief","smooth_heldout","biphasic_heldout")) {
    x<-summary[summary$scenario==scenario & summary$metric==metric,]
    graphics::plot(c(.7,5.3),range(c(0,x$lower,x$upper)),type="n",xaxt="n",
      xlab="",ylab=if(metric=="trial_deviation_r")"Trial-deviation correlation"else"Excess previous-source RMS",
      main=gsub("_"," ",scenario))
    graphics::axis(1,1:5,labels,las=2,cex.axis=.7)
    for(source in c("oracle","learned")) {
      z<-x[x$source==source,];z<-z[match(methods,z$method),]
      position<-1:5+if(source=="oracle")-.08 else .08
      graphics::segments(position,z$lower,position,z$upper,col=colors)
      graphics::points(position,z$mean,col=colors,pch=if(source=="oracle")1 else 19)
    }
    if(metric=="trial_deviation_r" && scenario=="brief")graphics::legend("bottomleft",
      c("True HRF/rho","Fitted HRF/rho"),pch=c(1,19),bty="n",cex=.65)
  }
  grDevices::dev.off()
  # First replicate is fixed by protocol, never selected for visual appeal.
  shard<-dirname(d$paths[vapply(d$protocols,function(p)"smooth_heldout"%in%names(p$scenarios),logical(1))])
  example<-readRDS(file.path(shard,"example_smooth_heldout.rds"))
  saveRDS(example,file.path(out_dir,"example_smooth_heldout.rds"),compress="xz")
  grDevices::png(file.path(out_dir,"example_waveforms.png"),width=2200,height=1500,res=170)
  graphics::par(mfrow=c(3,2),mar=c(4,4,2,1))
  truth<-example$data$truth[[1]];n<-nrow(truth);rows<-33:(n-33L)
  for(i in 1:3)for(source in c("oracle","learned")) {
    pair<-example$data$pairs[[1]][i,];at<-pair$onset+(-6:16)
    t<-wv_interpolate(truth,at)[,1]
    curves<-vapply(methods,function(method) {
      id<-d$protocols[[1]]$settings$id[d$protocols[[1]]$settings$method==method]
      activity<-example$activity[[paste(source,id,sep="|")]]
      wv_interpolate(wv_align_offset(activity,truth,rows),at)[,1]
    },numeric(length(at)))
    graphics::matplot(-6:16,curves,type="l",lty=1,col=colors,ylim=range(c(t,curves)),
      xlab="Seconds from target",ylab="Activity",main=paste(source,"HRF; gap",pair$gap,"s; null",pair$null))
    graphics::lines(-6:16,t,lwd=2,col="black");graphics::abline(v=c(-pair$gap,0),lty=3,col="gray50")
    if(i==1)graphics::legend("topright",c("Truth",labels),col=c("black",colors),lty=1,cex=.55,bty="n")
  }
  grDevices::dev.off()
  grDevices::png(file.path(out_dir,"null_target_examples.png"),width=2000,height=1300,res=170)
  graphics::par(mfrow=c(2,2),mar=c(4,4,2,1))
  for(scenario in c("brief","neural_overlap")) {
    shard<-dirname(d$paths[vapply(d$protocols,function(p)scenario%in%names(p$scenarios),logical(1))])
    ex<-readRDS(file.path(shard,paste0("example_",scenario,".rds")))
    saveRDS(ex,file.path(out_dir,paste0("example_",scenario,".rds")),compress="xz")
    truth<-ex$data$truth[[1]];rows<-33:(nrow(truth)-33L)
    pair<-ex$data$pairs[[1]][1,];at<-pair$onset+(-6:16)
    for(source in c("oracle","learned")) {
      t<-wv_interpolate(truth,at)[,1]
      curves<-vapply(methods,function(method) {
        id<-d$protocols[[1]]$settings$id[d$protocols[[1]]$settings$method==method]
        activity<-ex$activity[[paste(source,id,sep="|")]]
        wv_interpolate(wv_align_offset(activity,truth,rows),at)[,1]
      },numeric(length(at)))
      graphics::matplot(-6:16,curves,type="l",lty=1,col=colors,ylim=range(c(t,curves)),
        xlab="Seconds from silent target",ylab="Activity",main=paste(scenario,source,"HRF; gap",pair$gap))
      graphics::lines(-6:16,t,lwd=2,col="black");graphics::abline(v=c(-pair$gap,0),lty=3,col="gray50")
      graphics::legend("topright",c("Truth",labels),col=c("black",colors),lty=1,cex=.55,bty="n")
    }
  }
  grDevices::dev.off()
  writeLines(c(paste("Report MD5",unname(tools::md5sum("inst/validation/report_waveform_recovery.R"))),
    unlist(lapply(d$paths,function(p)c(paste("Input",dirname(p)),readLines(file.path(dirname(p),"environment.txt"))))),
    paste("All voxel fits converged:",all(solver$converged))),file.path(out_dir,"provenance.txt"))
  invisible(list(summary=summary,paired=paired))
}

if(sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)!=2L)stop("Supply input and output directories")
  report_waveform_recovery(args[1L],args[2L])
}
