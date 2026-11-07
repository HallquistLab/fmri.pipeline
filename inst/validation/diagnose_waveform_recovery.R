# Finite-amplitude perturbation maps for nonlinear continuous inverses and a
# reproducible CPU benchmark. These do not select tuning parameters.
source(file.path("inst","validation","waveform_recovery.R"),local=TRUE)

wv_plot_probes<-function(data,settings,out_dir) {
  ids<-settings$id
  colors<-c("#b35806","#2166ac","#1b7837","#762a83","#666666")
  grDevices::png(file.path(out_dir,"temporal_perturbation.png"),width=2000,height=1300,res=170)
  graphics::par(mfrow=c(2,2),mar=c(4,4,2,1))
  for(context in c("quiet","noisy"))for(amplitude in c(.25,1)) {
    x<-data[data$context==context & data$amplitude==amplitude & data$source_lag==0,]
    graphics::plot(range(x$target_lag),range(c(x$gain,0,1),finite=TRUE),type="n",xlab="Lag from injected impulse (s)",
      ylab="Recovered response / injected amplitude",main=paste(context,"background; amplitude",amplitude))
    graphics::segments(0,0,0,1,lty=3,lwd=2,col="black")
    for(i in seq_along(ids)) {z<-x[x$id==ids[i],];graphics::lines(z$target_lag,z$gain,col=colors[i],lwd=1.5)}
    graphics::abline(h=0,lty=3,col="gray60")
    graphics::legend("topright",ids,col=colors,lty=1,cex=.65,bty="n")
  }
  grDevices::dev.off()
}

wv_plot_benchmark<-function(data,out_dir) {
  ids<-unique(data$id);colors<-c("#b35806","#2166ac","#1b7837","#762a83","#666666")
  grDevices::png(file.path(out_dir,"runtime.png"),width=1900,height=900,res=170)
  graphics::par(mfrow=c(1,2),mar=c(4,4,2,1))
  for(mode in c("cold","warm")) {
    d<-data[data$mode==mode,];z<-stats::aggregate(elapsed~voxels+id,d,median)
    graphics::plot(range(z$voxels),pmax(range(z$elapsed),.001),type="n",log="xy",
      xlab="Voxels in a shared HRF/mask group",ylab="Median seconds",main=paste(mode,"fit, one BLAS thread"))
    for(i in seq_along(ids)) {x<-z[z$id==ids[i],];x<-x[order(x$voxels),];graphics::lines(x$voxels,pmax(x$elapsed,.001),type="b",col=colors[i])}
    graphics::legend("topleft",ids,col=colors,lty=1,pch=1,cex=.65,bty="n")
  }
  grDevices::dev.off()
}

wv_probe_recovery <- function(selection_path,out_dir) {
  selection<-readRDS(selection_path)
  if(!identical(selection$source_md5,tools::md5sum(wv_sources())))stop("Development sources changed")
  wv_load_solver(); dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  settings<-selection$settings; settings$method[settings$method=="diff2_fixed"]<-"diff2"
  set.seed(28000101L); n<-160L; center<-64L
  h<-wv_hrf(0:32,6);h<-h/sum(h)
  design<-wv_design(n,h,1:n,.5); forward<-design$forward[,length(h)-1L+seq_len(n)]
  reference<-matrix(100+wv_ar(n,.5,.04),n,1)
  maps<-list()
  for(i in seq_len(nrow(settings)))for(context in c("quiet","noisy")) {
    setting<-settings[i,,drop=FALSE];cache<-new.env(parent=emptyenv())
    background<-if(context=="quiet")matrix(100,n,1)else reference
    base<-wv_reconstruct(background,matrix(h,ncol=1),.5,rep(TRUE,n),setting,"parcel",cache,
      scale_reference=reference,max_iterations=60000L)
    for(amplitude in c(.25,1))for(lag in 0:8) {
      truth<-matrix(0,n,1);truth[center+lag+1L,1]<-amplitude
      y<-background+forward%*%truth
      answer<-wv_reconstruct(y,matrix(h,ncol=1),.5,rep(TRUE,n),setting,"parcel",cache,
        scale_reference=reference,max_iterations=60000L)
      converged<-all(answer$diagnostic$converged) && all(base$diagnostic$converged)
      # A constant neural offset is unidentified after intercept projection.
      delta<-wv_align_offset(answer$activity-base$activity,truth,33:(n-33L))/amplitude
      maps[[length(maps)+1L]]<-data.frame(id=setting$id,context=context,amplitude=amplitude,
        source_lag=lag,target_lag=-8:20,gain=if(converged)delta[center+(-8:20)+1L,1]else NA_real_,
        ideal=as.numeric(-8:20==lag),converged=converged,
        iterations=max(answer$diagnostic$iterations,base$diagnostic$iterations))
    }
  }
  result<-do.call(rbind,maps)
  utils::write.csv(result,file.path(out_dir,"temporal_perturbation.csv"),row.names=FALSE)
  wv_plot_probes(result,settings,out_dir)
  writeLines(c("Seed: 28000101; TR 1; n=160; true and fitted positive-gamma HRF shape=6; fixed rho=.5.",
    "Finite-amplitude conditional maps: [F(y+H*delta)-F(y)]/amplitude.",
    "BOLD scaling frozen to noisy reference for both contexts and every perturbation.",
    "Probe-only iteration cap 60000; primary study remains 12000; all other stopping settings unchanged.",
    "Nonconverged probe gains are missing and flagged explicitly in the CSV.",
    "Evaluation-only constant offset alignment; no time shift or gain correction.",
    "These are not a global linear resolution matrix or posterior uncertainty.",
    paste("Source MD5",unname(tools::md5sum("inst/validation/diagnose_waveform_recovery.R")))),file.path(out_dir,"perturbation_protocol.txt"))
  invisible(result)
}

wv_benchmark <- function(selection_path,out_dir,sizes=c(8L,64L,256L),repetitions=3L) {
  selection<-readRDS(selection_path)
  if(!identical(selection$source_md5,tools::md5sum(wv_sources())))stop("Development sources changed")
  wv_load_solver();dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  settings<-selection$settings;settings$method[settings$method=="diff2_fixed"]<-"diff2"
  set.seed(28000201L);n<-256L;h<-wv_hrf(0:32,6);h<-h/sum(h)
  mask<-rep(TRUE,n);mask[82:88]<-FALSE
  results<-list()
  for(v in sizes) {
    # Independent voxel noise and signed trial activity; shared kernel/mask.
    neural<-matrix(0,n,v);neural[c(40,70,110,150,190),]<-matrix(rnorm(5*v),5,v)
    f<-wv_design(n,h,which(mask),.5)$forward[,length(h)-1L+seq_len(n)]
    y<-100+f%*%neural+replicate(v,wv_ar(n,.5,.04))
    for(i in seq_len(nrow(settings)))for(rep in seq_len(repetitions)) {
      cache<-new.env(parent=emptyenv());setting<-settings[i,,drop=FALSE]
      for(mode in c("cold","warm")) {
        answer<-wv_reconstruct(y,matrix(rep(h,v),ncol=v),rep(.5,v),mask,setting,rep("parcel",v),cache)
        results[[length(results)+1L]]<-data.frame(id=setting$id,voxels=v,replicate=rep,mode=mode,
          elapsed=answer$elapsed,iterations=max(answer$diagnostic$iterations),convergence=mean(answer$diagnostic$converged))
      }
      utils::write.csv(do.call(rbind,results),file.path(out_dir,"benchmark.csv"),row.names=FALSE)
    }
  }
  writeLines(c(capture.output(sessionInfo()),paste("OPENBLAS_NUM_THREADS",Sys.getenv("OPENBLAS_NUM_THREADS")),
    paste("OMP_NUM_THREADS",Sys.getenv("OMP_NUM_THREADS")),"n=256, one HRF/mask/parcel group, no compilation in timing.",
    "Cold: preprocessing, design/factorization and solve. Warm: preprocessing and solve with cached factors.",
    "Dense validation prototype: not a whole-brain or 50,000-voxel benchmark.",
    paste("Source MD5",unname(tools::md5sum("inst/validation/diagnose_waveform_recovery.R")))),file.path(out_dir,"benchmark_environment.txt"))
  wv_plot_benchmark(do.call(rbind,results),out_dir)
  invisible(do.call(rbind,results))
}

if(sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)!=3L)stop("Supply probe|benchmark, selection.rds, output directory")
  if(args[1L]=="probe")wv_probe_recovery(args[2L],args[3L])else if(args[1L]=="benchmark")
    wv_benchmark(args[2L],args[3L])else stop("Unknown diagnostic")
}
