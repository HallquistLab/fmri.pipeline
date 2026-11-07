wv_solver_env <- new.env(parent = globalenv())
wv_test_env <- function() {
  withr::local_dir(pkg_root)
  if (!exists("wave_solve_cpp", envir = wv_solver_env, mode = "function", inherits = FALSE)) {
    sys.source("inst/validation/run_waveform_recovery.R", wv_solver_env)
    wv_solver_env$wv_load_solver(file.path(tempdir(), "waveform-test-cache"))
  }
  env <- new.env(parent = wv_solver_env)
  sys.source("inst/validation/run_waveform_recovery.R", env)
  env$wave_prepare_cpp <- wv_solver_env$wave_prepare_cpp
  env$wave_solve_cpp <- wv_solver_env$wave_solve_cpp
  env
}

test_that("continuous forward generation retains variable signed shapes and true neural overlap", {
  e<-wv_test_env(); config<-e$wv_scenarios()$spatial_heterogeneity
  config$voxels<-4L; d<-e$wv_generate(config,8531L)
  expect_identical(d,e$wv_generate(config,8531L))
  expect_true(any(d$gain<0)); expect_true(any(d$gain>0))
  for(r in seq_along(d$bold)) {
    expect_equal(d$truth[[r]],d$previous[[r]]+d$target[[r]]+d$background[[r]])
    expect_equal(unname(d$bold[[r]]),100+d$signal[[r]]+d$noise[[r]])
    expect_true(any(d$noise[[r]]!=0))
    expect_equal(nrow(d$bold[[r]]),config$n)
    expect_true(all(c(3,6,10)%in%d$pairs[[r]]$gap))
  }
  b<-e$wv_generate(e$wv_scenarios()$brief,8531L)
  o<-e$wv_generate(e$wv_scenarios()$neural_overlap,8531L)
  at<-b$pairs[[1]]$onset[1]+1
  expect_equal(e$wv_interpolate(b$previous[[1]],at),matrix(0,1,8))
  expect_gt(max(e$wv_interpolate(o$previous[[1]],at)),1)
  expect_true(any(e$wv_shape(seq(0,8,.1),8,"biphasic")<0))
  expect_true(all(e$wv_shape(c(-1,8,9),8,"smooth")==0))
})

test_that("forward convolution agrees with direct numerical integration and preserves causal alignment", {
  e<-wv_test_env(); dt<-.1; t<-seq(0,49,dt); x<-as.numeric(t>=10 & t<12)
  y<-e$wv_forward(x,6)
  direct<-vapply(0:49,function(at)sum(x*e$wv_hrf(at-t,6))*dt,numeric(1))
  expect_equal(y,direct,tolerance=1e-12)
  expect_lt(max(abs(y[1:10])),1e-12)
  expect_gt(max(y),.1)
  expect_equal(sum(e$wv_hrf(0:32,6)/sum(e$wv_hrf(0:32,6))),1)
})

test_that("masked nuisance-adjusted quadratic solver agrees with the package inverse", {
  e<-wv_test_env(); n<-64L; h<-hrf_library(1)[,1L]; rows<-setdiff(1:n,c(17:21,39))
  z<-cbind(1,seq(-1,1,length.out=n)); design<-e$wv_design(n,h,rows,.5,z)
  set.seed(8532); y<-matrix(rnorm(n*3),n,3)
  p<-e$wave_prepare_cpp(design$a,"diff2",e$wv_laplacian(3),0,.01)
  answer<-e$wave_solve_cpp(p,design$transform%*%y[rows,],.01)
  settings<-list(backend="cpp",lambda=.01,penalty="diff2",ridge_floor=1e-8,tail_padding=0L)
  reference<-fmri.pipeline:::hrf_inverse_operator(n,h,rows,.5,z,settings)%*%y[rows,]
  expect_equal(answer$activity[length(h)-1L+seq_len(n),],reference,tolerance=2e-6)
  expect_true(all(answer$converged==1)); expect_equal(answer$iterations,0)
  expect_error(e$wv_design(n,h,c(1,1,4),.5),"is.unsorted")
})

# A deliberately simple R ADMM reference uses explicit analysis matrices and
# fresh linear solves. Native code uses differences, cached inverses and blocks.
wv_reference_admm <- function(a,y,method,lambda,spatial=0,rho=.01,ridge=1e-8) {
  p<-ncol(a); v<-ncol(y); d<-diff(diag(p)); c<-diag(p)
  if(method=="mixture") {
    c<-cbind(diag(p),diag(p)); d<-rbind(cbind(diag(p),matrix(0,p,p)),cbind(matrix(0,p-1,p),d))
  }
  x<-matrix(0,ncol(c),v); z<-d%*%x; u<-z
  g<-crossprod(a%*%c); k<-g+rho*crossprod(d)+diag(ridge,nrow(g))
  eig<-eigen(diag(v)-matrix(1/v,v,v),symmetric=TRUE)
  rhs<-crossprod(a%*%c,y)
  inverses<-lapply(eig$values,function(l)solve(k+spatial*max(0,l)*crossprod(c)))
  for(i in 1:20000) {
    rotated<-(rhs+rho*crossprod(d,z-u))%*%eig$vectors
    x<-sapply(seq_len(v),function(j)inverses[[j]]%*%rotated[,j])%*%t(eig$vectors)
    previous<-z; value<-d%*%x+u; z<-sign(value)*pmax(abs(value)-lambda/rho,0); u<-u+d%*%x-z
    if(max(abs(d%*%x-z))<1e-9 && max(abs(z-previous))<1e-9)break
  }
  c%*%x
}

test_that("native TV, mixed and spatial updates agree with an independent reference", {
  e<-wv_test_env(); set.seed(8533)
  a<-matrix(rnorm(16*10),16,10)/4; y<-matrix(rnorm(16*2),16,2)
  for(method in c("tv","mixture"))for(spatial in c(0,.03)) {
    p<-e$wave_prepare_cpp(a,method,e$wv_laplacian(2),spatial,.03,.5)
    answer<-e$wave_solve_cpp(p,y,.03,20000,1e-8,1e-7)
    reference<-wv_reference_admm(a,y,method,.03,spatial,rho=.5)
    expect_true(all(answer$converged==1))
    expect_equal(answer$activity,reference,tolerance=2e-6)
    expect_true(all(answer$primal<=answer$primal_limit))
    expect_true(all(answer$dual<=answer$dual_limit))
    expect_true(is.finite(answer$objective))
  }
})

test_that("solver exposes failure and rejects invalid observations or spatial priors", {
  e<-wv_test_env(); set.seed(8534); a<-diag(12); y<-matrix(rnorm(24),12,2)
  p<-e$wave_prepare_cpp(a,"tv",e$wv_laplacian(2))
  fail<-e$wave_solve_cpp(p,y,.1,1,1e-12,1e-12)
  expect_false(all(fail$converged==1)); expect_equal(fail$iterations,1)
  expect_error(e$wave_solve_cpp(p,y*NA),"Invalid")
  expect_error(e$wave_prepare_cpp(a,"tv",-diag(2)),"positive semidefinite")
  expect_error(e$wave_prepare_cpp(a,"unknown",diag(2)),"Invalid")
  p_constant<-e$wave_prepare_cpp(a,"tv",e$wv_laplacian(2),admm_rho=1)
  constant<-e$wave_solve_cpp(p_constant,y,100,12000,1e-8,1e-7)
  expect_true(all(constant$converged==1))
  expect_equal(constant$activity,matrix(rep(colMeans(y)/(1+1e-8),each=12),12,2),tolerance=2e-5)
})

test_that("reconstruction is event blind, scale equivariant and preserves the acquisition clock", {
  e<-wv_test_env(); withr::local_dir(pkg_root)
  config<-e$wv_scenarios()$censored_heldout; config$voxels<-4L
  data<-e$wv_generate(config,8535L); k<-e$wv_kernels(data,"oracle")[[1]]
  setting<-e$wv_grid()[6,]; cache<-new.env(parent=emptyenv())
  first<-e$wv_reconstruct(data$bold[[1]],k$h,k$rho,data$masks[[1]],setting,data$groups,cache)
  changed<-data$bold[[1]]; changed[!data$masks[[1]],]<-NaN
  second<-e$wv_reconstruct(changed,k$h,k$rho,data$masks[[1]],setting,data$groups,cache)
  expect_equal(first$activity,second$activity)
  scaled<-e$wv_reconstruct(4*changed,k$h,k$rho,data$masks[[1]],setting,data$groups,cache)
  expect_equal(scaled$activity,4*first$activity,tolerance=1e-7)
  expect_equal(dim(first$activity),c(config$n,config$voxels))
  expect_true(all(is.finite(first$activity)))
  expect_false(any(c("events","rpe","truth")%in%names(formals(e$wv_reconstruct))))
  expect_error(e$wv_interpolate(first$activity,-1),"outside")
})

test_that("leakage metric separates continuing neural activity from reconstruction error", {
  e<-wv_test_env(); data<-e$wv_generate(e$wv_scenarios()$neural_overlap,8536L)
  exact<-list(activity=data$truth[[1]],diagnostic=data.frame(converged=TRUE,iterations=0),elapsed=0)
  counter<-exact; counter$activity<-exact$activity-data$previous[[1]]
  scores<-e$wv_metrics(exact,counter,data,1)
  expect_equal(scores$metrics[["neural_nmse"]],0,tolerance=1e-20)
  expect_equal(scores$metrics[["waveform_r"]],1)
  expect_equal(scores$metrics[["trial_deviation_r"]],1)
  expect_equal(scores$metrics[["trial_deviation_nmse"]],0,tolerance=1e-20)
  expect_lt(scores$metrics[["excess_leakage"]],1e-12)
  expect_gt(scores$metrics[["previous_activity"]],.2)
  expect_equal(scores$metrics[["previous_activity"]],scores$metrics[["recovered_previous_activity"]])
})

test_that("native plan cache handles changing HRF group sizes across runs", {
  e<-wv_test_env();set.seed(8537);n<-64L;h<-hrf_library(1)[,1L]
  y<-100+matrix(rnorm(n*4),n,4);cache<-new.env(parent=emptyenv())
  for(i in c(6L,14L)) {
    setting<-e$wv_grid()[i,]
    first<-e$wv_reconstruct(y,matrix(rep(h,4),ncol=4),rep(.5,4),rep(TRUE,n),setting,rep("parcel",4),cache)
    shared<-e$wv_reconstruct(y[,1:3],matrix(rep(h,3),ncol=3),rep(.5,3),rep(TRUE,n),setting,rep("parcel",3),cache)
    fresh<-e$wv_reconstruct(y[,1:3],matrix(rep(h,3),ncol=3),rep(.5,3),rep(TRUE,n),setting,rep("parcel",3))
    expect_equal(shared$activity,fresh$activity)
    expect_true(all(shared$diagnostic$converged))
    expect_equal(ncol(first$activity),4L)
  }
})

test_that("overlapping event-specific waveforms can be unidentifiable even with a known HRF", {
  e<-wv_test_env(); n<-64; h<-hrf_library(1)[,1]
  first<-second<-numeric(n); first[21:30]<-1; second[26:35]<-.7
  exchanged<-numeric(n); exchanged[26:30]<-.4
  hmat<-e$wv_design(n,h,1:n,0)$forward[,length(h)-1L+seq_len(n)]
  expect_equal(hmat%*%(first+second),hmat%*%((first+exchanged)+(second-exchanged)))
  expect_gt(sum(exchanged^2),0)
})

test_that("development selection excludes failed candidates and protects held-out families", {
  e<-wv_test_env(); withr::local_dir(pkg_root)
  grid<-e$wv_grid(); metrics<-grid; metrics$source<-"oracle"; metrics$neural_nmse<-seq_len(nrow(grid))
  metrics$convergence<-1; metrics$convergence[1]<-.9
  selected<-e$wv_select(metrics)
  expect_equal(selected$id[1],grid$id[2]); expect_equal(nrow(selected),5L)
  expect_true("diff2_fixed/0.010"%in%selected$id)
  tmp<-tempfile();dir.create(tmp);writeLines("frozen",file.path(tmp,"settings.rds"))
  expect_error(e$run_waveform_recovery(tmp,"development",1),"empty output")
  expect_error(e$run_waveform_recovery(tempfile(),"development",1,"smooth_heldout"),"Held-out")
  expect_error(e$run_waveform_recovery(tempfile(),"evaluation",1),"frozen development")
})

test_that("waveform summaries use independent datasets and preserve paired differences", {
  withr::local_dir(pkg_root)
  e<-new.env(parent=globalenv());sys.source("inst/validation/report_waveform_recovery.R",e)
  d<-data.frame(scenario="test",replicate=rep(1:4,2),seed=rep(1:4,2),source="oracle",
    id=rep(c("tv/.03","diff2/.01"),each=4),method=rep(c("tv","diff2_fixed"),each=4),
    lambda=rep(c(.03,.01),each=4),spatial=0,neural_nmse=c(.2,.3,.4,.5,.5,.5,.5,.5))
  p<-e$wv_paired(d)
  expect_equal(p$n,4L);expect_equal(p$mean,-.15)
  expect_equal(p$mc_se,sd(c(-.3,-.2,-.1,0))/2)
  expect_equal(p$lower,-.15-qt(.975,3)*p$mc_se)
  d$convergence<-1;d$convergence[1]<-.9
  filtered<-e$wv_paired(d,converged_only=TRUE)
  expect_equal(filtered$n[filtered$metric=="neural_nmse"],3L)
  expect_equal(filtered$mean[filtered$metric=="neural_nmse"],-.1)
})

test_that("waveform collector rejects missing records, wrong seeds and changed sources", {
  withr::local_dir(pkg_root)
  e<-new.env(parent=globalenv());sys.source("inst/validation/report_waveform_recovery.R",e)
  tmp<-tempfile();dir.create(tmp);config<-e$wv_scenarios()$brief;config$voxels<-2L
  setting<-e$wv_grid()[2,,drop=FALSE]
  protocol<-list(phase="development",repetitions=1L,scenarios=list(brief=config),settings=setting,
    source_md5=tools::md5sum(e$wv_sources()),seed_base=24000000L,sources="oracle")
  saveRDS(protocol,file.path(tmp,"settings.rds"))
  metrics<-data.frame(scenario="brief",replicate=1L,seed=24010001L,source="oracle",id=setting$id,
    method="diff2",lambda=.01,spatial=0,neural_nmse=.2)
  saveRDS(metrics,file.path(tmp,"replicate_metrics.rds"))
  diagnostic<-expand.grid(scenario="brief",replicate=1L,source="oracle",id=setting$id,
    run=1:3,voxel=1:2,fit=c("observed","counterfactual"),stringsAsFactors=FALSE)
  diagnostic$seed<-24010001L
  saveRDS(diagnostic,file.path(tmp,"solver_diagnostics.rds"))
  leakage<-data.frame(scenario="brief",replicate=1L,seed=24010001L,source="oracle",id=setting$id,
    run=rep(1:3,each=5),gap=c(3,6,10,3,6,6,10,3,6,3,10,3,6,3,6),
    null=rep(c(TRUE,FALSE,TRUE,FALSE,TRUE),3),excess_leakage=0)
  saveRDS(leakage,file.path(tmp,"leakage.rds"))
  expect_equal(nrow(e$wv_collect(tmp)$metrics),1L)
  saveRDS(rbind(metrics,metrics),file.path(tmp,"replicate_metrics.rds"))
  expect_error(e$wv_collect(tmp),"Incomplete or duplicate")
  metrics$seed<-99L;saveRDS(metrics,file.path(tmp,"replicate_metrics.rds"))
  expect_error(e$wv_collect(tmp),"Incorrect dataset seeds")
  metrics$seed<-24010001L;saveRDS(metrics,file.path(tmp,"replicate_metrics.rds"))
  saveRDS(diagnostic[-1,],file.path(tmp,"solver_diagnostics.rds"))
  expect_error(e$wv_collect(tmp),"Incomplete solver")
  saveRDS(diagnostic,file.path(tmp,"solver_diagnostics.rds"))
  protocol$source_md5[]<-"changed";saveRDS(protocol,file.path(tmp,"settings.rds"))
  expect_error(e$wv_collect(tmp),"Study sources changed")
})

test_that("resumed serial and parallel datasets reproduce identical recovery metrics", {
  skip_on_os("windows")
  withr::local_dir(pkg_root)
  e<-new.env(parent=wv_test_env());sys.source("inst/validation/resume_waveform_recovery.R",e)
  input<-tempfile();dir.create(input);config<-e$wv_scenarios()$brief
  config$n<-160L;config$voxels<-4L
  protocol<-list(phase="evaluation",repetitions=2L,scenarios=list(brief=config),
    settings=e$wv_grid()[3,,drop=FALSE],source_md5=tools::md5sum(e$wv_sources()),
    seed_base=26000000L,sources=c("oracle","learned"))
  saveRDS(protocol,file.path(input,"settings.rds"))
  serial<-tempfile();multi<-tempfile();e$resume_waveform_recovery(input,"brief",serial)
  e$resume_waveform_recovery(input,"brief",multi,workers=2L)
  a<-readRDS(file.path(serial,"replicate_metrics.rds"));b<-readRDS(file.path(multi,"replicate_metrics.rds"))
  expect_identical(a[setdiff(names(a),"elapsed")],b[setdiff(names(b),"elapsed")])
  expect_identical(readRDS(file.path(serial,"leakage.rds")),readRDS(file.path(multi,"leakage.rds")))
  expect_identical(readRDS(file.path(serial,"solver_diagnostics.rds")),readRDS(file.path(multi,"solver_diagnostics.rds")))
  expect_error(e$resume_waveform_recovery(input,"brief",multi),"empty output")
  expect_error(e$resume_waveform_recovery(input,"unknown",tempfile()),"Unknown evaluation")
})
