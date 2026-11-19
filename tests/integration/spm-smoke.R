# Run from the repository root with MATLAB on PATH:
# SPM_PATH=/path/to/spm12 Rscript tests/integration/spm-smoke.R
# SPM_TEST_OUTPUT optionally preserves outputs at a chosen location.
pkgload::load_all(quiet = TRUE)
source('tests/testthat/helpers-mock_factories.R')
set.seed(739)
root <- Sys.getenv('SPM_TEST_OUTPUT', tempfile('spm-smoke-', tmpdir = Sys.getenv('TMPDIR', '/tmp')))
spm_path <- Sys.getenv('SPM_PATH')
stopifnot(dir.exists(spm_path), nzchar(Sys.which('matlab')))
message('SPM smoke output: ', root)
dir.create(root, recursive = TRUE, showWarnings = FALSE)
files <- vapply(1:2, function(i) {
  dir.create(file.path(root, paste0('run',i)), showWarnings = FALSE)
  p <- file.path(root, paste0('run',i), 'func.nii.gz')
  time <- (0:49 - c(2,1)[i]) * 2
  signal <- vapply(time, function(t) sum(dgamma(t - c(4,16,30,46,66), 6) - dgamma(t - c(4,16,30,46,66), 16)/6), numeric(1))
  data <- array(100 + rep(40*signal, each=125) + rnorm(5*5*5*50, sd=.2), dim=c(5,5,5,50))
  RNifti::writeNifti(data, p)
  p
}, character(1))
reg <- cbind(trial=1:5, onset=c(4,16,30,46,66), duration=2, value=1)
attr(reg,'event') <- 'A'
bdm <- structure(list(design=array(list(reg, reg), dim=c(2,1), dimnames=list(c('run1','run2'),'A')),
  run_niftis=files, run_volumes=c(48,49), drop_volumes=c(2,1), tr=2), class=c('bdm','list'))
conf <- vapply(1:2, function(i) {
  p<-file.path(root,paste0('conf',i,'.txt'))
  write.table(matrix(rnorm(100),50,2),p,row.names=FALSE,col.names=FALSE)
  p
},character(1))
gpa <- create_mock_gpa(n_subjects=1,n_runs=2,include_l1_models=TRUE,output_directory=root)
gpa$glm_software <- 'spm'
gpa$glm_settings <- list(spm=list(spm_path=spm_path, fmri_t=16, fmri_t0=8, cvi='AR(1)', write_residuals=TRUE))
gpa$output_locations$spm_l1_directory <- file.path(root,'{id}','ses-{session}','{l1_model}')
gpa$parallel$compute_environment <- list(global=character(),spm=character(),r=character())
gpa$l1_models$models <- list(model1=structure(list(name='model1',contrasts=matrix(1,1,dimnames=list('A','A'))),class=c('l1_model_spec','list')))
out <- spm_l1_model(id='sub1',session=1L,model_name='model1',d_obj=bdm,run_nifti=files,l1_confound_files=conf,gpa=gpa,execute_spm=TRUE)
stopifnot(isTRUE(out$spm_complete))
manifest <- fmri.pipeline:::spm_read_contrast_manifest(out$spm_dir)
stopifnot(manifest$contrast_name == 'A', manifest$stat == 'T')
# A repeated complete run must recreate cleaned-up temporary inputs.
out2 <- spm_l1_model(id='sub1',session=1L,model_name='model1',d_obj=bdm,run_nifti=files,l1_confound_files=conf,gpa=gpa,execute_spm=TRUE)
stopifnot(isTRUE(out2$spm_complete))
stopifnot(length(list.files(out2$spm_dir, pattern='^Res_.*nii$')) == 97L)

# Exercise the L3 entry point and its validation with a small one-sample model.
l3_data <- data.frame(id=paste0('sub',1:4), session=1L)
l3_model <- fmri.pipeline:::create_new_hi_model(data=l3_data, level=3L, cur_model_names=NULL,
  spec_list=list(name='group',model_formula='~ 1'),lg=lgr::get_logger('spm_smoke'))
gpa$l3_models <- list(models=list(group=l3_model))
gpa$output_locations$spm_l3_directory <- file.path(root, 'group')
gpa$glm_settings$spm$write_residuals <- FALSE
l3_data$con_file <- vapply(1:4, function(i) {
  p <- file.path(root,paste0('subject_con_',i,'.nii'))
  RNifti::writeNifti(array(rnorm(125, mean=2, sd=.3),dim=c(5,5,5)), p)
  p
},character(1))
l3_data$l1_model <- 'model1'; l3_data$l1_cope_name <- 'A'; l3_data$l3_model <- 'group'
l3 <- spm_l3_model(l3_data, gpa, execute_spm=TRUE)
stopifnot(isTRUE(l3$spm_complete))
unlink(manifest$stat_file)
stopifnot(!get_spm_status(out$spm_dir)$spm_complete)
cat('SPM runtime smoke passed: real AR(1) fit, residuals, contrast, repeat after cleanup, L3 fit, missing-stat detection.\n')
