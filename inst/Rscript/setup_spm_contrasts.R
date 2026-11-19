#!/usr/bin/env Rscript

#read in command line arguments.
args <- commandArgs(trailingOnly = FALSE)

scriptpath <- dirname(sub("--file=", "", grep("--file=", args, fixed=TRUE, value=TRUE), fixed=TRUE))
argpos <- grep("--args", args, fixed=TRUE)
if (length(argpos) > 0L) {
  args <- args[(argpos+1):length(args)]
} else {
  args <- c()
}

condition_contrasts <- FALSE
unit_contrasts <- FALSE
effects_of_interest_F <- FALSE
spm_path <- NULL

argpos <- 1
while (argpos <= length(args)) {
  if (args[argpos] == "-mat_file") {
    matfile <- args[argpos + 1] #name of preprocessed fMRI data
    stopifnot(file.exists(matfile))
    argpos <- argpos + 2
  } else if (args[argpos] == "-condition_contrasts") {
    condition_contrasts <- as.logical(args[argpos+1])
    argpos <- argpos + 2
  } else if (args[argpos] == "-unit_contrasts") {
    unit_contrasts <- as.logical(args[argpos+1])
    argpos <- argpos + 2
  } else if (args[argpos] == "-effects_of_interest_F") {
    effects_of_interest_F <- as.logical(args[argpos+1])
    argpos <- argpos + 2
  } else if (args[argpos] == "-spm_path") {
    spm_path <- args[argpos+1]
    stopifnot(dir.exists(spm_path))
    argpos <- argpos + 2
  } else {
    stop("Not sure what to do with argument: ", args[argpos])
  }
}

library(R.matlab)
spm_info <- readMat(matfile)
mnames <- unlist(spm_info$mnames)
cpos <- as.vector(spm_info$cpos)
bpos <- as.vector(spm_info$bpos)
npos <- as.vector(spm_info$npos)

#output script to same location as matfile
output_dir <- normalizePath(dirname(matfile))

spm_preamble <- c(
  ifelse(is.null(spm_path), "", paste0("addpath('", spm_path, "');")),
  "spm('defaults', 'fmri');",
  "spm_jobman('initcfg');",
  ""    
)

#by default, only estimate contrasts of interest
minterest <- mnames[cpos]
reg_names <- sub("Sn\\(\\d+\\)\\s+", "", minterest, perl=TRUE)
reg_split <- split(minterest, reg_names)

#don't estimate the run intercept contrast (should be handled already by using cpos above)
#if (any(which_const <- names(reg_split) == "constant")) { reg_split <- reg_split[!which_const] }

cn <- 1 #contrast number
expected_names <- character(0)

m_string <- c(
  spm_preamble,
  "% SETUP BATCH JOB STRUCTURE",
  "contrast = struct;",
  "% spmmat",
  paste0("contrast.matlabbatch{1}.spm.stats.con.spmmat = { ['", output_dir, "' filesep 'SPM.mat']};"),
  ""
)

if (effects_of_interest_F) {
  #f contrast tutorial: https://www.youtube.com/watch?v=aqWtk4-ayg4
  cvec <- rep(0, length(mnames))
  cvec[cpos] <- 1 #unit height for all effects of interest
  m_string <- c(m_string,
    paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn, "}.fcon.name = 'Effects of interest';"),
    paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn, "}.fcon.weights = [", paste(cvec, collapse=", "), "];"),
    paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn, "}.fcon.sessrep = 'none';")
  )
  expected_names <- c(expected_names, "Effects of interest")
  cn <- cn + 1 #increment contrast number
}

if (condition_contrasts) {
  # compute average contrasts for conditions (potentially across multiple sessions/runs)

  cmat <- c()
  for (rr in 1:length(reg_split)) {
    nruns <- length(reg_split[[rr]])
    cweight <- 1/nruns
    cvec <- rep(0, length(mnames))
    cvec[ mnames %in% reg_split[[rr]] ] <- cweight
    cmat <- rbind(cmat, cvec)
  }
  rownames(cmat) <- paste0(make.names(names(reg_split)), ".ravg")
  colnames(cmat) <- make.names(mnames)

  for (ii in 1:nrow(cmat)) {
    m_string <- c(m_string,
      paste0("% consess ", cn + ii - 1),
      paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn + ii - 1, "}.tcon.name = '", rownames(cmat)[ii], "';"),
      paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn + ii - 1, "}.tcon.convec = [ ", paste(cmat[ii,], collapse=", "), " ];"),
      paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn + ii - 1, "}.tcon.sessrep = 'none';") #not using session replication at the moment...
    )
  }
  
  expected_names <- c(expected_names, rownames(cmat))
  cn <- cn + nrow(cmat) #increment contrast number

  m_string <- c(m_string,
    "",
    "% delete?",
    "contrast.matlabbatch{1}.spm.stats.con.delete = 1;",
    "% RUN BATCH JOB",
    "spm_jobman('run',contrast.matlabbatch);"
  )

  cat(m_string, file=file.path(output_dir, "estimate_glm_contrasts.m"), sep="\n")
  writeLines(expected_names, file.path(output_dir, "spm_expected_contrasts.txt"))

}
