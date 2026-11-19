#' R6 class for 3dClustSim automation
#' @importFrom tibble tibble
#' @importFrom stats qnorm qt
#' @export
afni_3dclustsim <- R6::R6Class("afni_3dclustsim",
  private = list(
    out_dir = getwd(),
    clustsim_df = tibble::tibble(),
    out_files = NULL,
    prefix = "3dclustsim",
    pvt_use_fwhmx_acf = FALSE,
    acf_params = setNames(rep(NA_real_, 3), c("a", "b", "c")),
    clustsim_batch = NULL,
    clustsim_call = NULL,
    clustsim_mask = NULL,
    clustsim_complete = FALSE,

    # controls the volume over which the 3dClustSim simulation occurs (see AFNI help)
    # 'mask': Uses the -mask dataset to specify the spatial extent of the simulation
    # 'inset': use 3dttest++ null datasets in BRIK/HEAD form as input to 3dClustsim
    # 'insdat': 3dttest++ null datasets in sdat form as input to 3dClustsim
    # 'xyz': Use the spatial extent specified by the xyz matrix size and the voxel size
    volume_method = NULL,
    nopad = "",
    pthr = ".01 .005 .002 .001 .0005 .0002 .0001",
    athr = ".05 .02 .01 .005 .002 .001 .0005 .0002 .0001",
    iter = 30000L,
    nodec = "",
    seed = 0,
    sim_string = "",
    scheduler = "slurm",
    ncpus = 8L,
    residuals_file = NULL,
    residuals_mask_file = NULL,
    residuals_njobs = 32, # split iter into this many jobs for permutations
    which_sim = function() {
      if (private$volume_method == "inset") {
        private$sim_string <- glue::glue("-inset {paste(self$inset_files, collapse=' ')}")
        if (!is.null(private$clustsim_mask)) {
          private$sim_string <- paste(private$sim_string, glue::glue("-mask {private$clustsim_mask}")) # -mask compatible with -inset
        }
      } else if (private$volume_method == "insdat") {
        private$sim_string <- glue::glue("-insdat {self$insdat_mask_file} {self$insdat_file}")

        # note that the clustsim mask controls the bounds of the simulation for FWE, 
        # while the insdat_mask is for identifying values in the insdat data file
        if (!is.null(private$clustsim_mask)) {
          private$sim_string <- paste(private$sim_string, glue::glue("-mask {private$clustsim_mask}")) # -mask compatible with -inset
        }
      } else if (private$volume_method == "mask") {
        private$sim_string <- glue::glue("-mask {private$clustsim_mask}")
      } else if (private$volume_method == "xyz") {
        # use voxel size and matrix size
        private$sim_string <- glue::glue("-nxyz {paste(private$nxyz, collapse=' ')} -dxyz {paste(private$dxyz, collapse=' ')}")
      } else {
        stop("Cannot sort out how to setup 3dClustSim volume based on method:", private$volume_method)
      }
    },
    build_call = function() {
      if (private$volume_method %in% c("inset", "insdat")) {
        acf_string <- "" # irrelevant
      } else if (isTRUE(private$pvt_use_fwhmx_acf)) {
        if (isFALSE(self$fwhmx_set$is_complete())) {
          warning("Cannot build 3dClustSim call because some 3dFWHMx runs are incomplete.")
          return(invisible(NULL))
        }
        acf_string <- glue::glue("-acf {paste(self$fwhmx_set$get_acf_average(), collapse=' ')}")
      } else { # assuming user passed acf params in manually
        acf_string <- glue::glue("-acf {paste(private$acf_params, collapse=' ')}")
      }

      # setup the simulation volume string specification
      private$which_sim()
      private$clustsim_call <- glue::glue(
        "3dClustSim {private$sim_string} -iter {private$iter} {acf_string} -prefix {private$prefix}",
        " -pthr {private$pthr} -athr {private$athr} -seed {private$seed}{private$nopad}{private$nodec}"
      )
    },
    build_df = function() {
      if (isTRUE(private$clustsim_complete) && nrow(private$clustsim_df) > 0L) {
        return(invisible(NULL)) # data frame has already been compiled
      }

      clustsim_files <- self$get_clustsim_output_files()
      if (length(clustsim_files) == 0L) {
        # most likely, the $submit method has not completed yet.
        return(invisible(NULL))
      } else if (length(clustsim_files) == 9L) {
        private$clustsim_complete <- TRUE
      } else {
        warning("Fewer than 9 clustsim files output files detected: ", paste(clustsim_files, collapse=", "))
      }

      read_clustsim_1d <- function(file) {
        checkmate::assert_file_exists(file)
        nn <- as.integer(sub("^.*NN([1-3])_.*", "\\1", file, perl = TRUE))
        sided <- sub("^.*NN[1-3]_(1|2|bi)sided.*", "\\1", file, perl = TRUE)
        data <- read.table(file, header = FALSE)
        header <- grep(pattern = "^\\s*#.*", readLines(file), value = TRUE)
        p_line <- gsub("[#|]", "", grep(pattern = "^\\s*#\\s*pthr.*", header, value = TRUE)) # strip # and | from header line
        c_names <- strsplit(trimws(p_line), "\\s+")[[1]]
        stopifnot(length(c_names) == ncol(data))
        names(data) <- c_names
        data_long <- data %>%
          tidyr::pivot_longer(cols = !pthr, names_to = "athr", values_to = "nvoxels") %>%
          dplyr::mutate(athr = as.numeric(athr), nn = !!nn, sided = !!sided) %>%
          dplyr::select("nn", "sided", "pthr", "athr", "nvoxels")
        return(data_long)
      }

      # read each clustsim output and compile into single indexed data.frame
      private$clustsim_df <- dplyr::bind_rows(lapply(clustsim_files, read_clustsim_1d))
    }
  ),
  public = list(
    #' @field fwhmx_set A afni_3dfwhmx_list object containing 3dFWHMx information for all fwhmx_input_files
    fwhmx_set = NULL,

    #' @field inset_files A character vector of files to use directly as volumes to threshold and clusterize
    inset_files = NULL,

    #' @field insdat_file An sdat file containing permutations to be passed to 3dClustSim through -insdat
    insdat_file = NULL,

    #' @field insdat_mask_file A mask file corresponding to the insdat_file data that indicates where each value is in space
    insdat_mask_file = NULL,

    #' @field null_3dttest_obj only used if residuals_file is passed in, this contains the object for running the permutations
    null_3dttest_obj = NULL,

    #' @description Create a new afni_3dclustsim object
    #' @param out_dir the intended output directory for 3dClustSim files
    #' @param prefix the prefix to be included in the names of 3dClustSim output files
    #' @param fwhmx_input_files A character vector of input files to be passed through 3dFWHMx (-ACF method)
    #' @param fwhmx_mask_files A character vector of masks containing the volume over which to estimate the ACF in 3dFWHMx. 
    #'   Must match 1:1 with \code{fwhmx_input_files}.
    #' @param residuals_file The filename of the group residuals file to be used in null dataset generation (permutation approach).
    #' @param residuals_mask_file The volume over which null datasets should be generated from the residuals
    #' @param residuals_njobs The number of independent jobs for splitting up the residuals permutations. If NULL, 32 jobs will be used.
    #' @param inset_files A character vector of dataset files to pass directly to 3dClustSim via -inset.
    #' @param insdat_file An sdat file containing null datasets to be passed to 3dClustSim via -insdat.
    #' @param insdat_mask_file A mask file corresponding to the insdat_file that maps values to space.
    #' @param dxyz the size of voxels in x y z (vector of 3 numbers)
    #' @param nxyz the number of voxels along x y z (vector of 3 positive integers)
    #' @param clustsim_mask This controls the volume over which to correct for FWE using 3dClustSim. If you give a whole-brain mask,
    #'   then your cluster thresholds reflect whole-brain FWE correction. If you give a smaller mask (e.g., a single region or network),
    #'   you are correcting only over that volume (i.e., a small-volume correction).
    #' @param acf_params a vector of 3 autocorrelation parameters (a, b, c) to be used to simulate smoothness. Usually produced by 3dFWHMx. 
    #'   The 'a' parameter must be between 0 and 1. The 'b' and 'c' parameters (scale radii) must be positive. The spatial autocorrelation function 
    #'   is given by: `ACF(r) = a * exp(-r*r/(2*b*b)) + (1-a)*exp(-r/c)`
    #' @param nopad If TRUE, disable 3dClustSim's default 'padding' slices along each face to allow for edge effects of the smoothing process. Default: FALSE
    #' @param pthr A vector of voxelwise p-values to be tested in 3dClustSim (-pthr). Can be a space-separated string or
    #'   a numeric vector. Default: ".01 .005 .002 .001 .0005 .0002 .0001".
    #' @param athr A vector of cluster p-values to be tested in 3dClustSim (-athr). Can be a space-separated string or
    #'   a numeric vector. Default: ".05 .02 .01 .005 .002 .001 .0005 .0002 .0001".
    #' @param iter The number of iterations to use in simulating null datasets in 3dClustSim. Default: 30000.
    #' @param nodec If TRUE, clusters will be printed without decimal places, rounding up (e.g., 27.2 becomes 28). Default: FALSE.
    #' @param seed The seed to use when starting 3dClustSim random number generation. Default: 0 (sets a random seed)
    #' @param scheduler The HPC scheduler to use. Can be 'local', 'slurm', or 'torque'.
    #' @param ncpus The number of cores to use in the 3dClustSim job. This sets OMP_NUM_THREADS in the 3dClustSim job to 
    #'   speed up computation.
    initialize = function(out_dir = NULL, prefix=NULL,
                          fwhmx_input_files = NULL, fwhmx_mask_files = NULL,
                          residuals_file = NULL, residuals_mask_file = NULL, residuals_njobs = NULL,
                          inset_files = NULL, insdat_file = NULL, insdat_mask_file = NULL,
                          dxyz = NULL, nxyz = NULL,
                          clustsim_mask = NULL, acf_params = NULL, nopad = NULL, pthr = NULL, athr = NULL, iter = NULL, nodec = NULL, 
                          seed = NULL, scheduler = NULL, ncpus = NULL) {

      # only one method can be used for the volume over which 3dClustSim simulates. (the first two are parts of the same method)
      n_passed <- sapply(list(dxyz, nxyz, clustsim_mask, residuals_file, inset_files, insdat_file, insdat_mask_file), function(x) !is.null(x))
      methods_passed <- 0
      if (!is.null(dxyz) || !is.null(nxyz)) methods_passed <- methods_passed + 1
      if (!is.null(inset_files)) methods_passed <- methods_passed + 1
      if (!is.null(insdat_file) || !is.null(insdat_mask_file)) methods_passed <- methods_passed + 1
      if (!is.null(clustsim_mask) && is.null(insdat_file) && is.null(inset_files)) methods_passed <- methods_passed + 1

      if (methods_passed != 1L) {
        stop("Only one method can be passed for the 3dClustSim volume settings. Use dxyz + nxyz, clustsim_mask, residuals_file, insdat_file, or inset_files")
      }

      if (!is.null(scheduler)) {
        checkmate::assert_string(scheduler)
        checkmate::assert_subset(scheduler, c("torque", "slurm", "local"))
        private$scheduler <- scheduler
      }

      if (!is.null(out_dir)) {
        checkmate::assert_directory_exists(out_dir)
        private$out_dir <- out_dir
      }

      if (!is.null(prefix)) {
        checkmate::assert_string(prefix)
        private$prefix <- prefix
      }

      if (!is.null(fwhmx_input_files)) {
        private$pvt_use_fwhmx_acf <- TRUE
        self$fwhmx_set <- afni_3dfwhmx_list$new(input_files = fwhmx_input_files, mask_files = fwhmx_mask_files, scheduler=scheduler)
        if (!is.null(acf_params)) {
          stop("Cannot pass fwhmx_input_files and acf_params since the ACF parameters are calculated by 3dFWHMx!")
        }
      }

      if (!is.null(inset_files)) {
        if (!is.null(fwhmx_input_files)) {
          stop("Cannot specify both inset_files and fwhmx_input_files. See 3dClustSim documentation!")
        }

        private$volume_method <- "inset"
        checkmate::assert_file_exists(inset_files)
        self$inset_files <- inset_files
      }

      if (!is.null(insdat_file)) {
        private$volume_method <- "insdat"
        checkmate::assert_string(insdat_file)
        checkmate::assert_file_exists(insdat_file)

        checkmate::assert_string(insdat_mask_file)
        checkmate::assert_file_exists(insdat_mask_file)
        self$insdat_file <- insdat_file
        self$insdat_mask_file <- insdat_mask_file
      }

      # if a residuals file is provided, we pass it to 3dttest++ for permutation testing
      if (!is.null(residuals_file)) {
        checkmate::assert_file_exists(residuals_file)
        private$residuals_file <- normalizePath(residuals_file)
        private$volume_method <- "residuals"
        if (is.null(residuals_mask_file)) {
          stop("At present, a residuals_mask_file must be passed with the residuals_file.")
        } else {
          checkmate::assert_file_exists(residuals_mask_file)
          private$residuals_mask_file <- residuals_mask_file
        }

        if (!is.null(residuals_njobs)) {
          checkmate::assert_integerish(residuals_njobs, lower = 1, upper = 1e5)
          private$residuals_njobs <- as.integer(residuals_njobs)
        }

        # if out_dir is not provided, default to same directory as residuals file
        if (is.null(out_dir)) {
          private$out_dir <- dirname(private$residuals_file)
        }

        # note that the number of permutations should always equal the number of iterations of 3dClustSim (AFNI resets this internally anyhow)
        self$null_3dttest_obj <- simulate_null_3dttest$new(
          residuals_file = residuals_file, mask_file = residuals_mask_file,
          n_permutations = private$iter, njobs = private$residuals_njobs
        )
      }

      if (!is.null(dxyz) || !is.null(nxyz)) {
        private$volume_method <- "xyz"
        if (!is.null(dxyz)) {
          private$dxyz <- check_nums(dxyz, lower = .01, upper = 100) # millimeters
        } else {
          stop("Spatial domain 3dClustSim specification requires both dxyz and nxyz. I only received nxyz.")
        }

        if (!is.null(nxyz)) {
          private$nxyz <- check_nums(nxyz, lower = 1, upper = 10000) # voxels
        } else {
          stop("Spatial domain 3dClustSim specification requires both dxyz and nxyz. I only received dxyz.")
        }
      }

      if (!is.null(clustsim_mask)) {
        checkmate::assert_string(clustsim_mask)
        checkmate::assert_file_exists(clustsim_mask)
        private$clustsim_mask <- clustsim_mask

        # if inputs have not triggered a method setting of inset, insdat, residuals, or xyz above, we must be
        # using the specified mask as the volume over which to simulate.
        if (is.null(private$volume_method)) {
          private$volume_method <- "mask"
        }
      } else {
        message("No clustsim_mask argument received. All voxels in volume will be used for cluster simulations!")
      }

      # user specification of ACF params (e.g., if 3dFWHMx done externally)
      if (!is.null(acf_params)) {
        checkmate::assert_numeric(acf_params, len=3)
      }

      if (!is.null(nopad)) {
        checkmate::assert_logical(nopad, len = 1L)
        if (isTRUE(nopad)) private$nopad <- " -nopad"
      }

      if (!is.null(pthr)) {
        private$pthr <- check_nums(pthr, lower = 1e-10, upper = .999)
      }

      if (!is.null(athr)) {
        private$athr <- check_nums(athr, lower = 1e-10, upper = .999)
      }

      if (!is.null(iter)) {
        checkmate::assert_integerish(iter, lower = 10, upper = 1e8)
        private$iter <- as.integer(iter)
      }

      if (!is.null(nodec)) {
        checkmate::assert_logical(nodec, len = 1L)
        if (isTRUE(nodec)) private$nodec <- " -nodec"
      }

      if (!is.null(seed)) {
        checkmate::assert_integerish(seed, len = 1L, lower=0)
        private$seed <- as.integer(seed)
      }

      if (!is.null(ncpus)) {
        checkmate::assert_integerish(ncpus, lower = 1, upper = 1e3, len = 1)
        private$ncpus <- ncpus
      }

      # always put the clustsim method in the filename for clarity
      private$prefix <- paste(private$prefix, private$volume_method, sep = "_")

      # populate 3dClustSim object at initialization (if output files are cached already)
      private$build_df()

      # always refresh the object at initialization so that if the underlying 3dFWHMx or 3dttest++ files are complete,
      # these are current in the created clustsim object.
      self$refresh()

    },

    #' @description submit the 3dClustSim compute job to the cluster
    #' @param force if TRUE, re-estimate 3dClustSim even though its output files already exist
    #' @details Note that if fwhmx_input_files are provided at the corresponding 3dFWHMx has not been run yet,
    #'   this job will be submitted and serve as a dependency for the 3dClustsim job. Likewise, if a residuals_file is
    #'   provided, then the 3dttest++ permutation step will be submitted as a job and 3dClustSim will be dependent on this finishing.
    submit = function(force=FALSE) {
      if (isTRUE(private$clustsim_complete) && isFALSE(force)) {
        message("3dClustSim already finished for these inputs. Use $submit(force=TRUE) to re-run or $get_clustsim_df() to retrieve results.")
        return(invisible(NULL))
      }

      # create batch job
      private$clustsim_batch <- R_batch_job$new(
        job_name = "run_3dclustsim", n_cpus = private$ncpus,
        wall_time = private$wall_time, scheduler = private$scheduler,
        input_objects = list(csim_obj = self), # pass the current object as input to the batch
        r_packages = "fmri.pipeline",
        r_code = c(
          "# refresh completeness of 3dFWHMx runs in case these were run by a preceding batch job",
          "csim_obj$refresh()", # any 3dFWHMx calculations or 3dttest++ permutations are now populated
          "setwd(csim_obj$get_out_dir())",
          "run_afni_command(csim_obj$get_call(), omp_num_threads = csim_obj$get_ncpus())"
        ),
        batch_directory = private$out_dir
      )

      # need to run 3dFWHMx before 3dClustSim can run
      if (isTRUE(private$pvt_use_fwhmx_acf) && isFALSE(self$fwhmx_set$is_complete())) {
        fwhmx_batch <- self$fwhmx_set$get_batch()
        batch_obj <- private$clustsim_batch
        batch_obj$depends_on_parents <- "run_3dfwhmx"
        clust_seq <- R_batch_sequence$new(fwhmx_batch, batch_obj)
        clust_seq$submit()

      # the volume method of 'residuals' specifies that we need to run 3dtest++ first before 3dClustSim
      } else if (isTRUE(private$volume_method == "residuals") && isFALSE(self$null_3dttest_obj$is_complete())) {
        null_batch <- self$null_3dttest_obj$get_batch()
        batch_obj <- private$clustsim_batch
        batch_obj$depends_on_parents <- "run_3dttest"
        clust_seq <- R_batch_sequence$new(null_batch, batch_obj)
        clust_seq$submit()

      # only 3dClustSim needs to run
      } else {
        private$clustsim_batch$submit()
      }

      return(invisible(self))
    },

    #' @description return a data.frame containing the results of the 3dClustSim for all NN levels and sidedness.
    get_clustsim_df = function() {
      private$build_df()
      private$clustsim_df
    },

    #' @description return a character vector of output text files generated by 3dClustSim
    get_clustsim_output_files = function() {
      list.files(path = private$out_dir, pattern = glue::glue("{private$prefix}.*sided\\.1D"), full.names = TRUE)
    },

    #' @description return the 3dClustSim call related to these inputs. This is what will be passed to the scheduler.
    get_call = function() {
      private$build_call()
      private$clustsim_call
    },

    #' @description return the output directory for 3dClustSim files
    get_out_dir = function() {
      private$out_dir
    },

    #' @description return the number of cpus (cores) to be used by 3dClustsim (via OMP_NUM_THREADS)
    get_ncpus = function() {
      private$ncpus
    },

    #' @description return TRUE/FALSE for whether this clustsim relies on the ACF estimates from 3dFWHMx
    use_fwhmx_acf = function() {
      private$pvt_use_fwhmx_acf
    },

    #' @description return whether 3dClustSim has completed for this input
    is_complete = function() {
      self$refresh() # always refresh object status in case it completed through a scheduled job
      private$clustsim_complete
    },

    #' @description Simple method to refresh the clustsim_df, fwhmx files, and 3dttest++ permutation files
    #' @details The refresh is useful if the object needs to be updated just in time to determine whether permutations 
    #'   or ACF params are available.
    refresh = function() {
      # read clustsim output files, if available
      private$build_df()
      if (isTRUE(self$use_fwhmx_acf())) {
        self$fwhmx_set$refresh()
      }

      # if 3dttest++ permutations completed, switch method from residuals to insdat and populate files
      if (private$volume_method == "residuals" && isTRUE(self$null_3dttest_obj$is_complete())) {
        if (isTRUE(self$null_3dttest_obj$get_use_sdat())) {
          private$volume_method <- "insdat"
          ofiles <- self$null_3dttest_obj$get_permutation_files()
          self$insdat_file <- ofiles["permutation_file"]
          self$insdat_mask_file <- ofiles["mask_file"]
        } else {
          stop("Not supported yet...")
          # need to switch to -inset files
        }
      }
    },

    #' @description method to apply 3dClustSim results to a statistic image to create an integer-valued cluster mask and/or a thresholded
    #'   statistic image
    #' @param statistic_nifti A filename to a NIfTI image containing voxelwise statistics that should be used to calculate the threshold for pthr
    #' @param NN The cluster definition from 3dClustSim to be applied when thresholding the image (1, 2, or 3)
    #' @param sided Whether to apply the cluster threshold for one-sided, two-sided, or bi-sided tests ('1', '2', or 'bi')
    #' @param athr The clusterwise p-value to be applied. Default: .05
    #' @param pthr The voxelwise threshold to be applied to the statistic_file. Default: .001.
    #' @param voxelwise_stat A list object specifying the statistic contained in the \code{statistic_nifti}.
    #'   At present, this consists of the \code{stat_type}: 'z', 't', or 'p'. If a 't' statistic is passed,
    #'   also include $df in the list specifying the degrees of freedom.
    #' @param output_cluster_mask If \code{TRUE}, output an integer-valued mask containing any whole brain-significant clusters
    #'   (calculated by 3dClusterize).
    #' @param output_thresholded_image If \code{TRUE}, output a copy of \code{statistic_nifti} that has been thresholded
    #'   according to the settings specified here.
    #' @param add_whereami Whether to also call AFNI whereami to get anatomical landmarks of interest. Default: TRUE
    #' @param whereami_atlases An optional character vector of atlases to be requested in whereami.
    #' @return an afni_3dclusterize object containing clusters in \code{statistic_nifti} that survive the 3dClustSim
    #'   correction requested.
    apply_clustsim = function(statistic_nifti = NULL, NN = 1, sided = "bi", athr = .05, pthr = .001,
      voxelwise_stat = list(stat_type = "z"), output_cluster_mask = TRUE, output_thresholded_image = FALSE,
      add_whereami = TRUE, whereami_atlases = NULL) {
      # eventually, would be nice to allow for multiple rows to be tolerated in sim_calc and to apply each to the data
      checkmate::test_file_exists(statistic_nifti)
      statistic_nifti <- normalizePath(statistic_nifti) #make sure it's a clear absolute path
      checkmate::assert_integerish(NN, len = 1L, lower = 1, upper = 3)

      bisided <- onesided <- twosided <- FALSE

      checkmate::assert_string(sided)
      sided <- dplyr::recode(tolower(sided), "one" = "1", "two" = "2") # consistent nomenclature
      checkmate::assert_subset(sided, c("1", "2", "bi"))

      if (sided == "bi") {
        bisided <- TRUE
      } else if (sided == "1") {
        onesided <- TRUE
      } else if (sided == "2") {
        twosided <- TRUE
      }

      # for bisided and twosided, divide threshold alpha by 2 (both tails)
      test_p <- ifelse(bisided || twosided, pthr / 2, pthr)

      if (!self$is_complete()) {
        stop("Cannot use $apply_clustsim() method because 3dClustSim is not complete. Need to use $submit() first and for this to finish!")
      }

      df <- self$get_clustsim_df()
      if (is.null(df)) {
        warning("Cannot apply 3dClustSim results because $get_clustsim_df() returns NULL. Has the $submit() method been called yet and did it complete successfully?")
        return(invisible(NULL))
      }

      if (!any(abs(athr - df$athr) <= abs(athr) * 1e-10)) { stop("Specified athr: ", athr, " not found in 3dClustSim results") }
      if (!any(abs(pthr - df$pthr) <= abs(pthr) * 1e-10)) { stop("Specified pthr: ", pthr, " not found in 3dClustSim results") }

      sim_calc <- df %>% dplyr::filter(nn == !!NN & sided == !!sided &
        abs(.data$athr - !!athr) <= abs(!!athr) * 1e-10 &
        abs(.data$pthr - !!pthr) <= abs(!!pthr) * 1e-10)
      if (nrow(sim_calc) == 0L) {
        stop("could not find")
      } else if (nrow(sim_calc) > 1L) {
        print(sim_calc)
        stop("More than one threshold found for this combination of settings")
      }

      # calculate the threshold value based on the statistic type
      stype <- tolower(as.character(voxelwise_stat$stat_type)[1L])
      if (stype == "z") {
        thresh_val <- stats::qnorm(test_p, lower.tail = FALSE)
      } else if (stype == "p") {
        thresh_val <- test_p
      } else if (stype == "t") {
        thresh_val <- stats::qt(test_p, df = voxelwise_stat$df, lower.tail = FALSE)
      } else if (stype %in% c("chisq", "chi2")) {
        thresh_val <- stats::qchisq(pthr, df = voxelwise_stat$df, lower.tail = FALSE)
      } else if (stype == "f") {
        df1_val <- voxelwise_stat$df1 %||% voxelwise_stat$df[1L]
        df2_val <- voxelwise_stat$df2 %||% voxelwise_stat$df[2L]
        thresh_val <- stats::qf(pthr, df1 = df1_val, df2 = df2_val, lower.tail = FALSE)
      } else {
        stop("Cannot interpret stat_type ", voxelwise_stat$stat_type)
      }

      if (bisided || twosided) {
        lower_thresh <- -1 * thresh_val
        upper_thresh <- thresh_val

        # not sure this really makes sense... like do we really want the most non-significant p-values?
        if (stype == "p") lower_thresh <- 1 + lower_thresh
      } else {
        one_thresh <- thresh_val
      }

      clusters_file <- NULL
      if (isTRUE(output_cluster_mask)) {
        clusters_file <- paste0(
          file_sans_ext(statistic_nifti),
          glue::glue("_3dC_clustmask_voxp{pthr}_clusp{athr}_NN{NN}_{sided}sided.nii.gz")
        )
      }

      thresholded_stat_file <- NULL
      if (isTRUE(output_thresholded_image)) {
        thresholded_stat_file <- paste0(
          file_sans_ext(statistic_nifti),
          glue::glue("_3dC_thresholded_voxp{pthr}_clusp{athr}_NN{NN}_{sided}sided.nii.gz")
        )
      }

      clust_nvox <- sim_calc %>%
        pull("nvoxels") %>%
        ceiling(.)

      arg_list <- list(
          threshold_file = statistic_nifti, bisided = bisided, onesided = onesided, twosided = twosided,
          NN = NN, clust_nvox = clust_nvox, pref_map = clusters_file, pref_dat = thresholded_stat_file
      )
      simulation_mask <- private$clustsim_mask %||% private$residuals_mask_file %||% self$insdat_mask_file
      if (!is.null(simulation_mask)) {
        ptfce_read_mask(simulation_mask, statistic_nifti)
        arg_list$mask <- simulation_mask
      }

      cat(glue("{clust_nvox} voxels needed for a whole-brain FWE cluster at voxelwise p = {pthr}, cluster p = {athr}\n\n", .trim=FALSE))

      if (bisided || twosided) {
        arg_list[["lower_thresh"]] <- lower_thresh
        arg_list[["upper_thresh"]] <- upper_thresh
      } else {
        arg_list[["one_thresh"]] <- one_thresh
        arg_list[["one_tail"]] <- voxelwise_stat$one_tail %||% "RIGHT_TAIL"
      }

      cobj <- do.call(afni_3dclusterize$new, arg_list)
      cobj$run(quiet = TRUE) # run 3dClusterize if needed

      if (isTRUE(add_whereami)) {
        cobj$add_whereami(atlases = whereami_atlases)
      }

      return(cobj)
    }
  )
)

#' R6 class for a list of 3dClustSim runs
#' @importFrom tibble tibble
#' @keywords internal
afni_3dclustsim_list <- R6::R6Class("afni_3dclustsim_list",
  private = list(
    clustsim_objs = NULL
  ),
  public = list(
    #' @description create a new afni_3dclustsim_list object
    #' @param obj_list A list of afni_3dclustsim objects.
    #' @param ... One or more afni_3dclustsim objects. Used if \code{obj_list} is NULL.
    initialize = function(obj_list=NULL, ...) {
      if (is.null(obj_list)) {
        # assume the ... contains a set of afni_3dclustsim objects
        obj_list <- list(...)
      }

      class_match <- sapply(obj_list, function(x) {
        checkmate::test_class(x, "afni_3dclustsim")
      })

      if (!all(class_match == TRUE)) {
        stop("At least one input is not a afni_3dclustsim object.")
      }

      private$clustsim_objs <- obj_list
    },

    #' @description submit all jobs in this list
    #' @param force If TRUE, pass force = TRUE to all the submit method of all subsidiary 3dClustSim objects
    submit = function(force = FALSE) {
      lapply(private$clustsim_objs, function(x) { x$submit(force = force) })
      return(invisible(self))
    },

    #' @description return the clustsim objects
    get_objs = function() {
      private$clustsim_objs
    }
  )
)
