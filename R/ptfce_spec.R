#' R6 class for pTFCE specification for one or more z-statistic images
#'
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_directory_exists assert_subset assert_string
#' @export
ptfce_spec <- R6::R6Class("ptfce_spec",
  private = list(
    z_files = NULL,
    mask_files = NULL,
    residuals_file = NULL,
    dof = NULL,
    fsl_smoothest_file = NULL,
    write_thresh_imgs = TRUE,
    smoothness_method = NULL,
    pvt_two_sided = TRUE,
    pvt_fwe_p = NULL,
    expect_list = NULL,
    scheduler = "slurm",
    time_per_zstat = "00:00:45",
    memgb_per_command = 1,
    pvt_cluster_objs = NULL,
    all_cope_df = NULL, # if gfeats are passed in, keep the contrast info for all z's as a data.frame
    set_fwep = function(vec) {
      checkmate::assert_numeric(vec, lower = 1e-10, upper = .9999, any.missing = FALSE, unique = TRUE, null.ok = FALSE)
      private$pvt_fwe_p <- vec
    },
    set_twosided = function(val) {
      checkmate::assert_logical(val, len = 1L, null.ok = FALSE)
      private$pvt_two_sided <- val
    },
    ptfce_outfile_from_zstat = function(z_file, what = "ptfce", fwe_p = NULL) {
      z_dir <- normalizePath(dirname(z_file))
      ext <- sub(".*(\\.nii(\\.gz)?)$", "\\1", z_file, perl = TRUE)
      base <- sub(ext, "", basename(z_file), fixed = TRUE)
      if (what == "ptfce") {
        if (is.null(fwe_p)) {
          return(file.path(z_dir, paste0(base, "_ptfce", ext)))
        } else if (!is.null(fwe_p)) {
          fwe_p <- format_fwe_probability(fwe_p)
          return(file.path(z_dir, paste0(base, "_ptfce_fwep_", fwe_p, ext)))
        }
      } else if (what == "csv") {
        return(file.path(z_dir, paste0(base, "_ptfce_zthresh.csv")))
      }
    },
    worker_args = function(i) {
      smoothness <- switch(private$smoothness_method,
        smoothest = c("--fsl_smoothest", private$fsl_smoothest_file[i]),
        residuals = c("--residuals", private$residuals_file[i], "--dof", private$dof[i]))
      c("--zstat", private$z_files[i], "--mask", private$mask_files[i],
        "--fwep", as.character(self$fwe_p), smoothness,
        if (private$write_thresh_imgs) "--write_thresh_imgs",
        if (self$two_sided) "--twosided" else "--onesided")
    },
    pending_maps = function(force) {
      checkmate::assert_flag(force)
      if (force) return(seq_along(private$z_files))
      complete <- vapply(private$expect_list, function(files) {
        if (!all(file.exists(unlist(files)))) return(FALSE)
        table <- tryCatch(utils::read.csv(files$z_csv), error = function(e) NULL)
        required <- c("p_value", "two_sided", "ptfce_version", "write_threshold_images")
        if (!all(required %in% names(table))) return(FALSE)
        # The submission host can prepare jobs without installing pTFCE locally.
        current_version <- if (requireNamespace("pTFCE", quietly = TRUE)) as.character(utils::packageVersion("pTFCE"))
        identical(sort(table$p_value), sort(as.numeric(self$fwe_p))) &&
          all(table$two_sided %in% self$two_sided) &&
          (is.null(current_version) || all(table$ptfce_version %in% current_version)) &&
          (!private$write_thresh_imgs || all(table$write_threshold_images %in% TRUE))
      }, logical(1))
      which(!complete)
    },
    set_expected_files = function() {
      f_list <- list()
      for (ii in seq_along(private$z_files)) {
        this_z <- private$z_files[ii]
        f_list[[ii]] <- list(
          # z_input = this_z,
          mask_input = private$mask_files[ii],
          z_ptfce = private$ptfce_outfile_from_zstat(this_z),
          z_csv = private$ptfce_outfile_from_zstat(this_z, what = "csv")
        )
        if (private$write_thresh_imgs) for (pp in private$pvt_fwe_p) {
          f_list[[ii]][[paste0("z_ptfce_fwep_", format_fwe_probability(pp))]] <- private$ptfce_outfile_from_zstat(this_z, fwe_p = pp)
        }
      }
      private$expect_list <- f_list
    }
  ),
  active = list(
    #' @field fwe_p a vector of p-values used for familywise error (FWE) z-statistic threshold calculations in pTFCE.
    fwe_p = function(vec) {
      if (missing(vec)) {
        return(private$pvt_fwe_p)
      } else {
        private$set_fwep(vec)
        private$set_expected_files()
      }
    },

    #' @field two_sided If \code{TRUE}, use native two-sided pTFCE on a shared threshold grid. If \code{FALSE},
    #'   only run pTFCE on the positive tail (z > 0).
    two_sided = function(val) {
      if (missing(val)) {
        return(private$pvt_two_sided)
      } else {
        private$set_twosided(val)
      }
    }
  ),
  public = list(
    #' @description Create a new ptfce_spec object
    #' @param gfeat_dir One or more .gfeat folders containing a higher-level FSL analysis. These will be used for zstat images,
    #'   mask files, and fsl residual smoothness estimates.
    #' @param zstat_numbers if a \code{gfeat_dir} is used, a vector of zstat numbers can also be provided to
    #'   subset the zstat images that are used in pTFCE correction. Ignored if gfeat_dir is not provided.
    #' @param fsl_smoothest_file The smoothness file created by the smoothest command in FSL. Created by FEAT automatically
    #'   or can be run manually with smoothest on the res4d file. Used by ptfce to compute resels.
    #' @param dof The degrees of freedom for the test of interest. Found in the dof file created by FSL and used by
    #'   ptfce for FWE correction.
    #' @param residuals_file The residuals file from the group analysis used by ptfce for FWE correction. In FSL, this
    #'   is the res4d.nii.gz file created by FEAT.
    #' @param z_files A vector of z-statistic filenames that should be corrected using pTFCE
    #' @param mask_files A vector of mask filenames that correspond to \code{z_files}. If this is of length 1,
    #'   then the mask file will be recycled for all zstat images.
    #' @param fwe_p A vector of p-values for which z-statistic thresholds will be calculated.
    #' @param two_sided If TRUE, p-values for \code{fwe_p} are treated as two-tailed (i.e., p-values are divided by 2 in
    #'   the TFCE z-threshold calculation.)
    #' @param write_thresh_imgs If TRUE, then pTFCE thresholds for each fwe_p will be applied to the TFCE image and saved
    #'   to the same folder as the z-statistic. These thresholded files let you look at the map at a given FWE threshold.
    #' @param scheduler Which scheduler to use for submitting jobs. Options are 'local', 'slurm', and 'torque'.
    #' @param time_per_zstat The amount of time to budget for each zstat to run through pTFCE in dd-hh:mm:ss format.
    #'   Default is 00:00:45 (45 seconds).
    #' @param memgb_per_command How many GB of memory should be requested for each pTFCE command/run. If not provided, defaults to 1.
    initialize = function(gfeat_dir = NULL, zstat_numbers = NULL, fsl_smoothest_file = NULL, dof = NULL, residuals_file = NULL,
                          z_files = NULL, mask_files = NULL, fwe_p = .05, two_sided = TRUE, write_thresh_imgs = TRUE,
                          scheduler = NULL, time_per_zstat = NULL, memgb_per_command = NULL) {
      if (!is.null(gfeat_dir)) {
        if (!is.null(z_files)) {
          stop("Cannot provide both z_files and gfeat_dir as inputs")
        }

        checkmate::assert_directory_exists(gfeat_dir)
        if (!is.null(zstat_numbers)) checkmate::assert_integerish(zstat_numbers)
        groups <- lapply(gfeat_dir, read_gfeat_dir, what = "stat_files", statistics = "zstat")
        maps <- lapply(groups, function(group) {
          lapply(group$cope_dirs, function(cope) {
            if (is.null(cope)) return(NULL)
            rows <- cope$cope_df[, c("cope_number", "contrast_name", "z"), drop = FALSE]
            if (!is.null(zstat_numbers)) rows <- rows[rows$cope_number %in% zstat_numbers, , drop = FALSE]
            rows <- rows[!is.na(rows$z), , drop = FALSE]
            if (!nrow(rows)) return(NULL)
            rows$mask <- cope$mask_file
            rows$smoothness <- cope$smoothness_file
            rows
          }) |> dplyr::bind_rows()
        }) |> dplyr::bind_rows()
        if (!nrow(maps)) stop("No matching z-statistic images in gfeat_dir.")
        private$all_cope_df <- maps[, c("cope_number", "contrast_name", "z"), drop = FALSE]
        private$z_files <- maps$z
        mask_files <- maps$mask
        fsl_smoothest_file <- maps$smoothness
      } else if (!is.null(z_files)) {
        checkmate::assert_file_exists(z_files)
        private$z_files <- z_files
      } else {
        stop("Either z_files or gfeat_dir must be provided as an input to ptfce_spec")
      }

      if (is.null(mask_files)) {
        stop("You must supply one or more mask files to be used in pTFCE calculations.")
      } else if (length(mask_files) == 1L && length(private$z_files) > 1L) {
        # replicate mask file for all inputs
        private$mask_files <- rep(mask_files, length(private$z_files))
      } else {
        stopifnot(length(mask_files) == length(private$z_files))
        private$mask_files <- mask_files
      }

      # verify that all mask files exist
      checkmate::assert_file_exists(private$mask_files)

      if (!is.null(residuals_file)) {
        checkmate::assert_file_exists(residuals_file)
        if (!length(residuals_file) %in% c(1L, length(private$z_files))) stop("Supply one residual file or one per map.")
        private$residuals_file <- rep_len(residuals_file, length(private$z_files))

        if (is.null(dof)) {
          stop("If a residuals file is provided, you also need to provide dof (either as an integer or a file containing that integer).")
        } else if (checkmate::test_integerish(dof)) {
          private$dof <- as.integer(dof)
        } else if (checkmate::test_file_exists(dof)) {
          private$dof <- as.integer(readLines(dof))
        } else {
          stop("Can't sort out this dof input")
        }

        if (!length(private$dof) %in% c(1L, length(private$z_files))) stop("Supply one dof or one per map.")
        private$dof <- rep_len(private$dof, length(private$z_files))
        private$smoothness_method <- "residuals"
      } else if (!is.null(fsl_smoothest_file)) {
        checkmate::assert_file_exists(fsl_smoothest_file)
        if (length(fsl_smoothest_file) == 1L) {
          # replicate smoothest for all inputs
          private$fsl_smoothest_file <- rep(fsl_smoothest_file, length(private$z_files))
        } else {
          stopifnot(length(fsl_smoothest_file) == length(private$z_files))
          private$fsl_smoothest_file <- fsl_smoothest_file
        }
        private$smoothness_method <- "smoothest"
      } else {
        message("Using smoothness from z-statistic image itself. Not recommended!")
        private$smoothness_method <- "zstat"
      }

      # populate fwep and twosided
      private$set_fwep(fwe_p)
      private$set_twosided(two_sided)

      if (!is.null(write_thresh_imgs)) {
        checkmate::assert_logical(write_thresh_imgs, len = 1L)
        private$write_thresh_imgs <- write_thresh_imgs
      }

      private$set_expected_files()

      if (!is.null(scheduler)) {
        checkmate::assert_string(scheduler)
        checkmate::assert_subset(scheduler, c("torque", "slurm", "local"))
        private$scheduler <- scheduler
      }

      if (!is.null(time_per_zstat)) {
        private$time_per_zstat <- validate_dhms(time_per_zstat)
      }

      if (!is.null(memgb_per_command)) {
        checkmate::assert_number(memgb_per_command)
        private$memgb_per_command <- memgb_per_command
      }
    },

    #' @description method to return calls to external ptfce_zstat.R script for each zstat
    #' @param include_complete if TRUE, return calls for zstats that already appear to have
    #'   pTFCE-corrected images in place. Default: FALSE.
    get_ptfce_calls = function(include_complete = FALSE) {
      script <- resolve_ptfce_worker_script()
      vapply(private$pending_maps(include_complete), function(i) {
        fwe_shell_command("Rscript", c("--vanilla", script, private$worker_args(i)))
      }, character(1L))
    },

    #' @description return the vector of expected output files
    get_expected_files = function() {
      private$expect_list
    },

    #' @description Run each pending map sequentially in the current R session.
    #' @param force if TRUE, re-run pTFCE on an existing output
    run = function(force = FALSE) {
      for (i in private$pending_maps(force)) ptfce_cli(private$worker_args(i))
      invisible(self)
    },

    #' method to submit all ptfce_zstat.R calls to a cluster based on the scheduler specified
    #' @param force if TRUE, re-run pTFCE for zstat images that already appear to have pTFCE-corrected outputs in place
    submit = function(force = FALSE) {
      checkmate::assert_flag(force)
      if (private$scheduler == "local") {
        self$run(force)
        return(invisible(character(0)))
      }
      # lookup ptfce calls for all inputs
      ptfce_calls <- self$get_ptfce_calls(include_complete = force)
      if (length(ptfce_calls) == 0L) {
        message("No pTFCE jobs need to be run for this input. If you want to recalculate completed outputs, use force=TRUE")
        return(invisible(character(0)))
      }

      # submit job to cluster
      child_job_ids <- fmri.pipeline::cluster_submit_shell_jobs(
        ptfce_calls,
        time_per_command = private$time_per_zstat,
        memgb_per_command = private$memgb_per_command, fork_jobs = TRUE,
        scheduler = private$scheduler,
        job_script_prefix = "job_ptfce"
      )
      return(invisible(child_job_ids))
    },

    #' @description returns \code{TRUE} when all expected outputs exist and their threshold tables match the requested settings
    is_complete = function() {
      length(private$pending_maps(FALSE)) == 0L
    },

    #' @description for each input file, obtain 3dClusterize objects that reflect the pTFCE-corrected clusters
    #' @param fwep the whole-brain familywise error rate to use (often .05)
    #' @param clust_nvox The minimum number of voxels to allow in a given cluster
    #' @param NN The cluster definition in AFNI terms. 1 = faces touch, 2 = edges touch, 3 = corners touch
    #' @param add_whereami if TRUE, lookup labels for each cluster
    #' @param whereami_atlases The atlases to request in the whereami lookup. If NULL, it uses the defaults
    #' @details Note that even though pTFCE enhances clusters, you still see some very small clusters in some cases.
    #'   Hence, there is no requirement to have clust_nvox > 1, but it may be a good idea for your sanity.
    get_clusters = function(fwep = .05, clust_nvox = 10, NN = 1L, add_whereami = TRUE, whereami_atlases = NULL) {
      if (!self$is_complete()) {
        stop("Cannot use get_cluster_mask until pTFCE has completed for all inputs.")
      }

      checkmate::assert_number(fwep, lower = 1e-10, upper = .999)
      checkmate::assert_integerish(clust_nvox, lower = 1, upper = 1e6)
      checkmate::assert_integerish(NN, len = 1L, lower = 1L, upper = 3L)
      checkmate::assert_logical(add_whereami, len = 1L)
      checkmate::assert_character(whereami_atlases, null.ok = TRUE)

      private$pvt_cluster_objs <- lapply(private$expect_list, function(input) {
        fwep_csv <- read.csv(input$z_csv)
        threshold <- ptfce_threshold_from_table(fwep_csv, fwep)
        bisided <- threshold$two_sided
        onesided <- !bisided
        z_thresh <- threshold$z
        lower_thresh <- upper_thresh <- one_thresh <- NULL

        if (isTRUE(bisided)) {
          lower_thresh <- -1*z_thresh
          upper_thresh <- z_thresh
        } else {
          one_thresh <- z_thresh
        }

        cobj <- afni_3dclusterize$new(
          threshold_file = input$z_ptfce, mask = input$mask_input, NN = NN, clust_nvox = clust_nvox,
          bisided = bisided, onesided = onesided, lower_thresh = lower_thresh,
          upper_thresh = upper_thresh, one_thresh = one_thresh
        )

        # x <- afni_whereami$new(afni_3dclusterize_obj = cobj)
        cobj$run(quiet = TRUE) # run 3dClusterize if needed
        if (isTRUE(add_whereami)) {
          cobj$add_whereami(atlases = whereami_atlases)
        }

        return(cobj)
      })

      # name elements by the original input
      names(private$pvt_cluster_objs) <- private$z_files

      if (!is.null(private$all_cope_df)) {
        tbl <- tibble(private$all_cope_df, cluster_obj = private$pvt_cluster_objs)
      } else {
        tbl <- tibble(cluster_obj = private$pvt_cluster_objs)
      }

      return(tbl)
    }
  )
)
