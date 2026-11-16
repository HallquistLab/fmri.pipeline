#' wrapper class for 3dClusterize
#' @importFrom tidyselect everything
#' @importFrom tidyr nest
#' @importFrom RNifti readNifti writeNifti
#' @export
afni_3dclusterize <- R6::R6Class("afni_3dclusterize",
  private = list(

    # private fields
    pvt_input_file = NULL, # internal record of the input file to be passed as -inset to 3dClusterize
    pvt_inset_file = NULL,
    pvt_mask = NULL,
    pvt_sided = "bi", # default to bi-sided (i.e., adjacent positive and negative clusters cannot merge)
    pvt_lower_thresh = NULL, # lower threshold value for clusterizing (two and bi)
    pvt_upper_thresh = NULL, # upper threshold value for clusterizing (two and bi)
    pvt_one_thresh = NULL, # single threshold value for clusterizing (one)
    pvt_one_tail = NULL, # which tail to test in case of one-sided (LEFT or RIGHT)
    pvt_mask_from_hdr = NULL,
    pvt_out_mask = NULL,
    pvt_threshold_file = NULL, # two-part alternative to inset file -- stitched together by this class before 3dClusterize
    pvt_data_file = NULL, # two-part alternative to inset file -- stitched together by this class before 3dClusterize
    pvt_ithr = NULL,
    pvt_idat = NULL,
    pvt_NN = 1L, # faces + edges + corners must touch
    pvt_clust_nvox = 10, # very loose lower bound by default (tiny clusters)
    pvt_pref_map = NULL,
    pvt_pref_dat = NULL,
    pvt_1Dformat = TRUE, # currently not used since parser depends on 1D format
    pvt_quiet = NULL,
    pvt_orient = "LPI", #AFNI default is RAI, which seems like a silly default to me
    pvt_binary = FALSE,
    pvt_combine_call = NULL, # If needed, the 3dTcat call for combining pvt_threshold_file and pvt_data_file
    pvt_clusterize_call = NULL,
    pvt_clusterize_output_file = NULL,
    pvt_dirty_call = FALSE, # whether object fields have changed since last build_call
    pvt_whereami = NULL,
    pvt_has_clusters = NULL, # cached TRUE/FALSE for whether 3dClusterize found clusters (NULL if no output)
    pvt_clust_df = NULL, # cached data.frame of parsed 3dClusterize output
    pvt_subclust_list = NULL, # cached list of afni_3dclusterize objects for each parent cluster larger than threshold
    pvt_subclust_df = NULL, # cached data.frame of subcluster information
    pvt_atlas_files = list(), # location of atlas overlap files from calls to subset_atlas_against_clusters
    pvt_subclust_files = list(), # location of subclustering outputs

    # private methods
    # ---------------
    # methods to set private fields (used by active bindings)
    set_lower_thresh = function(val) {
      checkmate::assert_number(val)
      private$pvt_lower_thresh <- val
    },
    set_upper_thresh = function(val) {
      checkmate::assert_number(val)
      private$pvt_upper_thresh <- val
    },
    set_one_thresh = function(val) {
      checkmate::assert_number(val)
      private$pvt_one_thresh <- val
    },
    set_mask = function(mask_file) {
      checkmate::assert_file_exists(mask_file)
      private$pvt_mask <- mask_file
    },
    set_sided = function(val) {
      checkmate::assert_string(val)
      checkmate::assert_subset(val, c("1", "one", "2", "two", "bi"))
      if (val == "1") {
        val <- "one" # convert to standard nomenclature
      } else if (val == "2") {
        val <- "two"
      }
      private$pvt_sided <- val
      private$pvt_dirty_call <- TRUE # will rebuild call
    },
    set_NN = function(val) {
      if (checkmate::test_string(val)) {
        val <- as.integer(val)
      }

      checkmate::assert_integerish(val, lower = 1, upper = 3, len = 1L, any.missing = FALSE)
      private$pvt_NN <- as.integer(val)
      private$pvt_dirty_call <- TRUE # will rebuild call
    },
    set_clust_nvox = function(val) {
      if (checkmate::test_string(val)) {
        val <- as.integer(val)
      }

      checkmate::assert_integerish(val, lower = 1, upper = 1e5, len = 1L, any.missing = FALSE)
      private$pvt_clust_nvox <- as.integer(val)
      if (!is.null(private$pvt_clust_vol)) {
        message("Resetting clust_vol to NULL because clust_nvox was set.")
        private$pvt_clust_vol <- NULL
      }
      private$pvt_dirty_call <- TRUE # will rebuild call
    },
    set_clust_vol = function(val) {
      if (checkmate::test_string(val)) {
        val <- as.numeric(val)
      }

      checkmate::assert_numeric(val, lower = 0.5, upper = 1e6, len = 1L, any.missing = FALSE)
      private$pvt_clust_vol <- as.numeric(val)
      if (!is.null(private$pvt_clust_nvox)) {
        message("Resetting clust_nvox to NULL because clust_vol was set.")
        private$pvt_clust_nvox <- NULL
      }
      private$pvt_dirty_call <- TRUE # will rebuild call
    },
    set_clusterize_output_file = function(val) {
      if (missing(val) || is.null(val)) { return(invisible(NULL)) } # nothing to do
      checkmate::assert_string(val)
      private$pvt_clusterize_output_file <- val
    },
    set_pref_map = function(val) {
      if (is.null(val)) {
        private$pvt_pref_map <- val
      } else {
        checkmate::assert_string(val)
        ext <- file_ext(val)
        if (is.na(ext)) {
          val <- paste0(val, ".nii.gz") # force nifti output
        } else if (grepl("(brik|head).*", ext)) {
          message("Forcing -pref_map extension to .nii.gz for consistency")
          val <- paste0(file_sans_ext(val), ".nii.gz") # force nifti output
        }
        private$pvt_pref_map <- val
      }
    },
    set_pref_dat = function(val) {
      if (is.null(val)) {
        private$pvt_pref_dat <- val
      } else {
        checkmate::assert_string(val)
        ext <- file_ext(val)
        if (is.na(ext)) {
          val <- paste0(val, ".nii.gz") # force nifti output
        } else if (grepl("(brik|head).*", ext)) {
          message("Forcing -pref_dat extension to .nii.gz for consistency")
          val <- paste0(file_sans_ext(val), ".nii.gz") # force nifti output
        }
        private$pvt_pref_dat <- val
      }
    },

    # merge clust_df against whereami and subcluster data
    merge_aux_data = function(include_whereami = TRUE, include_subclusters = TRUE, include_overlap = TRUE) {
      if (isFALSE(self$is_complete)) { return(invisible(NULL)) } # fail gracefully

      clust_df <- private$pvt_clust_df
      if (isTRUE(include_whereami) && private$has_whereami()) {
        clust_df <- clust_df %>% dplyr::left_join(self$whereami$get_whereami_df(), by = "roi_num")
      }

      # for now, return overlap data as a nested list-column in the clust_df, which can be unnested later, if valuable
      if (isTRUE(include_overlap) && private$has_whereami()) {
        overlap_df <- self$whereami$get_overlap_df() %>% group_by(roi_num) %>% nest(overlap = -roi_num)
        clust_df <- clust_df %>% dplyr::left_join(overlap_df, by="roi_num")
      }

      if (isTRUE(include_subclusters) && is.data.frame(private$pvt_subclust_df) && "roi_num" %in% names(private$pvt_subclust_df)) {
        clust_df <- clust_df %>%
          dplyr::bind_rows(private$pvt_subclust_df) %>%
          dplyr::select("roi_num", "subroi_num", everything()) %>%
          dplyr::arrange(roi_num, subroi_num)
      }

      return(clust_df)
    },

    # determine whether any clusters were found
    check_for_clusters = function(quiet = FALSE) {
      if (checkmate::test_file_exists(private$pvt_clusterize_output_file)) {
        clust_txt <- readLines(private$pvt_clusterize_output_file, n = 3)
        if (grepl("#** NO CLUSTERS FOUND ***", clust_txt[1L], fixed = TRUE)) {
          if (!quiet) warning("No clusters found by 3dClusterize")
          private$pvt_has_clusters <- FALSE
          private$pvt_clust_df <- data.frame()
        } else {
          private$pvt_has_clusters <- TRUE
        }
      } else {
        private$pvt_has_clusters <- NULL # output file not in place yet
      }
    },

    # whether a whereami object has already been added to this object
    has_whereami = function() {
      !is.null(private$pvt_whereami) && inherits(private$pvt_whereami, "afni_whereami")
    },

    get_combine_call = function() {

    },
    build_call = function() {
      if (isFALSE(private$pvt_dirty_call)) return(invisible(NULL)) # no need to rebuild call

      str <- glue("3dClusterize -overwrite -inset \"{private$pvt_input_file}\" -ithr {private$pvt_ithr} -NN {private$pvt_NN} -1Dformat")
      if (!is.null(private$pvt_mask)) str <- glue("{str} -mask \"{private$pvt_mask}\"")
      if (isTRUE(private$pvt_mask_from_hdr)) str <- glue("{str} -mask_from_hdr")
      if (!is.null(private$pvt_out_mask)) str <- glue("{str} -out_mask \"{private$pvt_out_mask}\"")

      # I do not support -within_range yet...

      if (!is.null(private$pvt_idat)) str <- glue("{str} -idat {private$pvt_idat}")
      if (private$pvt_sided == "one") {
        str <- glue("{str} -1sided {private$pvt_one_tail} {private$pvt_one_thresh}")
      } else if (private$pvt_sided == "two") {
        str <- glue("{str} -2sided {private$pvt_lower_thresh} {private$pvt_upper_thresh}")
      } else if (private$pvt_sided == "bi") {
        str <- glue("{str} -bisided {private$pvt_lower_thresh} {private$pvt_upper_thresh}")
      }

      if (!is.null(private$pvt_clust_nvox)) {
        str <- glue("{str} -clust_nvox {private$pvt_clust_nvox}")
      } else if (!is.null(private$pvt_clust_vol)) {
        str <- glue("{str} -clust_vol {private$pvt_clust_vol}")
      } else {
        stop("Both clust_nvox and clust_vol are missing... This shouldn't happen!")
      }

      if (!is.null(private$pvt_pref_map)) str <- glue("{str} -pref_map \"{private$pvt_pref_map}\"")
      if (!is.null(private$pvt_pref_dat)) str <- glue("{str} -pref_dat \"{private$pvt_pref_dat}\"")
      if (isTRUE(private$pvt_quiet)) str <- glue("{str} -quiet")
      if (!is.null(private$pvt_orient)) str <- glue("{str} -orient {private$pvt_orient}")
      if (isTRUE(private$pvt_binary)) str <- glue("{str} -binary")

      str <- glue(str, " > \"{private$pvt_clusterize_output_file}\"")

      private$pvt_clusterize_call <- str
    }
  ),

  # active bindings (these are the only things that can change after initialization)
  active = list(
    #' @field sided Whether to clusterize 1-sided ('one'), two-sided ('two'), or bi-sided ('bi')
    sided = function(val) {
      if (missing(val)) {
        return(private$pvt_sided)
      } else {
        private$set_sided(val)
      }
    },

    #' @field NN the cluster definition basis: 1 = faces touch; 2 = edges touch; 3 = corners touch
    NN = function(val) {
      if (missing(val)) {
        return(private$pvt_NN)
      } else {
        private$set_NN(val)
      }
    },

    #' @field clust_nvox The minimum number of voxels required for each cluster
    clust_nvox = function(val) {
      if (missing(val)) {
        return(private$pvt_clust_nvox)
      } else {
        private$set_clust_nvox(val)
      }
    },

    #' @field clust_vol The minimum volume (in microliters) required for each cluster
    clust_vol = function(val) {
      if (missing(val)) {
        return(private$pvt_clust_vol)
      } else {
        private$set_clust_vol(val)
      }
    },

    #' @field lower_thresh For two/bi-sided clusterizing, the lower threshold for the left tail of the distribution
    lower_thresh = function(val) {
      if (missing(val)) {
        return(private$pvt_lower_thresh)
      } else {
        private$set_lower_thresh(val)
      }
    },

    #' @field upper_thresh For two/bi-sided clusterizing, the upper threshold for the right tail of the distribution
    upper_thresh = function(val) {
      if (missing(val)) {
        return(private$pvt_upper_thresh)
      } else {
        private$set_upper_thresh(val)
      }
    },

    #' @field one_thresh For one-sided clusterizing, the threshold for the test statistic distribution
    one_thresh = function(val) {
      if (missing(val)) {
        return(private$pvt_one_thresh)
      } else {
        private$set_one_thresh(val)
      }
    },

    #' @field mask the mask within which 3dClusterize searches for clusters
    mask = function(val) {
      if (missing(val)) {
        return(private$pvt_mask)
      } else {
        private$set_mask(val)
      }
    },

    #' @field pref_map The name/location of the pref_map (aka cluster_map) file containing an integer-valued mask of identified clusters
    pref_map = function(val) {
      if (missing(val)) {
        return(private$pvt_pref_map)
      } else {
        private$set_pref_map(val)
      }
    },

    #' @field pref_dat The name/location of the pref_data (aka cluster_masked_data) file containing the input data masked by the clusters
    pref_dat = function(val) {
      if (missing(val)) {
        return(private$pvt_pref_dat)
      } else {
        private$set_pref_dat(val)
      }
    },

    #' @field clusterize_output_file The name of the text file containing the output of 3dClusterize (i.e., the table of clusters)
    clusterize_output_file = function(val) {
      if (missing(val)) {
        return(private$pvt_clusterize_output_file)
      } else {
        private$set_clusterize_output_file(val)
      }
    },

    #' @field whereami passthrough access to whereami object if that has been setup
    whereami = function(val) {
      if (missing(val)) {
        if (is.null(private$pvt_whereami)) {
          message("No whereami has been added to this object yet. Use $add_whereami() to do so.")
          return(invisible(NULL))
        } else {
          private$pvt_whereami
        }
      } else {
        checkmate::assert_class(val, "afni_whereami", null.ok = TRUE)
        private$pvt_whereami <- val
      }
    }

  ),

  public = list(
    #' @description initialization function for a new afni_3dclusterize object. Arguments largely mirror the 3dClusterize parameters.
    #' @param inset A 4D dataset containing the statistic to use for thresholding (ithr) and, optionally, the data value to output/retain
    #' @param mask If specified, the volume will be masked by \code{mask} prior to clusterizing
    #' @param threshold_file A 3D dataset containing the statistic to use for thresholding. Mutually exclusive with \code{inset}
    #'   If passed, \code{ithr} and \code{idat} are ignored because the \code{inset} file is generated internally.
    #' @param data_file A 3D dataset containing the data value to be retained in clusters post-thresholding.
    #'   Must be passed with \code{threshold_file} and will be stitched together with it internally. Mutually exclusive with \code{inset}.
    #' @param mask_from_hdr passes through as -mask_from_hdr
    #' @param out_mask passes through as -out_mask
    #' @param ithr sub-brik number for the voxelwise threshold. Passes through as -ithr
    #' @param idat sub-brik number for the voxelwise data to be output in cluster table. Passes through as -idat
    #' @param onesided if TRUE, clusterizing will be conducted on one tail of the statistic distribution (-ithr)
    #' @param twosided if TRUE, clusterizing will be conducted on both tails of the statistic distribution (-ithr)
    #' @param bisided if TRUE, clusterizing will be conducted on each tail of the distribution individually
    #' @param lower_thresh the lower tail cutoff for two/bi-sided testing
    #' @param upper_thresh the upper tail cutoff for two/bi-sided testing
    #' @param one_thresh The threshold value for one-sided testing
    #' @param one_tail For one-sided clusterizing, whether to threshold the LEFT or RIGHT tail of the distribution
    #' @param clust_nvox The minimum number of voxels allowed in a cluster. Passes through as -clust_nvol
    #' @param clust_vol The minimum volume in (microliters) allowed in a cluster (mutually exclusive with clust_nvox). 
    #'   Passes through as -clust_vol
    #' @param pref_map File name for the integer-valued mask containing each cluster, ordered by descending voxel size.
    #'   Passes through as -pref_map.
    #' @param pref_dat File name for the clusterized and thresholded data. Passes through as -pref_dat.
    #' @param NN 1, 2, 3. Default: 1. Passes through as -NN.
    #' @param quiet passes through as -quiet.
    #' @param orient passes through as -orient. 'RAI' or 'LPI'. Default is LPI.
    #' @param binary if TRUE, the pref_map (cluster mask) will be output as a 1/0 binary image instead of integer-valued. 
    #'   Passes through as -binary.
    #' @param clusterize_output_file The name/location of the 3dClusterize output file containing a table of identified clusters.
    #'   Defaults to adding the suffix '_clusters.1D' to the input image and placing the file in the same folder as the input.
    initialize = function(inset = NULL, mask = NULL, threshold_file = NULL, data_file = NULL, mask_from_hdr = NULL, out_mask = NULL, 
      ithr = NULL, idat = NULL, onesided = NULL, twosided = NULL, bisided = NULL, 
      lower_thresh = NULL, upper_thresh = NULL, one_thresh = NULL, one_tail = NULL,
      NN = NULL, clust_nvox = NULL, clust_vol = NULL, pref_map = "default", pref_dat = "default",
      quiet = NULL, orient = NULL, binary = NULL, clusterize_output_file = NULL) {

      if (!is.null(inset)) {
        checkmate::assert_file_exists(inset)
        private$pvt_inset_file <- inset #normalizePath(inset)
      } else if (!is.null(threshold_file)) {
        checkmate::assert_file_exists(threshold_file)
        private$pvt_threshold_file <- normalizePath(threshold_file)

        if (!is.null(data_file)) {
          checkmate::assert_file_exists(data_file)
          private$pvt_data_file <- normalizePath(data_file)
          ithr <- 0 # always put threshold image as first sub-brik in the internal file
          idat <- 1
        } else {
          ithr <- 0
          idat <- NULL
        }
      } else {
        stop("You must pass either inset (combined threshold + data file) or data_file and threshold_file.")
      }

      if (!is.null(mask)) {
        checkmate::assert_file_exists(mask)
        private$pvt_mask <- mask
      }

      if (!is.null(mask_from_hdr)) {
        checkmate::assert_logical(mask_from_hdr, len = 1L)
        private$pvt_mask_from_hdr <- mask_from_hdr
      }

      if (!is.null(out_mask)) {
        checkmate::assert_string(out_mask, len = 1L)
        private$pvt_out_mask <- out_mask
      }

      if (is.null(ithr)) { ithr <- 0 } # default to first sub-brik
      checkmate::assert_integerish(ithr, lower=0, upper=1e4)
      private$pvt_ithr <- ithr

      if (!is.null(idat)) {
        checkmate::assert_integerish(idat, lower = 0, upper = 1e4)
        private$pvt_idat <- idat
      }

      n_sided <- sum(sapply(list(onesided, twosided, bisided), function(x) isTRUE(x)))
      if (n_sided > 1L) {
        stop("Only one 'sided' option can be specified: onesided, twosided, or bisided.")
      } else if (isTRUE(onesided)) {
        checkmate::assert_logical(onesided, len = 1L)
        private$set_sided("one")
        if (is.null(one_thresh)) {
          stop("Must pass in a one_thresh value used for thresholding image prior to clusterizing")
        } else {
          checkmate::assert_number(one_thresh)
          private$pvt_one_thresh <- one_thresh
        }

        if (is.null(one_tail)) {
          stop("Must pass in a one_tail value -- LEFT or RIGHT.")
        } else {
          checkmate::assert_subset(one_tail, c("LEFT", "RIGHT", "LEFT_TAIL", "RIGHT_TAIL"), empty.ok = FALSE)
          private$pvt_one_tail <- one_tail
        }
      } else if (isTRUE(twosided)) {
        checkmate::assert_logical(twosided, len = 1L)
        private$set_sided("two")
      } else if (isTRUE(bisided)) {
        checkmate::assert_logical(bisided, len = 1L)
        private$set_sided("bi")
      }

      if (private$pvt_sided %in% c("two", "bi")) {
        private$set_lower_thresh(lower_thresh)
        private$set_upper_thresh(upper_thresh)
      }

      # set NN active binding
      if (!is.null(NN)) {
        private$set_NN(NN)
      }

      # set clust_nvox active binding
      if (is.null(clust_nvox) && is.null(clust_vol)) {
        stop("Both clust_nvox and clust_vol were not provided. Cannot move forward with clusterize!")
      } else if (!is.null(clust_nvox)) {
        private$set_clust_nvox(clust_nvox)
      } else if (!is.null(clust_vol)) {
        if (!is.null(clust_nvox)) {
          stop("Cannot pass both clust_nvox and clust_vol!")
        }
        private$set_clust_vol(clust_vol)
      }

      if (!is.null(quiet)) {
        checkmate::assert_logical(quiet, len = 1L)
        private$pvt_quiet <- quiet
      }

      if (!is.null(binary)) {
        checkmate::assert_logical(binary, len = 1L)
        private$pvt_binary <- binary
      }

      if (!is.null(orient)) {
        checkmate::assert_string(orient)
        orient <- toupper(orient)
        checkmate::assert_subset(orient, c("LPI", "RAI"))
        # should probably check LPI, RAI, etc.
        private$pvt_orient <- orient
      }

      # set output file for cluster table
      self$clusterize_output_file <- clusterize_output_file

      # sort out the -inset to be used in the 3dClusterize call (if a threshold and data file are passed)
      if (!is.null(private$pvt_threshold_file)) {
        default_wd <- dirname(private$pvt_threshold_file)
        if (is.null(private$pvt_data_file)) {
          private$pvt_input_file <- private$pvt_threshold_file # 3D input to clusterize
        } else {
          # build tmp dataset that combined threshold and data file
          private$pvt_input_file <- tempfile(pattern = "clustinput", fileext = ".nii.gz")
          private$pvt_combine_call <- glue("3dTcat -output {private$pvt_input_file} {private$pvt_threshold_file} {private$pvt_data_file}")
        }
      } else {
        default_wd <- dirname(private$pvt_inset_file)
        private$pvt_input_file <- private$pvt_inset_file
      }

      # set default working directory for outputs if not specified
      # private$set_default_wd()

      # default to naming clusters file according to the inset file
      if (is.null(private$pvt_clusterize_output_file)) {
        private$pvt_clusterize_output_file <- paste0(basename(file_sans_ext(private$pvt_input_file)), "_clusters.1D")
      }

      # if output files are provided without any path specification, use the directory of the input file
      private$pvt_clusterize_output_file <- R.utils::getAbsolutePath(private$pvt_clusterize_output_file,
        workDirectory = default_wd
      )

      
      # set the integer-valued cluster mask output name
      if (checkmate::test_string(pref_map)) {
        if (pref_map == "default") {
          # output a default pref map based on the input file name
          self$pref_map <- paste0(basename(file_sans_ext(private$pvt_input_file)), "_clusters.nii.gz")
        } else {
          self$pref_map <- pref_map
        }
      } else if (!is.null(pref_map)) {
        stop("Don't know how to interpret pref_map input")
      }

      # set the data output, masked by the clusterize result
      if (checkmate::test_string(pref_dat)) {
        if (is.null(private$pvt_idat)) {
          # no specific data sub-brik provided for the -pref_dat. Assume that it should be the same sub-brik as the threshold (-ithr).
          private$pvt_idat <- private$pvt_ithr
        }

        if (pref_dat == "default") {
          # output a default pref dat based on the input file name
          self$pref_dat <- paste0(basename(file_sans_ext(private$pvt_input_file)), "_cluster_data.nii.gz")
        } else {
          self$pref_dat <- pref_dat
        }
      } else if (!is.null(pref_dat)) {
        stop("Don't know how to interpret pref_dat input")
      }

      if (!is.null(private$pvt_pref_map)) {
        private$pvt_pref_map <- R.utils::getAbsolutePath(private$pvt_pref_map, workDirectory = default_wd)
      }

      if (!is.null(private$pvt_pref_dat)) {
        private$pvt_pref_dat <- R.utils::getAbsolutePath(private$pvt_pref_dat, workDirectory = default_wd)
      }

      # look at whether 3dClusterize has already run and, if so, populate whether clusters were found
      private$check_for_clusters()

    },

    #' @description run the 3dClusterize command relevant to this object
    #' @param force if TRUE, 3dClusterize will be re-run
    #' @param quiet if TRUE, don't output messages as object is run or checked
    run = function(force = FALSE, quiet = FALSE) {
      checkmate::assert_logical(force, len = 1L)
      checkmate::assert_logical(quiet, len = 1L)

      if (self$is_complete() && isFALSE(force)) {
        if (isFALSE(quiet)) {
          message("We will not re-run 3dClusterize because it is already complete. Use $run(force = TRUE) to re-run.")
        }
        return(invisible(NULL))
      }

      private$build_call()
      run_afni_command(private$pvt_clusterize_call, echo = !quiet, ignore.stderr = quiet)
      self$reset_cache() # if we have run 3dClusterize, nullify the cached whereami and clust_df objects so that these do not persist
      private$check_for_clusters() # populate whether clusters were found after run  
    },

    #' @description return the 3dClusterize table of clusters as a data.frame
    #' @details This function will return an empty data.frame if the 3dClusterize output file cannot be found.
    #' @param include_whereami If TRUE and if $add_whereami() is already complete, merge the whereami data
    #'   into the cluster data.frame that is returned by this function.
    #' @param include_subclusters If TRUE and if $generate_subclusters() is already complete, merge the subcluster
    #'   data into the cluster data.frame that is returned by this function.
    #' @param include_overlap If TRUE and if $add_whereami() is already complete, merge the mask overlap data into
    #'   the cluster data.frame as a nested list-column called overlap (since there are many rows for each ROI)
    get_clust_df = function(include_whereami = TRUE, include_subclusters = TRUE, include_overlap = TRUE) {
      if (!checkmate::test_file_exists(private$pvt_clusterize_output_file)) {
        warning(glue("The expected 3dClusterize output does not exist: {private$pvt_clusterize_output_file}"))
        return(invisible(data.frame()))
      } else if (is.data.frame(private$pvt_clust_df)) { # use cached object
        return(private$merge_aux_data(include_whereami, include_subclusters, include_overlap))
      }

      # populate whether clusters exist based on output file
      private$check_for_clusters()

      if (isFALSE(private$pvt_has_clusters)) {
        return(data.frame())
      }

      clust_txt <- readLines(private$pvt_clusterize_output_file)
      clust_df <- read.table(private$pvt_clusterize_output_file)
      comment_lines <- grep("^\\s*#", clust_txt, perl = TRUE)
      table_lines <- grep("^\\s*#", clust_txt, perl = TRUE, invert = TRUE)

      # header
      cols <- clust_txt[min(table_lines) - 2L] # header is 2 lines before table starts
      cols <- make.names(strsplit(sub("^#", "", cols), "\\s{2,}")[[1L]]) # columns are separated by 2+ spaces
      stopifnot(length(cols) == ncol(clust_df))
      names(clust_df) <- cols

      clust_df$roi_num <- seq_len(nrow(clust_df))
      clust_df <- clust_df %>% dplyr::select("roi_num", tidyselect::everything()) # place roi_num first

      parse_header_rows <- function(lines, clust_df) {
        attrib <- lapply(lines, function(x) {
          # don't allow errant lines that lack the comment #; also only parse lines with key - value pairing (separated by =)
          if (!grepl("^#", x) || !grepl("=", x, fixed=TRUE)) {
            return(NULL)
          } else {
            keyval_split <- strsplit(x, "=", fixed = TRUE)[[1L]]
            if (length(keyval_split) != 2) {
              warning("Cannot sort out how to parse 3dClusterize header line: ", x)
              return(NULL)
            } else {
              key <- trimws(sub("^\\s*#\\s*\\[\\s*(.+)", "\\1", keyval_split[[1]], perl = TRUE))
              value <- trimws(sub("\\s*(.*)\\]\\s*$", "\\1", keyval_split[[2]], perl = TRUE))
              return(value %>% setNames(key))
            }
          }
        })

        # add each attribute to the data.frame
        for (aa in attrib) {
          if (!is.null(aa)) {
            attr(clust_df, names(aa)) <- aa
          }
        }

        #attr(clust_df, key) <- value
        return(clust_df)
      }

      comment_txt <- clust_txt[comment_lines]

      # There are a couple of weird lines for threshold value and nvoxel threshold. These allow for multiple = signs in the values.
      # Here, we pre-parse these lines and make them into single key-value pairs per line by adding one line for each.
      eq_lens <- sapply(gregexpr("=", text = comment_txt, fixed = TRUE), length)
      if (any(eq_lens > 1L)) {
        mult_eq_txt <- comment_txt[eq_lens > 1]

        # remove these lines from comment_txt (will be replaced momentarily by expanded copies)
        comment_txt <- comment_txt[eq_lens == 1L]
        for (line in mult_eq_txt) {
          firsteq <- gregexpr("=", text = line, fixed = TRUE)[[1]][1]
          stem <- trimws(substr(line, 1, firsteq - 1))
          values <- substr(line, firsteq + 1, nchar(line))
          values_split <- trimws(strsplit(values, ";", fixed = TRUE)[[1]])
          has_subkeys <- all(grepl("=", values_split, fixed = TRUE)) # TRUE for lines like left-tail thr=-3.000;  right-tail thr=3.000
          if (isTRUE(has_subkeys)) {
            sub_split <- strsplit(values_split, "=", fixed=TRUE)
            for (ii in seq_along(sub_split)) {
              # combine overall stem with subkey-specific stem, separated by comma. Add " ]" for all items except the last (has it already)
              comment_txt <- c(
                comment_txt,
                paste0(stem, ", ", sub_split[[ii]][1], " = ", sub_split[[ii]][2], ifelse(ii == length(sub_split), "", " ]"))
              )
            }
          } else {
            # really just two lines have been piled into one, like #[ Nvoxel threshold    = 35;  Volume threshold = 425.845 ]
            # split these into two complete lines
            for (ii in seq_along(values_split)) {
              if (ii == 1L) { # for first split, key has been pulled into stem
                comment_txt <- c(comment_txt, paste0(stem, " = ", values_split[ii], " ]"))
              } else { # for other splits, key is still to the left of the equals sign
                comment_txt <- c(comment_txt, paste0("#[ ", values_split[ii]))
              }
            }
          }
        }
      }

      # add header row information as attributes to data.frame, cache parsed data.frame object
      private$pvt_clust_df <- parse_header_rows(comment_txt, clust_df)

      return(private$merge_aux_data(include_whereami, include_subclusters, include_overlap))
    },

    #' @description method to read and return the integer-valued clusterized mask (aka -pref_map) as an oro.nifti object
    #' @return an oro.nifti object containing the clusterized mask from 3dClusterize
    get_cluster_map_nifti = function() {
      if (checkmate::test_file_exists(private$pvt_pref_map)) {
        oro.nifti::readNIfTI(private$pvt_pref_map, reorient=FALSE)
      } else {
        NULL # return NULL if file is missing
      }
    },

    #' @description returns the 3dClusterize call for this specification
    get_call = function() {
      private$build_call()
      return(private$pvt_clusterize_call)
    },

    #' @description Provides a vector of expected output files that correspond to this 3dClusterize setup
    #' @param exclude_missing if TRUE (default), any output file that cannot be found will be returned as NA.
    #' @return a named vector of output files related to this 3dClusterize setup
    get_outputs = function(exclude_missing = TRUE) {

      cluster_map <- private$pvt_pref_map
      cluster_masked_data <- private$pvt_pref_dat
      cluster_table <- private$pvt_clusterize_output_file

      if (isTRUE(exclude_missing)) {
        if (!checkmate::test_file_exists(cluster_map)) cluster_map <- NA_character_
        if (!checkmate::test_file_exists(cluster_masked_data)) cluster_masked_data <- NA_character_
        if (!checkmate::test_file_exists(cluster_table)) cluster_table <- NA_character_
      }

      atlas_files <- private$pvt_atlas_files
      subcluster_files <- private$pvt_subclust_files

      named_list(cluster_table, cluster_map, cluster_masked_data, atlas_files, subcluster_files)

    },

    #' @description returns the orientation code for this 3dClusterize call (LPI or RAI)
    get_orient = function() {
      private$pvt_orient
    },

    #' @description Add's an afni_whereami object to this class in the $whereami slot. The corresponding
    #'   whereami command is also run when this is added so that coordinates and labels can be obtained
    #'   immediately. To access the whereami object and its methods, use $whereami()
    #' @param atlases An optional character vector of atlases to be requested in whereami.
    add_whereami = function(atlases=NULL) {
      if (private$has_whereami()) {
        message("whereami object already added to this clusterize object. Use $whereami() to access it.")
      } else if (!self$is_complete()) {
        message("Cannot add whereami to 3dclusterize object until clusterization is run!")
      } else {
        self$whereami <- afni_whereami$new(
          afni_3dclusterize_obj = self, atlases = atlases
        )

        self$whereami$run(force = TRUE)
      }
    },

    #' @description returns TRUE if all expected output files exist for this 3dClusterize call
    is_complete = function() {
      ofiles <- self$get_outputs(exclude_missing = TRUE)
      file_status <- c(
        table = checkmate::test_file_exists(ofiles$cluster_table),
        pref_map = checkmate::test_file_exists(ofiles$cluster_map),
        pref_dat = checkmate::test_file_exists(ofiles$cluster_masked_data)
      )
      required <- c(table=TRUE, pref_map=FALSE, pref_dat=FALSE)
      if (!is.null(private$pvt_pref_map)) required["pref_map"] <- TRUE
      if (!is.null(private$pvt_pref_dat)) required["pref_dat"] <- TRUE

      # if there are no clusters, then the object is complete with a .1D alone
      if (!isTRUE(private$pvt_has_clusters)) {
        required[c("pref_map", "pref_dat")] <- FALSE
      }

      all(file_status[required])
    },

    #' @description break up large clusters into subclusters
    #' @param break_nvox Break up any clusters larger than this value into subclusters. Default: 400
    #' @param min_subclust_nvox The smallest number of voxels allowed in a subcluster. Default: 25.
    #' @param max_subclust_nvox The largest numver of voxels allowed in a subcluster. If NULL, no upper limit is set.
    #' @param min_n_subclust The smallest number of subclusters that will be allowed. Must be 2 or greater. Default: 2
    #' @param max_n_subclust The maximum number of subclusters that will be allowed. If NULL, no upper limit is set.
    #' @param step_size The step size used to change the threshold values in the test statistic map being clusterized. Default: 0.1.
    #' @param max_iter The maximum number of steps to be taken for subcluster search. Default: 50.
    #' @param add_whereami If TRUE, whereami will be run for each subcluster. Default: TRUE
    #' @param whereami_atlases Passes through to afni_whereami for specifying which atlases to use in lookup
    #' @param print_progress If TRUE, the user will see the thresholds being used to subcluster each region.
    generate_subclusters = function(break_nvox = 400, min_subclust_nvox = 25, max_subclust_nvox = NULL,
      min_n_subclust = 2, max_n_subclust = NULL, step_size = .1, max_iter = 50,
      add_whereami = TRUE, whereami_atlases = NULL, print_progress = FALSE) {

      if (is.null(private$pvt_has_clusters)) {
        warning("Did not find expected 3dClusterize output file to use for clusters. Have you used $run() yet?")
        return(invisible(self))
      } else if (isFALSE(private$pvt_has_clusters)) {
        warning("No clusters were found for this object. Cannot create subclusters.")
        return(invisible(self))
      }

      if (is.null(max_subclust_nvox)) {
        max_subclust_nvox <- Inf # no limit on subcluster size
      }

      checkmate::assert_integerish(break_nvox, lower = 10, upper = 1e5, len = 1L)
      checkmate::assert_integerish(min_subclust_nvox, lower = 1, upper = 1e4, len = 1L)
      if (!is.infinite(max_subclust_nvox)) checkmate::assert_integerish(max_subclust_nvox, lower = 5, upper = 1e6, len = 1L)

      if (is.null(max_n_subclust)) {
        max_n_subclust <- Inf # makes the upper limit unbounded
      }

      checkmate::assert_integerish(min_n_subclust, lower = 2L, upper = 1000L)
      if (!is.infinite(max_n_subclust)) checkmate::assert_integerish(max_n_subclust, lower = 2L)

      cdf <- self$get_clust_df(include_whereami = FALSE, include_subclusters = FALSE)
      big_clusters <- cdf %>% dplyr::filter(Volume >= !!break_nvox)

      if (nrow(big_clusters) > 0L) {
        private$pvt_subclust_list <- lapply(seq_len(nrow(big_clusters)), function(ii) {

          if (max_subclust_nvox > big_clusters$Volume[ii]) {
            this_max_nvox <- big_clusters$Volume[ii]
          } else {
            this_max_nvox <- max_subclust_nvox
          }

          roi_val <- big_clusters$roi_num[ii]

          # clone object to get same starting point on thresholds, NN, and the like
          sobj <- self$clone(deep = TRUE)
          sobj$reset_cache()

          # create temp mask file that just contains statistics within the big ROI of interest
          temp_mask <- tempfile(pattern = "tmpmask", fileext = ".nii.gz")
          temp_clust_1d <- tempfile(pattern = "tmpclust", fileext = ".1D")
          temp_clust_nii <- tempfile(pattern = "tmpclust", fileext = ".nii.gz")

          run_fsl_command(glue("fslmaths {private$pvt_pref_map} -uthr {roi_val} -thr {roi_val} -bin {temp_mask}"))

          # clusterize within mask
          sobj$mask <- temp_mask
          sobj$clusterize_output_file <- temp_clust_1d
          sobj$pref_map <- temp_clust_nii
          sobj$pref_dat <- NULL # not used for subclustering
          sobj$whereami <- NULL # clear out

          # run subclustering algorithm within this mask
          best <- sobj$run_subclustering(
            min_n_subclust, max_n_subclust, min_subclust_nvox, this_max_nvox, step_size,
            max_iter = max_iter, print_progress = print_progress
          )

          # if subclustering fails, best will be NULL
          if (!is.null(best)) {
            if (isTRUE(add_whereami) && best$has_clusters()) {
              best$add_whereami(whereami_atlases)
            }

            subcluster_df <- best$get_clust_df(include_subclusters = FALSE) %>%
              dplyr::rename(subroi_num = "roi_num") %>%
              dplyr::mutate(roi_num = !!roi_val) %>%
              dplyr::select("roi_num", "subroi_num", everything())

            # move mask and cluster table into same folder as input map
            dest_table <- glue("{file_sans_ext(private$pvt_input_file)}_roi{roi_val}_subclusters.1D")
            dest_map <- glue("{file_sans_ext(private$pvt_input_file)}_roi{roi_val}_subclusters.nii.gz")

            subcluster_files <- c(cluster_table = dest_table, cluster_map = dest_map)

            file.copy(temp_clust_1d, dest_table, overwrite = TRUE)
            file.copy(temp_clust_nii, dest_map, overwrite = TRUE)
            subcluster_files
          } else {
            subcluster_df <- NULL # no subclusters to return
            subcluster_files <- c(cluster_table = NA_character_, cluster_map = NA_character_)
          }

          ret_list <- list(roi_num = roi_val, subcluster_df = subcluster_df, subcluster_files = subcluster_files) # , subcluster_obj = best)
          private$pvt_subclust_files[[paste0("roi", roi_val)]] <- subcluster_files
          unlink(c(temp_mask, temp_clust_1d, temp_clust_nii)) # cleanup tmp files

          return(ret_list)
        })
      } else {
        message(glue("No clusters exceeded the maximum number of voxels ({break_nvox}) that would lead to subclustering."))
        return(invisible(self))
      }

      private$pvt_subclust_df <- dplyr::bind_rows(lapply(private$pvt_subclust_list, "[[", "subcluster_df"))
      return(invisible(self))
    },

    #' @description runs a subclustering algorithm on this object, increasing the thresholds until the desired constraints are satisfied
    #' @details this is intended to be used internally
    #' @param min_clust The minimum number of clusters that will be accepted
    #' @param max_clust The maximum number of clusters that will be accepted
    #' @param min_nvox The minimum number of voxels in a subcluster that will be accepted
    #' @param max_nvox The maximum number of voxels in a subcluster that will be accepted
    #' @param step_size The increments in the threshold values from one step to the next.
    #' @param max_iter The maximum number of increment steps that will be taken before giving up. Default: 50.
    #' @param refine_steps The number of steps backward from a winning solution. This maximizes the subcluster sizes. Default: 5.
    #' @param print_progress If TRUE, the user will see the thresholds being used to subcluster each region.
    run_subclustering = function(min_clust = NULL, max_clust = NULL, min_nvox = NULL, max_nvox = NULL, step_size = NULL,
    refine_steps = 5, max_iter = 50, print_progress = TRUE) {

      checkmate::assert_integerish(min_clust, lower = 2, len = 1L)
      checkmate::assert_integerish(min_nvox, lower = 1, len = 1L)
      checkmate::assert_number(step_size, lower = 1e-10)

      # get current threshold as starting values
      if (self$sided == "bi" || self$sided == "2") {
        vals <- c(lower_thresh = self$lower_thresh, upper_thresh = self$upper_thresh)
      } else {
        vals <- c(one_thresh = self$one_thresh)
      }

      delta_thresh <- function(vals, step_size, incr = 1) {
        step_size <- incr * step_size # determine the current step forward or backward
        if (length(vals) == 1L) {
          if (vals < 0) {
            return(vals - step_size) # make more negative
          } else {
            return(vals + step_size) # make more positive
          }
        } else if (length(vals) == 2L) {
          # make more extreme
          return(c(vals[1L] - step_size, vals[2L] + step_size))
        } else {
          stop("Cannot figure out vals")
        }
      }

      self$clust_nvox <- min_nvox # use subclustering settings for defining lower bound on clusters
      search_finished <- FALSE
      incr <- 0 # how much to walk forward or backward
      starting_vals <- vals # where did we start
      iter <- 0
      solution_found <- FALSE # pin a satisfactory solution while we look around
      refine_iter <- 0 # number of steps in refinement of solution
      best_vals <- NULL
      best_obj <- NULL # object of best 3dClusterize attempt

      while (search_finished == FALSE) {
        iter <- iter + 1

        # handle solution refinement steps (final search)
        if (isTRUE(solution_found)) {
          refine_iter <- refine_iter + 1
          incr <- -1 # always backup by the smaller step size (overrides incrs below from initial search)
          if (refine_iter > refine_steps) {
            search_finished <- TRUE
            break
          }
        }

        # get current threshold values based on increment and step size
        vals <- delta_thresh(vals, step_size, incr)
        if (incr < 0) {
          # if we have backed up all the way to the starting values, the search has failed
          if (abs(vals[length(vals)] - starting_vals[length(starting_vals)]) < 1e-5) {
            message("Could not backup further")
            search_finished <- TRUE
          }
        }

        if (self$sided == "one") {
          self$one_thresh <- vals[1]
        } else {
          self$lower_thresh <- vals[1]
          self$upper_thresh <- vals[2]
        }

        self$run(force = TRUE, quiet = TRUE)
        dd <- self$get_clust_df()
        n_subclust <- nrow(dd)

        if (n_subclust > 0L) {
          biggest_subclust <- max(dd$Volume)
          smallest_subclust <- min(dd$Volume)
        } else {
          biggest_subclust <- NA_integer_
          smallest_subclust <- NA_integer_
        }

        good_solution <- FALSE # whether the current values satisfy the constraints

        if (n_subclust == 0L) {
          # we have walked all the way out in the thresholds algorithm and have found no subclusters
          # consider this a search failure, give up
          search_finished <- TRUE
        } else if (n_subclust < min_clust) {
          # press on with higher thresholds
          incr <- 1
        } else if (biggest_subclust > max_nvox) {
          # press on with higher threshold
          incr <- 1
        } else if (n_subclust > max_clust) {
          # need to backup to get the max down
          incr <- -0.5
        } else if (smallest_subclust < min_nvox) {
          # need to backup to get larger subclusters
          incr <- -0.5
        } else {
          # a solution has been found
          # if it's the first solution, walk backwards in smaller steps
          if (isFALSE(solution_found)) {
            step_size <- step_size / (refine_steps + 1)
          }

          solution_found <- TRUE # these settings satisfy constraint -- refine, if relevant
          good_solution <- TRUE # this solution works (for printing progress)
          best_vals <- vals # on a solution refinement, if we get here, the values are good
          best_obj <- self$clone() # keep the best object in the backward search
          #message("Solution found")
        }

        if (isTRUE(print_progress)) {
          print(
            data.frame(t(vals), biggest = biggest_subclust, smallest = smallest_subclust, nclust = n_subclust, good_solution = good_solution),
            row.names = FALSE
          )
        }

        if (iter > max_iter && isFALSE(solution_found)) {
          message(glue::glue("Could not find solution after {max_iter} iterations. Perhaps adjust step size?"))
          search_finished <- TRUE
        }
      }

      if (isTRUE(solution_found)) {
        # was attempting to just update the object itself, but this does not overwrite the entire object with the cloned best one
        # self <- best_obj

        # need to re-run 3dClusterize once more at the best params to ensure that the output files reflect the final settings
        # the alternative would be to use new files with every iteration, which I may implement in future
        best_obj$run(force = TRUE, quiet = TRUE)
      } else {
        message("No subclusters found that satisfy constraints")
      }

      # instead of assigning self, we just return the best cloned object
      return(best_obj)
    },

    #' @description check whether clusteres were found
    #' @details returns TRUE if clusters were found, FALSE if they were not found, and NULL if the expected cluster
    #'   output file does not exist (e.g., if 3dClusterize has not been run yet)
    has_clusters = function() {
      private$pvt_has_clusters
    },

    #' @description not intended to be called by user, this resets the cluster data.frame and whereami objects to NULL
    #' @details this is used internally when cloning the parent clusterize object for subclustering
    reset_cache = function() {
      private$pvt_whereami <- NULL
      private$pvt_clust_df <- NULL
      private$pvt_subclust_df <- NULL
      private$pvt_subclust_list <- NULL
      return(invisible(self))
    },

    #' @description return list of subcluster details for each large ROI broken up by generate_subclusters()
    get_subclust_list = function() {
      private$pvt_subclust_list
    },

    #' @description Compares the clusters generated by 3dClusterize to an atlas of interest, then returns the subset of the
    #'   atlas that overlaps sufficiently with the map from 3dClusterize
    #' @param atlas_file A NIfTI file containing parcels or perhaps meta-analytic statistics
    #' @param atlas_lower_threshold Only retain values greater than this threshold in the comparison against the clusters. Default: 0
    #' @param atlas_upper_threshold Only retain values less than this threshold in the comparison against the clusters.
    #'   Default: Inf (retain all ROIs above the lower threshold)
    #' @param minimum_overlap The proportion overlap of an atlas parcel with a cluster required for the parcel to be retained.
    #'   Default: 0.8. If an integer > 1 is provided, then the function treats this as the minimum number of *voxels* that
    #'   must overlap between the parcel and the data-derived cluster.
    #' @param mask_by_overlap If TRUE, only voxels in the atlas that overlapped with a cluster are retained. In essence,
    #'   this erodes the retained atlas parcels to only include voxels that were in a cluster. Default: FALSE
    #' @param output_atlas the name/location of the file to output containing the subset of parcels in \code{atlas_file} that
    #'   are retained by this function. If \code{'default'} or \code{TRUE}, the subset atlas will be placed in the same
    #'   folder as the 3dClusterize input image, with a filename that combines the atlas file name with the input/threshold
    #'   image name. To disable creation of this file, set \code{output_atlas = FALSE}.
    #' @param roi_stats A character vector of summary stats to extract from the cluster-masked data file. Supported values:
    #'   \code{"mean"}, \code{"min"}, \code{"max"}. Set to NULL to skip extraction.
    subset_atlas_against_clusters = function(atlas_file = NULL, atlas_lower_threshold = 0, atlas_upper_threshold = Inf,
      minimum_overlap = 0.8, mask_by_overlap = FALSE, output_atlas = "default", roi_stats = c("mean", "max", "min")) {

      if (isFALSE(self$is_complete)) {
        message("Cannot run subset_atlas_against_clusters because 3dClusterize has not been run successfully yet. Use $run()")
        return(invisible(self))
      }

      if (checkmate::test_number(minimum_overlap, lower = 0.01, upper = 1.0)) {
        use_proportion <- TRUE
      } else if (checkmate::test_integerish(minimum_overlap, lower=2, len=1L)) {
        use_proportion <- FALSE
        minimum_overlap <- as.integer(minimum_overlap)
      } else {
        stop("minimum_overlap must be either a 0.01 -- 1.0 proportion or a > 1 number of voxels")
      }
      checkmate::assert_logical(mask_by_overlap, len=1L)
      clust_file <- self$get_outputs()$cluster_map
      if (is.na(clust_file)) {
        message("Cannot find a cluster file for the 3dClusterize object. Make sure you include a pref_map when you setup the object.")
        return(invisible(self))
      } else if (!checkmate::test_file_exists(clust_file)) {
        message(glue("Cannot find the cluster map {clust_file}"))
        return(invisible(self))
      }

      # figure out whether to output the subset atlas and if so, where
      # output_atlas becomes TRUE/FALSE and output_atlas_file is the location
      checkmate::assert_file_exists(atlas_file)
      if (checkmate::test_logical(output_atlas, len = 1L) && isTRUE(output_atlas)) {
        output_atlas_file <- "default"
      } else if (is.null(output_atlas)) {
        output_atlas <- FALSE
      } else if (checkmate::test_string(output_atlas)) {
        output_atlas_file <- output_atlas # name of output
        output_atlas <- TRUE # switch to logical
      }

      atlas_nii <- RNifti::readNifti(atlas_file)
      clust_nii <- RNifti::readNifti(clust_file)
      if (!identical(dim(atlas_nii), dim(clust_nii))) {
        cat("Atlas dimensions: ", paste(dim(atlas_nii), collapse=", "), "\n")
        cat("Cluster file dimensions: ", paste(dim(clust_nii), collapse = ", "), "\n")
        stop("Cannot proceed with subsetting because the dimensions are different!")
      }

      uvals <- sort(unique(as.vector(atlas_nii)))
      if (!checkmate::test_integerish(uvals)) {
        stop("At present, only integer-valued atlases are allowed.")
      } else {
        uvals <- uvals[uvals != 0] # never retain 0 as an atlas value
      }

      checkmate::assert_number(atlas_lower_threshold)
      checkmate::assert_number(atlas_upper_threshold)
      stopifnot(atlas_upper_threshold >= atlas_lower_threshold)
      if (!is.infinite(atlas_lower_threshold)) {
        atlas_nii[atlas_nii < atlas_lower_threshold] <- 0 # nullify voxels below threshold
      }

      if (!is.infinite(atlas_upper_threshold)) {
        atlas_nii[atlas_nii > atlas_upper_threshold] <- 0 # nullify voxels above threshold
      }

      meets_criteria <- rep(NA, length(uvals))
      prop_overlap <- rep(NA_real_, length(uvals))
      nvox_overlap <- rep(NA, length(uvals))
      nvox_parcel <- rep(NA, length(uvals))

      for (ii in seq_along(uvals)) {
        at <- 1L * (atlas_nii == uvals[ii]) # convert to 1/0 image
        clust_bin <- (1L * (clust_nii != 0L)) # convert to 1/0 image
        clust_match <- clust_bin * at
        nvox_parcel[ii] <- sum(at)
        nvox_overlap[ii] <- sum(clust_match)
        prop_overlap[ii] <- nvox_overlap[ii] / sum(at)

        if ((isTRUE(use_proportion) && prop_overlap[ii] >= minimum_overlap) ||
          (isFALSE(use_proportion) && nvox_overlap[ii] >= minimum_overlap)) {
          meets_criteria[ii] <- TRUE
        } else {
          meets_criteria[ii] <- FALSE
        }
      }

      good_vals <- uvals[meets_criteria]
      bad_vals <- uvals[!meets_criteria]

      if (length(good_vals) > 0L) {
        at_mod <- atlas_nii
        at_mod[!at_mod %in% good_vals] <- 0
        cat(glue("The following atlas values were retained: {paste(good_vals, collapse=', ')}"), "\n")
        cat(glue("The following atlas values were excluded: {paste(bad_vals, collapse=', ')}"), "\n\n")

        overlap_df <- data.frame(
          roi_val = uvals, prop_overlap = prop_overlap, retained = meets_criteria,
          nvox_overlap = nvox_overlap, nvox_parcel = nvox_parcel
        )

        if (isTRUE(mask_by_overlap)) {
          at_mod <- at_mod * clust_bin # mask out retained atlas voxels that did not overlap with a cluster
        }

        # add voxels retained (if masked by overlap)
        overlap_df$nvox_retained <- sapply(overlap_df$roi_val, function(x) sum(at_mod == x))

        extract_roi_stats <- FALSE
        stats_file <- self$get_outputs()$cluster_masked_data
        if (!is.null(roi_stats) && checkmate::test_file_exists(stats_file)) {
          stats_nii <- RNifti::readNifti(stats_file)

          # NA zero values in cluster-masked statistics since these will bias/invalidate stastistic
          # We want to get mean, min, max of voxels that survive in the cluster-masked data, not those masked out (0)
          stats_nii[abs(stats_nii) < 1e-5] <- NA
          extract_roi_stats <- TRUE
          roi_data <- list()

          for (ff in roi_stats) {
            if (ff == "mean") {
              roi_data[["mean"]] <- tapply(stats_nii, at_mod, mean, na.rm=TRUE)[-1L] # remove 0 values
            } else if (ff == "min") {
              roi_data[["min"]] <- tapply(stats_nii, at_mod, min, na.rm=TRUE)[-1L] # remove 0 values
            } else if (ff == "max") {
              roi_data[["max"]] <- tapply(stats_nii, at_mod, max, na.rm=TRUE)[-1L] # remove 0 values
            }
          }

          roi_data <- as.data.frame(roi_data)
          roi_data$roi_val <- good_vals

          overlap_df <- overlap_df %>% left_join(roi_data)
        }

        cat("Summary of overlap in atlas parcels and clusters\n")
        print(overlap_df, row.names = FALSE)

        if (isTRUE(output_atlas)) {
          if (output_atlas_file == "default") {
            atlas_name <- basename(file_sans_ext(atlas_file))
            output_atlas_file <- paste0(file_sans_ext(private$pvt_input_file), "_", atlas_name, "_overlap.nii.gz")
          }

          output_overlap_file <- paste0(file_sans_ext(output_atlas_file), ".csv")

          message(glue("Writing atlas subset image to file: {output_atlas_file}"))
          message(glue("Writing atlas subset data summary to file: {output_overlap_file}"))
          RNifti::writeNifti(image = at_mod, file = output_atlas_file)
          readr::write_csv(overlap_df, file=output_overlap_file)
          avec <- c(image = output_atlas_file, summary = output_overlap_file)
          attr(avec, "rois_retained") <- good_vals
          private$pvt_atlas_files[[atlas_name]] <- avec
        }
      } else {
        message("No atlas parcel overlapped sufficiently")
      }

      return(invisible(self))
    },

    #' @description method to delete any/all files generated by this object
    #' @param prompt if TRUE, user will have to confirm deletion of each file. If FALSE, files are deleted without prompting.
    delete_outputs = function(prompt = FALSE) {
      checkmate::assert_logical(prompt, len = 1L)
      f_list <- self$get_outputs()
      f_exists <- sapply(f_list, checkmate::test_file_exists)
      f_list <- f_list[f_exists]
      if (length(f_list) == 0L) {
        message("No output files found for removal.")
      } else {
        if (isTRUE(prompt)) {
          for (ff in f_list) {
            cat("File:", ff, "\n")
            remove <- askYesNo("Delete?")
            if (isTRUE(remove)) {
              cat(sprintf("Removing: %s", ff), sep = "\n")
              unlink(ff)
            } else if (is.na(remove)) {
              # cancel
              return(invisible(NULL))
            }
          }
        } else {
          cat(sprintf("Removing: %s", f_list), sep="\n")
          unlink(f_list)
        }

      }
      self$reset_cache()
      return(invisible(self))
    }
  )
)
