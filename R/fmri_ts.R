#essential structure
# - time data: time, signal, unit (voxel/region) ([voxels + time] x signals)
# - metadata: id, trial, events (trials x signals)
# - time -> unit mapping: vnum, atlas_value, region etc. (units x grouping)

#ts_data
# $time: time index
# $signal: only used if there is more than one value column (multivariate ts)
# $value: DV value of time series
# $key: 

#' R6 class representing a multivariate time series object for fMRI analysis
#'
#' @importFrom R6 R6Class
#' @export
fmri_ts <- R6::R6Class("fmri_ts",
  private=list(
    kvars=NULL, #internal names of key variables
    vvars=NULL, #internal names of value variables
    vmvec=NULL,  #mapping between input names and internal names
    update_variable_tracking=function(vm) {
      private$vmvec <- unlist(vm) #yields key1, key2, etc.
      if ("key" %in% names(vm)) {
        # if there is more than one key, unlist adds a number to each element of the vector
        # but if there is only one key, unlist keeps it as 'key'
        private$kvars <- if (length(vm$key) > 1L) paste0("key", seq_along(vm$key)) else "key"
      }
      if ("value" %in% names(vm)) {
        private$vvars <- if (length(vm$value) > 1L) paste0("value", seq_along(vm$value)) else "value"
      }
    },
    names_to_original=function(dt) {
      setnames(dt, names(private$vmvec), private$vmvec, skip_absent=TRUE)
    },
    names_to_internal=function(dt) {
      #look for naming collisions
      poss_conf <- private$vmvec[ intersect(names(private$vmvec), names(dt)) ]
      poss_conf <- poss_conf[poss_conf != names(poss_conf)]

      if (length(poss_conf) > 0L) {
        new_names <- if (length(poss_conf) > 1L) { paste0(".aux", seq_along(poss_conf)) } else { ".aux" }
        cat("To avoid naming conflicts, renaming these columns:", paste(names(poss_conf), collapse=", "),
          "to:", paste0(new_names, collapse=", "), "\n")
        self$vm[[".aux"]] <- names(poss_conf)
        setnames(dt, self$vm[[".aux"]], new_names, skip_absent=FALSE)
        private$vmvec <- unlist(self$vm) #yields key1, key2, etc.
      }

      #handle internal renaming to make programming with these objects easy
      setnames(dt, private$vmvec, names(private$vmvec), skip_absent=TRUE)
    }

  ),
  public=list(
    #' @field ts_data time x signals data.table
    ts_data=NULL,

    #' @field ts_keys RLE-encoded keying variables
    ts_keys=NULL,

    #' @field event_data trial x event data, used for aligning time series with events
    event_data=NULL,

    #' @field vm a list of variable mappings between internal constructs and input variable names
    vm=NULL,

    #' @field tr the repetition time (TR) of the fMRI sequence in seconds
    tr=NULL,

    #' @description Create a new fmri_ts object
    #' @param ts_data a data.frame or data.table containing time series
    #' @param event_data a data.frame containing trial-level events that occurred in the time period represented by \code{ts_data}
    #' @param vm a list of variable names used in \code{ts_data} and \code{event_data} that map onto internal constructs
    #' @param tr the sampling rate (in seconds) of the fMRI data
    initialize = function(ts_data=NULL, event_data=NULL, vm=NULL, tr=NULL) {
      if (is.null(tr)) { stop("tr must be provided at object initialization") }
      checkmate::assert_numeric(tr, lower=1e-2, null.ok=FALSE)

      default_vm <- list(id="id", run="run", trial="trial", run_trial="run_trial", time="time", value="value")
      if ("key" %in% names(ts_data)) { default_vm[["key"]] <- "key" } #add default key mapping if keyed
      for (nn in names(default_vm)) { #populate default variable mappings if not provided in input vector
        if (!nn %in% names(vm)) { vm[nn] <- default_vm[nn] }
      }

      checkmate::assert_data_frame(ts_data)
      checkmate::assert_data_frame(event_data, null.ok=TRUE)

      #setup standardized naming
      private$update_variable_tracking(vm)

      if (!is.null(event_data)) {
        if (!is.null(vm$id)) {
          checkmate::assert_string(vm$id)
          stopifnot(vm$id %in% names(event_data))
          if (length(unique(event_data[[vm$id]])) > 1L) {
            stop("fmri_ts objects only support single runs of data for single IDs.",
            " You can combine fmri_ts objects using combine_ts().")
          }
        }
        if (!is.null(vm$run)) {
          checkmate::assert_string(vm$run)
          stopifnot(vm$run %in% names(event_data))
          if (length(unique(event_data[[vm$run]])) > 1L) {
            stop("fmri_ts objects only support single runs of data for single IDs.",
            " You can combine fmri_ts objects using combine_ts().")
          }
        }

        checkmate::assert_string(vm$trial) #singular string
        stopifnot(vm$trial %in% names(event_data))

        event_data <- data.table(event_data)

        #handle internal renaming to make programming with these objects easy
        private$names_to_internal(event_data)

        data.table::setorderv(event_data, "trial") #event names are now internal
      }

      # verify presence of required columns
      sapply(vm[c("time", "value", "key")], function(x) { stopifnot(all(x %in% names(ts_data))) })

      # convert to data table for speed, memory management
      ts_data <- data.table(ts_data)

      #handle internal renaming to make programming with these objects easy
      private$names_to_internal(ts_data)

      if (!is.null(private$kvars)) {
        setorderv(ts_data, private$kvars)
        ts_keys <- sapply(private$kvars, function(x) {
          if (is.factor(ts_data[[x]])) ts_data[[x]] <- as.character(ts_data[[x]]) # have to convert to character prior to rle
          rle(ts_data[[x]])
        }, simplify = FALSE)
        ts_data[, private$kvars := NULL] # drop key columns
        self$ts_keys <- ts_keys
      }

      self$ts_data <- ts_data
      self$event_data <- event_data
      self$vm <- vm
      self$tr <- tr
    },

    #' @description method to get rehydrated time series object with key values
    #' @param orig_names boolean indicating whether to return data.table with original naming scheme. Default: FALSE
    get_ts = function(orig_names=FALSE) {
      tsd <- data.table::copy(self$ts_data) #ensure that we copy the object to avoid altering $ts_data
      for (kk in seq_along(self$ts_keys)) { tsd[, names(self$ts_keys)[kk] := inverse.rle(self$ts_keys[[kk]])] }
      data.table::setcolorder(tsd, private$kvars) #put keying variables first in object
      if (isTRUE(orig_names)) { private$names_to_original(tsd) }
      return(tsd)
    },

    #' @description Get a copy of the event data, optionally with original names
    #' @param orig_names Whether to restore the input variable names
    get_events = function(orig_names=FALSE) {
      events <- data.table::copy(self$event_data)
      if (isTRUE(orig_names) && !is.null(events)) private$names_to_original(events)
      events
    },

    #' @description method to add a variable in ts_data to the set of keying variables for further use
    #' @param kv A vector of one or more variables to RLE-encode and add as keys the object
    add_keys = function(kv) {
      checkmate::assert_character(kv)
      internal <- names(private$vmvec)[match(kv, private$vmvec)]
      internal[is.na(internal)] <- kv[is.na(internal)]
      stopifnot(all(internal %in% names(self$ts_data)))
      nkeys <- length(private$kvars) #number of existing keys
      for (vv in seq_along(kv)) {
        vname <- internal[vv]
        # convert factors to character before rle (loses levels/detail)
        if (is.factor(self$ts_data[[vname]])) self$ts_data[[vname]] <- as.character(self$ts_data[[vname]])
        self$ts_keys[[nkeys+vv]] <- rle(self$ts_data[[vname]]) # RLE-compress the new key
      }

      self$vm[["key"]] <- c(self$vm[["key"]], kv)
      private$update_variable_tracking(self$vm) #refresh mapping to internal names
      names(self$ts_keys) <- private$kvars #name ts_keys according to key variables to let get_ts() work

      self$ts_data[, (internal) := NULL] #drop new keys from rectangular data
    },

    #' @description method to replace one or more variable mappings in the object
    #' @param ... a set of arguments, each one of which replaces a field in
    #'    the variable mapping with a new specification
    #  TODO: create a private validate_vm method and run both initialize and replace through it
    replace_vm = function(...) {
      replist <- list(...)
      repfields <- names(replist)
      stopifnot(all(repfields %in% names(self$vm))) #must be replacement

      #revert to original names before we modify
      private$names_to_original(self$ts_data)

      for (rr in seq_along(replist)) {
        this_field <- repfields[rr]
        self$vm[[this_field]] <- replist[[rr]]
      }

      #update internal renaming scheme
      private$vmvec <- unlist(self$vm) #yields key1, key2, etc.

      #convert to internal naming scheme
      private$names_to_internal(self$ts_data)
    },

    #' @description return names of key variables
    get_kvars = function() { #simple get method to allow access to key variables
      return(private$kvars)
    },

    #' @description return variable mapping information
    get_vmvec = function() {
      return(private$vmvec)
    },

    #' @description not currently used
    #' @param filename for writing out data
    export = function(filename) {

    }
  )
)
