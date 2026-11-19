# Build a task fMRI design from events and signals

Builds run-level timing and convolved regressors used by the FSL, AFNI,
and SPM model writers. FSL and AFNI timing-file exports are also
available. The SPM writer reads events and modulators directly from the
returned design object.

## Usage

``` r
build_design_matrix(
  events = NULL,
  signals = NULL,
  tr = NULL,
  center_values = TRUE,
  hrf_parameters = c(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35),
  baseline_coef_order = -1L,
  baseline_parameterization = "Legendre",
  run_data = NULL,
  run_nifti_drop_applied = TRUE,
  drop_volumes = 0L,
  runs_to_output = NULL,
  plot = TRUE,
  write_timing_files = NULL,
  output_directory = "run_timing",
  keep_empty_regressors = TRUE,
  convolve_wi_run = TRUE,
  high_pass = NULL,
  iti_post = 12,
  ts_multipliers = NULL,
  additional_regressors = NULL,
  lg = NULL
)
```

## Arguments

- events:

  a data.frame that includes a column for the event type (e.g. outcome
  vs. cue), run number (1:n), trial number (nested within run, 1:x),
  onset, duration of event

- signals:

  expects a list of list. The first level of lists reflects each signal
  (e.g. pe vs values). In the second level, BDM expects values and event
  (e.g. cue vs outcome). Values is a `data.frame` with run number (1:n),
  trial (1:x), and signal.

- tr:

  The repetition time of your fMRI sequence in seconds. You must specify
  this. This is important to specify correctly to get temporal filtering
  correct.

- center_values:

  A logical (`TRUE/FALSE`) indicating whether to center parameteric
  regressors prior to convolution. This is usually a good idea to remove
  collinearity between parametric and task indicator regressors. The
  default is `TRUE`.

- hrf_parameters:

  A named vector of HRF parameters passed to `fmri.stimulus` internally.
  The parameters are a1, a2, b1, b2, cc. Equation is (x/d1)^a1 \*
  exp(-(x - d1)/b1) - c \* (x/d2)^a2 \* exp(-(x - d2)/b2). Defaults and
  descriptions are: a1 = 6. Controls the time delay to the peak of the
  positive response a2 = 12. Controls the time delay to the (negative)
  peak of the undershoot b1 = 0.9. Controls the dispersion (breadth) of
  the positive response b2 = 0.9. Controls the dispersion (breadth) of
  the undershoot cc = 0.35. Controls the relative scaling of the
  positive response versus the undershoot.

  Note: These defaults come from Glover 1999.

  Note. FSL double gamma has a1 = 6, a2 = 16, cc = 1/6. Not yet sure
  about b1 and b2.

- baseline_coef_order:

  Default -1 (no baseline). If \>= 0, then design will include
  polynomial trends within each run (e.g. baseline_coef_order = 1
  includes both an intercept and a linear trend as regressors)

- baseline_parameterization:

  Defaults to "Legendre". This adds Legendre polynomials up to
  `baseline_coef_order` (e.g., 2). The alternative is
  "orthogonal_polynomials", which uses `fmri.design` from the `fmri`
  package (now internal to fmri.pipeline) to add polynomial regressors
  that are orthogonal to substantive design factors.

- run_data:

  a data.frame containing metadata about the runs for which we want to
  model the task design. This data.frame should contain the columns
  run_number, run_volumes, run_nifti, and drop_volumes. If run_nifti is
  provided, but run_volumes is not, then the number of volumes is looked
  up from the NIfTI header.

- run_nifti_drop_applied:

  Logical indicating whether the NIfTI inputs in `run_data$run_nifti`
  already have `drop_volumes` removed. If `TRUE` (default), reported
  NIfTI lengths are treated as post-drop lengths; if `FALSE`,
  `drop_volumes` are subtracted when deriving modeled volumes.

- drop_volumes:

  By default, all volumes are retained. If specified, this can be a
  vector of the number of volumes that will be removed from the
  *beginning* of each convolved regressor. If you pass a single number
  (e.g., 3), this number of volumes will be dropped from each run. This
  is useful if you have dropped the first n volumes of your MR data, for
  example to handle problems with steady state magnetization. If you
  include drop_volumes in the `run_data` data.frame, this will be used
  over this setting.

- runs_to_output:

  A numeric vector of run numbers to be output. By default, all runs are
  preserved. This is used to model only a subset such as `c(1, 2, 6)`.

- plot:

  By default (`TRUE`), `build_design_matrix` will plot the design matrix
  in the plot window of your R session. If `FALSE`, the plot is not
  displayed, but the ggplot object is still provided in the
  \$design_plot field.

- write_timing_files:

  When NULL (the default), the function does not write timing files to
  disk. This argument accepts a character vector that specifies whether
  to write "AFNI", "FSL", or "convolved" timing files to disk (in
  `output_directory`). AFNI files follow the dmBLOCK convention of
  TIME\*PARAMETER:DURATION for each event. FSL files follow the
  three-column format of onset, duration, value. And convolved files
  represent a given signal convolved with the HRF such that the 1-column
  output is in volumes (i.e., one row per volume).

- output_directory:

  Where to output the timing files. By default, the function will output
  the timing files to a folder called "run_timing" in the current
  working directory. If such a folder does not exist, it will make a
  folder in your R session's current working directory.

- keep_empty_regressors:

  If TRUE, then a regressor that contains no events for a given run is
  still retained as an all-zeros regressor in the convolved timing files
  and the FSL 3-column files. This is useful if you want to include
  regressors that only occur for some subjects, but not others. That
  said, handling contrasts involving these regressors is a more
  complicated matter.

- convolve_wi_run:

  If `TRUE` (the default), convolution – and any corresponding mean
  centering and normalization of the heights of parametric signals – is
  applied to each run separately. If FALSE, the events across runs are
  concatenated before convolution is applied (i.e., treating it as one
  long time series).

- high_pass:

  By default, `NULL`. If desired, pass in a number in Hz that specifies
  the high pass filter cutoff. In this case a FIR-based high-pass filter
  will be applied to remove any low-frequency fluctuations in the design
  signals after convolution. This can be useful and necessary if the MRI
  have been filtered, but the regressors have not. It is important that
  the frequency content of both MR and design signals matches. Some
  programs, including FEAT, ensure that equivalent filtering is applied
  to bot the Y and X sides of this equation, but you should be clear
  whether your program does so, too. If it doesn't, probably best to use
  this argument to filter things yourself. For example, 3dDeconvolve
  only handles filtering via a set of drift (polort) regressors. If you
  have used another tools such as fslmaths -bptf to filter the data, the
  polort will not necessarily result in the same removal of drift from
  regressors as was applied to the MR data.

- iti_post:

  By default, 12. Assumes 12 volumes after each run. Only necessary to
  specify if not supplying run_volumes and expecting function to use
  events information to calculate run_volumes. Wouldn't recommend this,
  just a default here.

- ts_multipliers:

  By default, `NULL`. If specified, expects either a vector of character
  strings for different .txt files for the time-based regressors OR a
  list of data.frames (1 df per run). These data.frames contain time
  series regressors (i.e., having the same length as run_volumes) and
  can be used as multipliers on a stimulus signal prior to convolution.
  This is primarily intended for PPI analysis, where the stimulus
  regressor is multiplied by a time series from a seed/candidate region.
  To use columns of `ts_multiplier`, you include a specifier for an
  element in the signals list: `ts_multiplier="vmPFC"` where "vmPFC" is
  column name in `ts_multipliers`. If you use a list of character
  vectors (i.e., read from .txt files), make sure that you have headers
  in the text files that will become the names of the ts_multiplier
  signals.

- additional_regressors:

  By default, `NULL`. If additional regressors specified, either expects
  character vector for different .txt files for the additional
  regressors (1 txt file per run) OR it expects a list of data.frames (1
  df per run). These values are tacked onto design_convolved (and not
  convolved with HRF), so each regressor should be length of the number
  of run_volumes within that run. If you pass in a vector of .txt files
  containing additional regressors, these will be read into R, truncated
  to run_volumes, and column-wise concatenated with substantive
  regressors.

- lg:

  An lgr object used for logging messages

## Value

A list of containing different aspects of the design matrix:

- `$design`: A runs x signals list containing events before convolution.
  Each element is a 2-D matrix containing, minimally, "trial", onset",
  "duration", and "value" columns. Onsets and durations are specified in
  seconds, consistent with FSL's 3-column format. Within each matrix,
  the onset, duration and value of the signal is specified.

- `$design_convolved`: The convolved design matrices for each run. Each
  element in the list contains a run. Within each design matrix, each
  column contains a regressor, encompassing substantive regressors,
  additional signals, and polynomial baseline regressors, if specified.
  Each row reflects the predicted value for each volume.

- `$design_unconvolved`: The unconvolved design matrices for each run.
  Same structure as `$design_convolved`, but prior to convolution with
  the HRF.

- `$collin_events`: A list containing information about the collinearity
  of regressors before convolution. At the highest level of the list,
  each element contains a run. At the second level of the list, the
  first element contains the correlation matrix of the regressors and
  the second element provides the variance inflation factor (VIF)
  associated with each regressor. Example:
  `design$collin_events$run1$vif`

- `$collin_convolved`: A list containing information about collinearity
  of the convolved regressors, including substantive signals, additional
  regressors, and polynomial regressors. Follows the same structure as
  `$collin_events`.

- `$concat_onsets`: A list containing concatenated event onset times for
  each signal. Each signal is an element of the list containing a vector
  of all onset times across runs. That is, the total time of run1 is
  added to onsets of run2 to support a combined analysis of all runs,
  which is common in AFNI (e.g., using -concat or multiple files to
  -input).

- `$run_volumes`: A vector containing the total number of volumes
  modeled for each run.

- `$design_plot`: A ggplot object showing the design matrix. This is
  generated by `visualize_design_matrix`.

## Details

The function outputs a list containing key aspects of the fMRI design,
including the unconvolved and convolved regressors, collinearity
diagnostics, the number of volumes modeled in each run, and a plot of
the design.

The basic logic of the inputs to build_design_matrix is that
task-related fMRI designs are organized around a set of events that
occur in time and have a specific duration. Furthermore, for a given
event, it could be a 0/1 non-occurrence versus occurrence
representation, *or* the event could be associated with a specific
parametric value such as working memory load, reward prediction error,
or expected value. These parametric effects are aligned in time with an
event, but there may be multiple predictions for a given event. For
example, we may align a 0/1 regressor and a reward prediction error the
outcome phase of a task.

Thus, the function abstracts timing-related information into `events`
and signals, whether parametric or binary, into the `signals`.

The `events` argument expects a `data.frame` that has, minimally, the
following structure:


     > print(events)
         event run_number trial onset duration
           cue          1     1     4        2
           cue          1     2     7        2
       outcome          1     1     6      0.5
       outcome          1     2   9.5      0.5
           cue          2     1   1.2        2
           cue          2     2    12        2
       outcome          2     1     6      0.5
       outcome          2     2   9.5      0.5

Note that you can tack on other columns to `events` if it useful to you.
Furthermore, if you want to test different durations (e.g., RT-convolved
versus fixed duration versus instantaneous), you can add these as
additional columns (e.g., `duration_1s`, `duration_instant`, etc.). To
make use of these in the design, specify the column name in events in
the `$duration` element of a given signal in the `signals` list. If you
do not specify the `$duration` element in `signals`,
`build_design_matrix` will assume that the relevant duration is stored
in the `$duration` column of `events`.

The `signals` argument expects a list where each element is a given
signal that should be aligned with an event and that has some height
(e.g., 0/1 or a parametric value) prior to convolution. The signals list
should be named by signal and each element should be a list itself, such
as the following:


      signals <- list(
        cue=list(event="cue", duration=0, value=1, normalization="none")
      )

The `event` element specifies the mapping between a given signal and the
corresponding timing in the `events` `data.frame`. In essence, this is
used to merge the event and signal data together. Here, we specify that
the cue signal is aligned in time with the cue event.

The `duration` element can be:

1.  A single number, in which case this fixed duration is used for all
    events

2.  A name of the column to be used in the `events` `data.frame` (e.g.,
    "duration_rtshift")

3.  Omitted altogether, in which case `build_design_matrix` will default
    to the "duration" column of `events`.

The `value` element can be a single number (e.g., 1) in which case this
height is used for all corresponding occurrences of a given event. Most
commonly, a fixed value is useful for modeling a 'taskness' regressor,
which captures a 0/1 representation of whether an event is occurring at
a given moment in time. In conventional task-reated fMRI, this task
indicator representation is then convolved with the HRF to model
expected BOLD activity due to the occurrence of an event. Alternatively,
`value` can be a data.frame containing `$run_number`, `$trial`, and
`$value` columns that specify the height of the regressor at each trial.
This specification is more useful for a parametric regressor, as in
model-based fMRI.

Optionally, one or more within-subject factor columns can be included in
the value data.frame and noted in the `wi_factors` element of the
signal. In this case, regressors for each level of the wi_factors will
be generated as separate regressors.

Here is an example:


      signals <- list(
        pe=list(event="outcome", normalization="none", convmax_1=TRUE, wi_factors="trustee",
        value=data.frame(
          run_number=rep(1,5),
          trial=1:5,
          value=c(10.2, -11.1, 6, 2.4, 1.5),
          trustee=c("Good", "Good", "Bad", "Bad", "Neutral")
        )
      )

Here, the parametrically varying prediction error signal will be aligned
at the "outcome" event, have a duration copied from the `$duration`
column of `events`, and will have parametrically varying heights (e.g.,
10.2 at trial 1) prior to convolution. Note that the value `data.frame`
need not have an entry for every trial in the run. For example, if a
given signal is only relevant or only occurs for some "outcome" events,
the trial column might be something like `c(2, 6, 10)`, indicating that
the parametric modulator is only modeled at those trials. This is
achieved by joining `events` with the relevant signal using `trial` as a
key.

The `$normalization` element handles the normalization of the HRF for
each regressor. This can be:

1.  `durmax_1`: pre-convolution, normalize the HRF max to 1.0 for long
    events (15+ sec) such that height of HRF is modulated by duration of
    event but maxes at 1. This is identical to dmUBLOCK(0).

2.  `evtmax_1`: pre-convolution, normalize the HRF max to 1.0 for each
    stimulus regardless of duration. This is identical to dmUBLOCK(1).

3.  `none`: No normalization of the HRF is performed prior to
    convolution.

The optional `$convmax_1` element handles rescaling the *convolved*
regressor to a maximum height of 1.0. If TRUE for a given signal, the
convolved regressor will be divided by its max, leading to a max of 1.0
across both runs (assuming `convolve_wi_run` is `TRUE`) and subjects.
This may be useful for scaling the regression coefficients in voxelwise
regression across subjects. For example, if the parametric signal
captures similar dynamics within subjects over the experiment, but the
scaling varies substantially between subjects, `convmax_1` can help to
place the betas on an equivalent scale across subjects (assuming the MR
data are also scaled similarly between subjects).

The optional `$demean_convolved` element handles whether to demean a
convolved regressor.

The optional `$beta_series` element handles whether to convert a given
signal into a set of regressors, one per event. This results in a wide
design matrix in which the beta series regressors each reflect a single
convolved event. All other signal arguments apply such as HRF
normalization. `$beta_series` defaults to FALSE.

The optional `ts_multipliers` element can be added to a given signal in
order to multiply the event in the design by a continuous time series
such as an ROI. This is primarily used to compute interaction terms
between events in the design and activity in a given region in order to
examine connectivity using a psychophysiological interaction (PPI)
approach. The `build_design_matrix` function will mean center the time
series prior to multiplying it by the relevant design regressor, then
convolve the result with the HRF. Thus, it is typical to provide a
\*deconvolved\* time series to `ts_multipliers`, though some packages
(e.g., FSL) don't typically use deconvolution in this way.

Finally, the optional `add_deriv` element determines whether the
temporal derivative of a regressor is added to the design matrix after
convolution. Following FSL, the derivatives are computed by a
first-order difference and are then residualized for other regressors in
the matrix. That is, the derivatives are orthogonalized with respect to
substantive regressors. By default, derivatives are not added, but if
`TRUE` for a given signal, this will be added to the convolved design
matrix.

This function was adapted from the fitclock package
(https://github.com/PennStateDEPENdLab/fitclock.git) to allow for more
general creation of design matrices for fMRI analyses.

## Author

Michael Hallquist

Alison Schreiber

## Examples

``` r

if (FALSE) { # \dontrun{
  data(example_events)
  data(example_signals)

  #basic convolved design matrix
  d <- build_design_matrix(
    events = example_events,
    signals = example_signals,
    tr = 1.0,
    plot = FALSE
  )

  data(example_nuisrun1) #load demo additional signals
  data(example_nuisrun2)

  #design matrix with 0,1,2 polynomial baseline regressors and a set of additional regressors
  #this does not contain a 'taskness' regressor for cue or outcome
  dnuis <- build_design_matrix(
    events = example_events,
    signals = example_signals,
    tr = 1.0,
    additional_regressors = list(example_nuisrun1, example_nuisrun2),
    baseline_coef_order = 2
  )

  #tweak the design to add temporal derivatives for both ev and pe
  example_signals$pe$add_deriv <- TRUE
  example_signals$ev$add_deriv <- TRUE

  #also use the evtmax_1 normalization method for both parametric signals
  example_signals$pe$normalization <- "evtmax_1"
  example_signals$ev$normalization <- "durmax_1"

  #finally, add a taskness regressor for cue and outcome
  example_signals$cue_evt <- list(value=1, event="cue", normalization="none")
  example_signals$outcome_evt <- list(value=1, event="outcome", normalization="none")

  #include up to quadratic drift terms in the design, drop 3 volumes from the beginning,
  #and write timing files in AFNI, FSL, and convolved formats (to the "run_timing" directory).
  d_modified <- build_design_matrix(
    events = example_events,
    signals = example_signals,
    tr = 1.0,
    baseline_coef_order = 2,
    drop_volumes = 3,
    write_timing_files = c("convolved", "AFNI", "FSL")
  )

  #show unconvolved design
  plot(visualize_design_matrix(concat_design_runs(d_modified, convolved=FALSE)))

  #beta series approach for cues: 1 regressor per cue event
  example_signals$cue_evt$beta_series <- TRUE
  example_signals$ev <- NULL
  example_signals$pe <- NULL

  d_beta <- build_design_matrix(
    events = example_events,
    signals = example_signals,
    tr = 1.0,
    plot = FALSE,
    baseline_coef_order = 2,
    drop_volumes = 3,
    write_timing_files = c("convolved", "FSL")
  )

  #PPI example

  events <- rbind(
    data.frame(event="cue", run_number=1, trial=1:5, onset=c(5, 20, 50, 100, 150), duration=4),
    data.frame(event="cue", run_number=2, trial=1:5, onset=c(15, 35, 60, 90, 105), duration=4)
  )

  signals <- list(
    load=list(value=rbind(
      data.frame(run_number=1, trial=1:5, value=1:5),
      data.frame(run_number=2, trial=1:5, value=1:5)
    ), event="cue", normalization="none"),
    cue_evt=list(value=1, event="cue", normalization="none"),
    cue_evt_ppi=list(value=1, event="cue", normalization="none", ts_multipliers="vs")
  )

  library(neuRosim)
  forppi1 <- simTSrestingstate(nscan=328, base=1, TR=1, SNR=6)
  forppi2 <- simTSrestingstate(nscan=242, base=1, TR=1, SNR=6)
  test_ppi <- list(run1=data.frame(vs=forppi1), run2=data.frame(vs=forppi2))

  #In PPI, the time series itself should typically be included as a regressor
  d <- build_design_matrix(events = events, signals = signals, tr=0.5, center_values=TRUE,
                write_timing_files = c("convolved", "AFNI", "FSL"),
                drop_volumes = 0, baseline_coef_order = 2, run_volumes = c(328, 242),
                ts_multipliers = test_ppi, additional_regressors = test_ppi)
} # }
```
