# R6 class representing a multivariate time series object for fMRI analysis

R6 class representing a multivariate time series object for fMRI
analysis

## Public fields

- `ts_data`:

  time x signals data.table

- `ts_keys`:

  RLE-encoded keying variables

- `event_data`:

  trial x event data, used for aligning time series with events

- `vm`:

  a list of variable mappings between internal constructs and input
  variable names

- `tr`:

  the repetition time (TR) of the fMRI sequence in seconds

## Methods

### Public methods

- [`fmri_ts$new()`](#method-fmri_ts-initialize)

- [`fmri_ts$get_ts()`](#method-fmri_ts-get_ts)

- [`fmri_ts$get_events()`](#method-fmri_ts-get_events)

- [`fmri_ts$add_keys()`](#method-fmri_ts-add_keys)

- [`fmri_ts$replace_vm()`](#method-fmri_ts-replace_vm)

- [`fmri_ts$get_kvars()`](#method-fmri_ts-get_kvars)

- [`fmri_ts$get_vmvec()`](#method-fmri_ts-get_vmvec)

- [`fmri_ts$export()`](#method-fmri_ts-export)

- [`fmri_ts$clone()`](#method-fmri_ts-clone)

------------------------------------------------------------------------

### `fmri_ts$new()`

Create a new fmri_ts object

#### Usage

    fmri_ts$new(ts_data = NULL, event_data = NULL, vm = NULL, tr = NULL)

#### Arguments

- `ts_data`:

  a data.frame or data.table containing time series

- `event_data`:

  a data.frame containing trial-level events that occurred in the time
  period represented by `ts_data`

- `vm`:

  a list of variable names used in `ts_data` and `event_data` that map
  onto internal constructs

- `tr`:

  the sampling rate (in seconds) of the fMRI data

------------------------------------------------------------------------

### `fmri_ts$get_ts()`

method to get rehydrated time series object with key values

#### Usage

    fmri_ts$get_ts(orig_names = FALSE)

#### Arguments

- `orig_names`:

  boolean indicating whether to return data.table with original naming
  scheme. Default: FALSE

------------------------------------------------------------------------

### `fmri_ts$get_events()`

Get a copy of the event data, optionally with original names

#### Usage

    fmri_ts$get_events(orig_names = FALSE)

#### Arguments

- `orig_names`:

  Whether to restore the input variable names

------------------------------------------------------------------------

### `fmri_ts$add_keys()`

method to add a variable in ts_data to the set of keying variables for
further use

#### Usage

    fmri_ts$add_keys(kv)

#### Arguments

- `kv`:

  A vector of one or more variables to RLE-encode and add as keys the
  object

------------------------------------------------------------------------

### `fmri_ts$replace_vm()`

method to replace one or more variable mappings in the object

#### Usage

    fmri_ts$replace_vm(...)

#### Arguments

- `...`:

  a set of arguments, each one of which replaces a field in the variable
  mapping with a new specification

------------------------------------------------------------------------

### `fmri_ts$get_kvars()`

return names of key variables

#### Usage

    fmri_ts$get_kvars()

------------------------------------------------------------------------

### `fmri_ts$get_vmvec()`

return variable mapping information

#### Usage

    fmri_ts$get_vmvec()

------------------------------------------------------------------------

### `fmri_ts$export()`

not currently used

#### Usage

    fmri_ts$export(filename)

#### Arguments

- `filename`:

  for writing out data

------------------------------------------------------------------------

### `fmri_ts$clone()`

The objects of this class are cloneable with this method.

#### Usage

    fmri_ts$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
