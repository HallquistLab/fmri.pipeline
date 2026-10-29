# This function pauses execution of an R script while a scheduled qsub job is not yet complete.

It is intended to give you control over job dependencies within R when
the formal PBS depend approach is insufficient, especially in the case
of a script that spawns child jobs that need to be scheduled or complete
before the parent script should continue.

## Usage

``` r
wait_for_job(
  job_ids,
  repolling_interval = 60,
  max_wait = 60 * 60 * 24,
  scheduler = "local",
  quiet = TRUE,
  stop_on_timeout = TRUE
)
```

## Arguments

- job_ids:

  One or more job ids of existing PBS or slurm jobs, or process ids of a
  local process for `scheduler="sh"`.

- repolling_interval:

  How often to recheck the job status, in seconds. Default: 30

- max_wait:

  How long to wait on the job before giving up, in seconds. Default: 24
  hours (86,400 seconds)

- scheduler:

  What scheduler is used for job execution. Options: c("torque", "qsub",
  "slurm", "sbatch", "sh", "local")

- quiet:

  If `TRUE`, `wait_for_job` will not print out any status updates on
  jobs. If `FALSE`, the function prints out status updates for each
  tracked job so that the user knows what's holding up progress.

- stop_on_timeout:

  Whether to stop if `max_wait` is reached.

## Value

Invisibly returns TRUE when all jobs complete successfully, or FALSE
when any job fails or is cancelled (or a non-stopping timeout occurs).

## Details

Note that for the `scheduler` argument, "torque" and "qsub" are the
same; "slurm" and "sbatch" are the same, and "sh" and "local" are the
same.

## Author

Michael Hallquist

## Examples

``` r
if (FALSE) { # \dontrun{
# example on qsub/torque cluster
wait_for_job("7968857.torque01.util.production.int.aci.ics.psu.edu", scheduler = "torque")

# example of waiting for two jobs on slurm cluster
wait_for_job(c("24147864", "24147876"), scheduler = "slurm")

# example of waiting for two jobs on local machine
wait_for_job(c("9843", "9844"), scheduler = "local")
} # }
```
