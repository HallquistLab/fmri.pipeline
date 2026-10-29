# This function submits a single script to a high-performance cluster using a scheduler (Slurm or TORQUE). It accepts a vector of arguments to be passed to the scheduler and a vector of environment variables that should be passed to the compute node at job execution.

The function returns the jobid of the scheduled job.

## Usage

``` r
cluster_job_submit(
  script,
  scheduler = "slurm",
  sched_args = NULL,
  env_variables = NULL,
  export_all = FALSE,
  echo = TRUE,
  fail_on_error = FALSE,
  wait_jobs = NULL,
  wait_signal = "afterok",
  repolling_interval = 60,
  tracking_sqlite_db = NULL,
  tracking_args = list()
)
```

## Arguments

- script:

  A script that should be executed by the scheduler. This can contain
  scheduler directives, but in the case of conflicts, the directives
  passed with `sched_args` will take precedence.

- scheduler:

  Which scheduler to use for job submission. Options are 'qsub',
  'torque', 'sbatch', 'slurm', or 'sh'. The terms 'qsub' and 'torque'
  are aliases (where 'torque' submits via the qsub command). Likewise
  for 'sbatch' and 'slurm'. The scheduler 'sh' does not submit to any
  scheduler at all, but instead executes the command immediately via sh.

- sched_args:

  A character vector of arguments to be included in the scheduling
  command. On TORQUE, these will typically begin with '-l' such as '-l
  wall_time=10:00:00'.

- env_variables:

  A named character vector containing environment variables and their
  values to be passed to the `script` at execution time. This is handled
  by the '-v' directive on TORQUE clusters and by '–export' on Slurm
  clusters. The names of this vector are the environment variable names
  and the values of the vector are the environment variable values to be
  passed in. If you want to propagate the current value of an
  environment variable to the compute node at runtime, use NA as the
  value of the element in `env_variables`. See examples.

- export_all:

  Whether to export all environment variables to the compute node at
  runtime. Default: FALSE

- echo:

  Whether to echo the job submission command to the terminal at the time
  it is scheduled. Default: TRUE.

- fail_on_error:

  Whether to stop execution of the script (TRUE), or issue a warning
  (FALSE) if the job submission fails. Defaults to FALSE (i.e., issue a
  warning).

- wait_jobs:

  a character string of jobs or process ids that should complete before
  this job is executed

- wait_signal:

  on torque or slurm clusters, the signal that should indicate that
  parent jobs have finished.

- repolling_interval:

  The number of seconds to wait before rechecking job status (used only
  for local scheduler)

- tracking_sqlite_db:

  Optional SQLite database used for job tracking.

- tracking_args:

  Optional job-tracking metadata passed through to the tracking
  database.

## Value

A character string containing the jobid of the scheduled job, or NULL
after a failed submission when fail_on_error is FALSE.

## Author

Michael Hallquist

## Examples

``` r
if (FALSE) { # \dontrun{
  #simple PBS submission
  cluster_job_submit(
    "myscript.bash",
    scheduler = "torque",
    sched_args = c("-l walltime=10:00:00", "-l nodes=1:ppn=20"),
    env_variables = c(RUN_INDEX = 2, MODEL_NAME = "FSE21")
  )

  #To forward environment variables without explicitly providing values. Note that these must
  #  be in R's system environment (cf. Sys.getenv) at execution time to forward as expected.
  cluster_job_submit('myscript.sbatch', scheduler="slurm",
     sched_args=c('-p general', '-N 1', '-n 12', '--mem=10g', '-t 02-00:00:00'),
     env_variables=c(RUN_INDEX=2, R_HOME=NA, JAVA_HOME=NA))
} # }
```
