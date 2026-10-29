# helper function to submit a set of shell jobs that are independent of one another

helper function to submit a set of shell jobs that are independent of
one another

## Usage

``` r
cluster_submit_shell_jobs(
  job_list,
  commands_per_cpu = 1L,
  cpus_per_job = 8L,
  memgb_per_command = 8,
  time_per_job = "1:00:00",
  time_per_command = NULL,
  fork_jobs = TRUE,
  pre = NULL,
  post = NULL,
  sched_args = NULL,
  env_variables = NULL,
  wait_jobs = NULL,
  scheduler = "slurm",
  job_out_dir = getwd(),
  job_script_prefix = "job",
  log_file = "cluster_submit_jobs.csv",
  debug = FALSE
)
```

## Arguments

- job_list:

  a list or character vector where each element represents an
  independent job to execute in a shell environment

- commands_per_cpu:

  how many elements from `job_list` are executed by each core within a
  single job

- cpus_per_job:

  how many cpus/cores are requested for each job

- memgb_per_command:

  amount of memory (RAM) requested for each command (in GB)

- time_per_job:

  amount of time requested for each job

- time_per_command:

  amount of time requested for each individual command

- fork_jobs:

  if TRUE, all jobs within a single batch will run simultaneously using
  the fork (&) approach. Each task runs in its own subshell with \`set
  -e\`. The batch waits for every task and returns the first nonzero
  task status in submission order, or the cleanup status if all tasks
  succeed.

- pre:

  user-specified code to include in the job script prior to job_list
  (e.g., module load commands)

- post:

  user-specified cleanup code run after all tasks finish; its success
  cannot mask a task failure.

- sched_args:

  scheduler directives passed to
  [`cluster_job_submit()`](https://hallquistlab.github.io/fmri.pipeline/reference/cluster_job_submit.md).

- env_variables:

  named environment variables passed to
  [`cluster_job_submit()`](https://hallquistlab.github.io/fmri.pipeline/reference/cluster_job_submit.md).

- wait_jobs:

  jobs that must complete before submitted jobs run.

- scheduler:

  scheduler backend used for submission.

- job_out_dir:

  the directory where job scripts should be written

- job_script_prefix:

  the filename prefix for each job script

- log_file:

  a csv log file containing the job ids and commands that were executed

- debug:

  a logical indicating whether to actually submit the jobs (TRUE) or
  just create the scripts for inspection (FALSE)
