# Execute a compiled FWE plan

Execute ready tasks locally or submit them through Slurm or TORQUE.
Completed tasks are skipped unless `force = TRUE`. Local execution is
synchronous and returns a refreshed plan; scheduler execution returns
submitted job IDs and can later be refreshed with
[`refresh_fwe_plan()`](https://hallquistlab.github.io/fmri.pipeline/reference/refresh_fwe_plan.md).

## Usage

``` r
run_fwe_plan(
  plan,
  scheduler = NULL,
  force = FALSE,
  dry_run = FALSE,
  allow_blocked = FALSE,
  worker_script = NULL,
  commands_per_cpu = NULL,
  cpus_per_job = NULL,
  memgb_per_command = NULL,
  time_per_command = NULL,
  sched_args = NULL,
  env_variables = NULL,
  wait_jobs = NULL,
  stop_on_error = TRUE,
  persist_results = TRUE,
  tracking_sqlite_db = NULL,
  tracking_args = list()
)
```

## Arguments

- plan:

  an `fwe_plan` or path to an RDS plan.

- scheduler:

  execution backend. Defaults to the scheduler compiled into the plan.

- force:

  rerun tasks whose required artifacts already exist.

- dry_run:

  compile commands without creating directories, executing, or
  submitting jobs.

- allow_blocked:

  skip blocked tasks instead of rejecting the run.

- worker_script:

  optional method worker override. Primarily useful for development and
  testing.

- commands_per_cpu, cpus_per_job, memgb_per_command, time_per_command:

  resource settings passed to
  [`cluster_submit_shell_jobs()`](https://hallquistlab.github.io/fmri.pipeline/reference/cluster_submit_shell_jobs.md)
  for Slurm and TORQUE execution.

- sched_args, env_variables, wait_jobs:

  optional scheduler settings passed to
  [`cluster_submit_shell_jobs()`](https://hallquistlab.github.io/fmri.pipeline/reference/cluster_submit_shell_jobs.md).
  Named `env_variables` are also applied during local execution.

- stop_on_error:

  raise an `fwe_execution_error` if a local command fails, does not
  produce all required artifacts, or a scheduler submission fails. The
  condition contains the run report in its `run` field.

- persist_results:

  if `TRUE`, real executions and submissions update the persistent
  artifact manifest and RDS result snapshot. Dry runs never write result
  files.

- tracking_sqlite_db, tracking_args:

  optional SQLite database and metadata for registering submitted
  correction jobs in pipeline job tracking.

## Value

an `fwe_run` containing the execution report and refreshed `fwe_plan`.
