# Diagnosing a run of fmri.pipeline

## Overview of `diagnose_pipeline` Function

The
[`diagnose_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/diagnose_pipeline.md)
function provides the user quick access to job tracking data, output
files, and log files to help diagnose faulty pipeline runs. To use the
function, supply it with either a ‘gpa’ object or a path string to a
top-level project directory (containing a ‘project_config.json’ file, a
‘scheduler_scripts’ folder, an SQLite database, and other files used for
an analysis). Note that when a ‘gpa’ object is supplied, the function
will try to access the output locations specified in the object; if a
top-level directory path is supplied, the function will search in
conventionally named folders for applicable information.

In the example below, it is assumed a ‘gpa’ object called `gpa` has
already been constructed in the R environment. To start a diagnosis,
simply execute:

    diagnose_pipeline(gpa)

You should receive a message that says

    Running pipeline diagnosis...
    (Hit ESC at any time to cancel)

indicating the diagnosis has launched successfully. Note that the ESC
key can be used to exit the diagnosis at any time. In the unusual case
that a specific element necessary to diagnose the pipeline (e.g., the
SQLite database) is missing, the function will throw an informative
error and cancel automatically.

### Submission History

The diagnosis function will first print a history of the submissions of
the current analysis including date/time, user, and steps submitted. The
example below is a multi-run FSL submission with a separate L2 step. An
AFNI or SPM submission using its own first-level contrasts goes directly
from L1 to L3.

*Example*:

    ── Analysis: 23june2025 ───────────────────────────────────────────
    • Output directory: /proj/23june2025
    • Number of runs: 1
    • History:
      1. Run on 07/31/25 at 08:58:49 AM by user (Latest)
        Finalize Config. → L1 Mods. → L2 Mods. → L3 Mods. → Clean Up

The function will then ask which run you would like to diagnose. Enter
the integer of the run you would like to examine more closely, or hit
enter to accept the default (the latest submission).

*Example*:

    Which pipeline run would you like to diagnose?
    Enter an integer (Press enter to accept default: 1) > 1

### Status Print Out

Next, the diagnosis function will print the status of each job (and a
summary of child job statuses) in the sequence:

- “\[✔\] COMPLETED” means the job has successfully completed.
- “\[↗\] QUEUED” means the job was submitted but has not started yet.
- “\[⋯\] STARTED” means the job has started but has not completed yet.
- “\[✗\] FAILED” means the job has failed, either because of an error
  within the step or because a parent/upstream job failed. In the latter
  case, this information will be provided.

*Note*: There are two types of “child” jobs in the pipeline: steps that
are downstream of another but will not start until the previous step is
completed (e.g., run L1 models will submit after setup L1 models has
completed), or smaller, “sub-” jobs that are part of the step itself
(e.g., separate feat jobs). These are specified by a child “level” – 1
for the former type and 2 for the latter type – in the tracking
database. All jobs within a series of steps that are dependent on
another (no matter the child level) are called a “chain,” collectively.
A collection of chains/jobs in a single analysis is called a “sequence”
and is given a unique sequence ID in the tracking database.

At this point, you can either examine one of the steps more closely to
determine what went wrong in the case of failure, or return the full
tracking tree if you initially assigned the diagnosis function to an
object (which was not done in the example; more on this second choice
below).

*Example*:

    The run you selected had 1 independent chain of steps:
    1. Finalize Configuration [✔] COMPLETED
    2. Set Up L1 Models [✔] COMPLETED
    3. Run L1 Models [✔] COMPLETED
      Child Jobs:
      └─ 100 of 100 [✔] COMPLETED
    4. Set Up/Run L2 Models [✔] COMPLETED
    5. Set Up/Run L3 Models [✔] COMPLETED
    6. Clean Up Pipeline [✗] FAILED

### Examining a Step

If you’d like to examine a step more closely, type ‘yes’ when prompted.
You will receive a message asking which chain (if applicable) and then
which step you’d like to examine. Enter the integer of the corresponding
chain/step you’d like to examine. The example below had a single chain
so the user only has to make the choice of which step.

*Example*:

    Would you like to examine any of these steps more closely?
    ℹ Type 'yes' to continue, or 'no' to return full job tree and exit now.
    (yes/no) > yes

    Which step would you like to examine more closely?
    ℹ There are 6 steps in chain 1
    Enter an integer > 6

The function will print out detailed information about the selected step
(queue/start/end times, if applicable, downstream jobs, child jobs).

*Example*:

    ── Job `cleanup_fsl`...
    • ...was submitted at 2025-07-31 14:50:14
    • ...and was started at 2025-07-31 16:18:58
    • ...failed before finishing
    • ...had no downstream (dependent) jobs
    • ...had no child jobs

### Performing Actions

The diagnosis function will then list a series of actions that can
possibly be performed to further diagnose any problems with the step.
These actions could include some or all of the following:

- `View output file in console`: This will print the R output file from
  the step in the console. The user can specify how many lines from the
  end to view (default: 15).
- `Return output file as character object`: This will return the R
  output file from the step as a character object if the function was
  assigned to an object initially.
- `View child job details (statuses, etc.)`: If the step had level 2
  child jobs (smaller, sub-jobs that are part of the step), this will
  print out a more detailed status table for each of those child jobs.
- `View feat status table`: This will display which feat jobs have been
  completed and which have failed, with their corresponding ID.
- `View feat log`: This will print a user-specified number (default: 15)
  of lines from the end of the corresponding log file for that step. It
  will also ask if the user wishes to view any individual-level logs by
  specifying the ID.
- `Exit diagnosis and return full tracking tree`: This will end the
  diagnosis function and return the full tracking tree (see below)

### Check the fit behind a scheduler job

The job tree records what the scheduler ran. To check whether a fit
produced the contrast needed downstream, load the batch’s saved `gpa`
after setup and refresh its status from disk. For example, for an SPM
first-level job:

``` r

gpa <- refresh_glm_status(gpa, level = 1L, glm_software = "spm")
gpa$l1_model_setup$spm[, c("id", "session", "spm_dir", "spm_complete", "spm_failed")]
```

Use the same call with `glm_software = "fsl"` or `"afni"` for those
engines. A failed or incomplete row points to the fit to investigate:

| Engine | Status columns in its L1 setup table | What to check |
|----|----|----|
| FSL | `feat_complete`, `feat_failed` | Open the `.fsf` at `feat_fsf` and inspect the FEAT report and `stats/` under `feat_dir`. |
| AFNI | `afni_complete`, `afni_failed` | Read `afni_glm.log` in `afni_dir` for the failing command. Check `contrasts.tsv` for the expected coefficient image. For a group mixed model, use `3dLMEr.log` in the mixed-model directory. |
| SPM | `spm_complete`, `spm_failed` | Check the job output for MATLAB errors, then inspect `SPM.mat` and `spm_contrasts.tsv` in `spm_dir` for the requested contrast. |

For a group job, refresh `level = 3L` and inspect the engine’s table in
`gpa$l3_model_setup`. If group setup finds no inputs, check the producer
in the submission preflight report and refresh the level that supplied
its contrast images. The [engine
guide](https://hallquistlab.github.io/fmri.pipeline/articles/glm-engines.md)
describes the supported combinations.

## Recover a stalled or failed analysis

Start with the first failing stage in the selected submission. Later
steps may be blocked by the same upstream failure. Check the step’s
output and the relevant child job before changing the model or
resources.

| Symptom | Where to look |
|----|----|
| L1 setup did not produce a design | Inspect the run’s timing, volume count, nuisance columns, and exclusion fields in the setup logs. See [run-level design review](https://hallquistlab.github.io/fmri.pipeline/articles/whole-brain.html#build-and-inspect-the-first-level-designs). |
| A fit started and stopped | Compare scheduler exit status and memory/time usage with the engine’s fit log. For FEAT, `job_manifest.tsv` in the fit directory links the artifact to its job and log. |
| L3 setup found no eligible inputs | Check the resolved producer, completed lower-level maps, and requested contrasts. Audit the recorded group inputs with the [FEAT output inventory](https://hallquistlab.github.io/fmri.pipeline/articles/whole-brain.html#find-the-group-map-and-its-inputs). |
| L3 completed and `run_fwe` failed | Inspect the inference worker log and artifact manifest. Retry only correction using [`run_fwe_correction()`](https://hallquistlab.github.io/fmri.pipeline/articles/fwe.html#run-ptfce-automatically-after-l3). |

After correcting a transient fit or resource problem, use the same
explicit model names in
[`run_glm_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_glm_pipeline.md)
and inspect the new submission. For FEAT, check the incomplete/failed
folder actions in `gpa$glm_settings$fsl` if you need to preserve failed
fit directories and their logs. When changing an input, design, or
hypothesis, use a separate analysis name to retain a clear record of
which specification produced each result.

## Function Output/Tracking Tree

It can be beneficial to assign the initial function call to an object as
in the example below.

    output <- diagnose_pipeline(gpa)

In this case, the diagnosis function will be able to return helpful
information corresponding to what the user specifies in the command line
interface. One such piece information is the full job tracking tree, a
readable version of the job tracking information stored in the SQLite
database. The tree comes in a data structure from the `data.tree`
package, an R package for working with nested data.

Printing the tree yields a readable structure with the top level being a
sequence, then chains, then steps/jobs, and finally sub-jobs/level 2
child jobs.

*Example:*

                                      levelName
    1  seq_ae2c17c7-95d8-4dc8-a8c8-e4532334ef47
    2   °--chain1
    3       ¦--finalize_configuration
    4       ¦--setup_l1
    5       ¦--run_l1
    6       ¦   ¦--featsep_l1_1
    7       ¦   ¦--featsep_l1_2
    8       ¦   ¦--featsep_l1_3
    9       ¦   ¦--featsep_l1_4
    10      ¦   °--featsep_l1_5
    11      ¦--setup_run_l2
    12      ¦--setup_run_l3
    13      °--cleanup_fsl

Accessible fields in the tree for each job include: `batch_directory`,
`batch_file`, `child_level`, `code_file`, `compute_file`, `head`, `id`,
`job_id`, `job_name`, `job_obj`, `mem_total`, `n_cpus`, `n_nodes`,
`scheduler`, `sequence_id`, `status`, `time_ended`, `time_started`,
`time_submitted`, `wall_time`, `parent_id`, `scheduler_options`, and
`mem_per_cpu`. For example, simply execute `print(output, "status")` to
view the status of each job in the tree.

Read the
[`data.tree`](https://cran.r-project.org/web/packages/data.tree/vignettes/data.tree.html)
documentation for more information about working with this easy-to-use
data structure.
