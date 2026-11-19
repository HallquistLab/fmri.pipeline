if (identical(.Platform$OS.type, "windows")) {
  test_that("FLAME runner behavioral tests require POSIX shell semantics", {
    skip("requires executable shell scripts and POSIX process semantics")
  })
} else {

test_that("flame_runner uses parallel and clamps to allocation", {
  tmp <- withr::local_tempdir()
  stub_dir <- file.path(tmp, "stubs_parallel")
  dir.create(stub_dir, recursive = TRUE, showWarnings = FALSE)

  make_stub <- function(name, lines) {
    path <- file.path(stub_dir, name)
    writeLines(lines, path)
    Sys.chmod(path, "0755")
    path
  }

  make_stub("sh", c(
    "#!/bin/bash",
    "exec /bin/sh \"$@\""
  ))
  make_stub("env", c(
    "#!/bin/bash",
    "exec /usr/bin/env \"$@\""
  ))
  make_stub("grep", c(
    "#!/bin/bash",
    "exec /usr/bin/grep \"$@\""
  ))

  parallel_log <- file.path(tmp, "parallel.log")
  make_stub("parallel", c(
    "#!/bin/bash",
    "jobs=1",
    "while [[ $# -gt 0 ]]; do",
    "  case \"$1\" in",
    "    --jobs|-j) jobs=\"$2\"; shift 2;;",
    "    *) break;;",
    "  esac",
    "done",
    "echo \"$jobs\" > \"$PARALLEL_LOG\"",
    "cmd=(\"$@\")",
    "while IFS= read -r line; do",
    "  if [ -n \"${line//[[:space:]]/}\" ]; then",
    "    run=(\"${cmd[@]}\")",
    "    replaced=0",
    "    for i in \"${!run[@]}\"; do",
    "      if [[ \"${run[$i]}\" == \"{}\" ]]; then",
    "        run[$i]=\"$line\"",
    "        replaced=1",
    "      fi",
    "    done",
    "    if [[ \"$replaced\" -eq 1 ]]; then",
    "      \"${run[@]}\"",
    "    else",
    "      \"${run[@]}\" \"$line\"",
    "    fi",
    "  fi",
    "done"
  ))

  flame_list <- file.path(tmp, "flame.list")
  out_file <- file.path(tmp, "out.txt")
  writeLines(c(
    paste("echo A >>", out_file),
    paste("echo B >>", out_file),
    paste("echo C >>", out_file)
  ), flame_list)

  withr::local_envvar(c(
    PATH = stub_dir,
    PARALLEL_LOG = parallel_log,
    SLURM_CPUS_PER_TASK = NA_character_,
    SLURM_NTASKS = "2",
    SLURM_CPUS_ON_NODE = NA_character_,
    PBS_NP = NA_character_,
    NCPUS = NA_character_,
    PBS_NUM_PPN = NA_character_
  ))

  flame_runner <- source_tree_file("inst", "bin", "flame_runner")
  system2(flame_runner, c(flame_list, "999"), stdout = TRUE, stderr = TRUE)

  expect_equal(readLines(parallel_log), "2")
  expect_true(file.exists(out_file))
  expect_equal(length(readLines(out_file)), 3)
})

test_that("flame_runner uses xargs when parallel is unavailable", {
  tmp <- withr::local_tempdir()
  stub_dir <- file.path(tmp, "stubs_xargs")
  dir.create(stub_dir, recursive = TRUE, showWarnings = FALSE)

  make_stub <- function(name, lines) {
    path <- file.path(stub_dir, name)
    writeLines(lines, path)
    Sys.chmod(path, "0755")
    path
  }

  make_stub("sh", c(
    "#!/bin/bash",
    "exec /bin/sh \"$@\""
  ))
  make_stub("env", c(
    "#!/bin/bash",
    "exec /usr/bin/env \"$@\""
  ))
  make_stub("grep", c(
    "#!/bin/bash",
    "exec /usr/bin/grep \"$@\""
  ))

  xargs_log <- file.path(tmp, "xargs.log")
  make_stub("xargs", c(
    "#!/bin/bash",
    "if [[ \"$1\" == \"--help\" ]]; then",
    "  echo \"  -d delim\"",
    "  exit 0",
    "fi",
    "p=1",
    "while [[ $# -gt 0 ]]; do",
    "  case \"$1\" in",
    "    -P) p=\"$2\"; shift 2;;",
    "    -d) shift 2;;",
    "    -n) shift 2;;",
    "    *) break;;",
    "  esac",
    "done",
    "echo \"$p\" > \"$XARGS_LOG\"",
    "cmd=(\"$@\")",
    "while IFS= read -r line; do",
    "  if [ -n \"${line//[[:space:]]/}\" ]; then",
    "    \"${cmd[@]}\" \"$line\"",
    "  fi",
    "done"
  ))

  flame_list <- file.path(tmp, "flame.list")
  out_file <- file.path(tmp, "out.txt")
  writeLines(c(
    paste("echo A >>", out_file),
    paste("echo B >>", out_file)
  ), flame_list)

  withr::local_envvar(c(
    PATH = stub_dir,
    XARGS_LOG = xargs_log
  ))

  flame_runner <- source_tree_file("inst", "bin", "flame_runner")
  system2(flame_runner, c(flame_list, "4"), stdout = TRUE, stderr = TRUE)

  expect_equal(readLines(xargs_log), "2")
  expect_true(file.exists(out_file))
  expect_equal(length(readLines(out_file)), 2)
})

test_that("flame_runner retries failed FLAME12 slices with FLAME1", {
  tmp <- withr::local_tempdir()
  stub_dir <- file.path(tmp, "stubs_flame1_retry")
  dir.create(stub_dir, recursive = TRUE, showWarnings = FALSE)

  make_stub <- function(name, lines) {
    path <- file.path(stub_dir, name)
    writeLines(lines, path)
    Sys.chmod(path, "0755")
    path
  }

  make_stub("grep", c(
    "#!/bin/bash",
    "exec /usr/bin/grep \"$@\""
  ))
  make_stub("sh", c(
    "#!/bin/bash",
    "exec /bin/sh \"$@\""
  ))
  make_stub("xargs", c(
    "#!/bin/bash",
    "if [[ \"$1\" == \"--help\" ]]; then",
    "  echo \"xargs help\"",
    "  exit 0",
    "fi",
    "exit 1"
  ))

  fake_flameo <- file.path(tmp, "fake_flameo")
  attempt_log <- file.path(tmp, "attempts.log")
  fallback_log <- file.path(tmp, "fallbacks.tsv")
  writeLines(c(
    "#!/bin/bash",
    "printf '%s\\n' \"$*\" >> \"$FAKE_FLAMEO_LOG\"",
    "case \" $* \" in",
    "  *' --runmode=flame12 '*) exit 42;;",
    "  *' --runmode=flame1 '*)",
    "    case \" $* \" in",
    "      *' --nj='*|*' --bi='*|*' --se='*|*' --zlt='*|*' --zut='*|*' --fm'*) exit 43;;",
    "    esac",
    '    exit "${FAKE_FLAME1_STATUS:-0}";;',
    "esac",
    "exit 44"
  ), fake_flameo)
  Sys.chmod(fake_flameo, "0755")

  flame_list <- file.path(tmp, "flame.list")
  writeLines(
    paste(
      fake_flameo,
      "--cope=tmpcope0007 --vc=tmpvarcope0007 --mask=tmpmask0007",
      "--ld=stats0007 --dm=design.mat --cs=design.grp --tc=design.con",
      "--runmode=flame12 --nj=10000 --bi=500 --se=1 --fm --zlt=3.05 --zut=3.45"
    ),
    flame_list
  )

  withr::local_envvar(c(
    PATH = stub_dir,
    FAKE_FLAMEO_LOG = attempt_log,
    FLAME_RUNNER_FALLBACK_LOG = fallback_log
  ))

  flame_runner <- source_tree_file("inst", "bin", "flame_runner")
  status <- system2(flame_runner, c(flame_list, "1"), stdout = TRUE, stderr = TRUE)

  expect_false(inherits(status, "status"))
  attempts <- readLines(attempt_log)
  expect_length(attempts, 2)
  expect_true(grepl("--runmode=flame12", attempts[[1]], fixed = TRUE))
  expect_true(grepl("--runmode=flame1", attempts[[2]], fixed = TRUE))
  expect_false(grepl("--nj=|--bi=|--se=|--fm|--zlt=|--zut=", attempts[[2]]))

  fallback_records <- readLines(fallback_log)
  audit_file <- file.path(tmp, "flame_runner_fallbacks.tsv")
  audit <- readLines(audit_file)
  expect_identical(sub("\t[^\t]+$", "", audit), fallback_records)
  expect_match(audit[[2]], "\t0$")
  expect_length(fallback_records, 2)
  expect_equal(
    fallback_records[[1]],
    "timestamp_utc\tpid\tflame12_status\tld_dir\tfailed_ld_dir\tfallback_command"
  )
  expect_true(grepl("\t42\tstats0007\t", fallback_records[[2]], fixed = TRUE))
  expect_true(grepl("--runmode=flame1", fallback_records[[2]], fixed = TRUE))

  withr::local_envvar(c(FAKE_FLAME1_STATUS = "43"))
  failed <- suppressWarnings(system2(flame_runner, c(flame_list, "1"), stdout = TRUE, stderr = TRUE))
  expect_true(attr(failed, "status") > 0L)
  expect_identical(utils::read.delim(audit_file)$fallback_status, 43L)
})

test_that("flame_runner reports unrecovered non-FLAME12 failures", {
  tmp <- withr::local_tempdir()
  stub_dir <- file.path(tmp, "stubs_failure")
  dir.create(stub_dir, recursive = TRUE, showWarnings = FALSE)

  make_stub <- function(name, lines) {
    path <- file.path(stub_dir, name)
    writeLines(lines, path)
    Sys.chmod(path, "0755")
    path
  }

  make_stub("grep", c(
    "#!/bin/bash",
    "exec /usr/bin/grep \"$@\""
  ))
  make_stub("sh", c(
    "#!/bin/bash",
    "exec /bin/sh \"$@\""
  ))
  make_stub("xargs", c(
    "#!/bin/bash",
    "if [[ \"$1\" == \"--help\" ]]; then",
    "  echo \"xargs help\"",
    "  exit 0",
    "fi",
    "exit 1"
  ))

  flame_list <- file.path(tmp, "flame.list")
  writeLines("exit 9", flame_list)

  withr::local_envvar(c(PATH = stub_dir))

  flame_runner <- source_tree_file("inst", "bin", "flame_runner")
  expect_warning(
    status <- system2(flame_runner, c(flame_list, "1"), stdout = TRUE, stderr = TRUE),
    "had status"
  )

  expect_equal(attr(status, "status"), 1)
})

test_that("flame_runner falls back when xargs lacks -d support", {
  tmp <- withr::local_tempdir()
  stub_dir <- file.path(tmp, "stubs_fallback")
  dir.create(stub_dir, recursive = TRUE, showWarnings = FALSE)

  make_stub <- function(name, lines) {
    path <- file.path(stub_dir, name)
    writeLines(lines, path)
    Sys.chmod(path, "0755")
    path
  }

  make_stub("sh", c(
    "#!/bin/bash",
    "exec /bin/sh \"$@\""
  ))
  make_stub("grep", c(
    "#!/bin/bash",
    "exec /usr/bin/grep \"$@\""
  ))
  make_stub("xargs", c(
    "#!/bin/bash",
    "if [[ \"$1\" == \"--help\" ]]; then",
    "  echo \"xargs help\"",
    "  exit 0",
    "fi",
    "exit 1"
  ))

  flame_list <- file.path(tmp, "flame.list")
  out_file <- file.path(tmp, "out.txt")
  writeLines(c(
    paste("echo A >>", out_file),
    paste("echo B >>", out_file)
  ), flame_list)

  withr::local_envvar(c(
    PATH = stub_dir
  ))

  flame_runner <- source_tree_file("inst", "bin", "flame_runner")
  system2(flame_runner, c(flame_list, "10"), stdout = TRUE, stderr = TRUE)

  expect_true(file.exists(out_file))
  expect_equal(length(readLines(out_file)), 2)
})
}
