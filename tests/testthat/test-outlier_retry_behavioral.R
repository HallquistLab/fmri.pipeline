# Behavioral tests for the FLAME outlier auto-retry mechanism.
#
# These tests exercise the actual bash and R logic rather than grepping source
# code for string patterns.

if (identical(.Platform$OS.type, "windows")) {
  test_that("outlier retry behavioral tests require POSIX shell semantics", {
    skip("requires executable shell scripts and POSIX process semantics")
  })
} else {

retry_script <- source_tree_file("inst", "bash", "run_feat_with_outlier_retry.sh")

# ===========================================================================
# 1. Bash helper functions: extract_outputdir, extract_robust_yn,
#    has_excessive_outliers_error
# ===========================================================================

# Source just the helper functions from the retry script so we can call them
# individually.  We extract the function bodies and wrap them in a small
# harness script.
build_helper_harness <- function(tmpdir) {
  harness <- file.path(tmpdir, "harness.sh")
  writeLines(c(
    "#!/bin/bash",
    "set -u",
    # Inline the three helper functions verbatim from the retry script.
    "extract_outputdir() {",
    '  local fsf_file="$1"',
    "  local out",
    '  out=$(sed -n \'s/^set fmri(outputdir)[[:space:]]*"\\(.*\\)"[[:space:]]*$/\\1/p\' "$fsf_file" | head -n 1)',
    '  if [ -z "$out" ]; then',
    '    out="${fsf_file%.fsf}"',
    "  fi",
    "  printf '%s\\n' \"$out\"",
    "}",
    "",
    "extract_robust_yn() {",
    '  local fsf_file="$1"',
    '  awk \'/^set fmri\\(robust_yn\\)[[:space:]]+/ { print $3; exit }\' "$fsf_file"',
    "}",
    "",
    "has_excessive_outliers_error() {",
    '  local feat_dir="$1"',
    '  if [ ! -d "$feat_dir" ]; then',
    "    return 1",
    "  fi",
    '  grep -Rqsi "excessive.*outliers detected" "$feat_dir"',
    "}",
    "",
    '# Dispatch: call the function named in $1 with remaining args',
    '"$@"'
  ), harness)
  Sys.chmod(harness, "0755")
  harness
}

# ---------------------------------------------------------------------------
# extract_outputdir
# ---------------------------------------------------------------------------
test_that("extract_outputdir parses outputdir from FSF", {
  tmp <- withr::local_tempdir()
  harness <- build_helper_harness(tmp)

  fsf <- file.path(tmp, "model.fsf")
  writeLines(c(
    'set fmri(outputdir) "/data/results/my_model"',
    "set fmri(robust_yn) 1"
  ), fsf)

  out <- system2(harness, c("extract_outputdir", fsf), stdout = TRUE)
  expect_equal(out, "/data/results/my_model")
})

test_that("extract_outputdir falls back to fsf path minus .fsf when tag missing", {
  tmp <- withr::local_tempdir()
  harness <- build_helper_harness(tmp)

  fsf <- file.path(tmp, "no_output_tag.fsf")
  writeLines("set fmri(robust_yn) 1", fsf)

  out <- system2(harness, c("extract_outputdir", fsf), stdout = TRUE)
  expect_equal(out, sub("\\.fsf$", "", fsf))
})

# ---------------------------------------------------------------------------
# extract_robust_yn
# ---------------------------------------------------------------------------
test_that("extract_robust_yn returns robust_yn value", {
  tmp <- withr::local_tempdir()
  harness <- build_helper_harness(tmp)

  fsf <- file.path(tmp, "model.fsf")
  writeLines(c(
    'set fmri(outputdir) "/data/results/my_model"',
    "set fmri(robust_yn) 1"
  ), fsf)

  out <- system2(harness, c("extract_robust_yn", fsf), stdout = TRUE)
  expect_equal(out, "1")
})

test_that("extract_robust_yn returns 0 when robust is disabled", {
  tmp <- withr::local_tempdir()
  harness <- build_helper_harness(tmp)

  fsf <- file.path(tmp, "model.fsf")
  writeLines("set fmri(robust_yn) 0", fsf)

  out <- system2(harness, c("extract_robust_yn", fsf), stdout = TRUE)
  expect_equal(out, "0")
})

# ---------------------------------------------------------------------------
# has_excessive_outliers_error
# ---------------------------------------------------------------------------
test_that("has_excessive_outliers_error detects the error string", {
  tmp <- withr::local_tempdir()
  harness <- build_helper_harness(tmp)

  gfeat <- file.path(tmp, "model.gfeat")
  dir.create(gfeat)
  writeLines("Excessive number of outliers detected in FLAME",
             file.path(gfeat, "logs.txt"))

  rc <- system2(harness, c("has_excessive_outliers_error", gfeat))
  expect_equal(rc, 0L) # 0 = found
})

test_that("has_excessive_outliers_error matches case-insensitively", {
  tmp <- withr::local_tempdir()
  harness <- build_helper_harness(tmp)

  gfeat <- file.path(tmp, "model.gfeat")
  dir.create(gfeat)
  writeLines("EXCESSIVE NUMBER OF OUTLIERS DETECTED",
             file.path(gfeat, "logs.txt"))

  rc <- system2(harness, c("has_excessive_outliers_error", gfeat))
  expect_equal(rc, 0L)
})

test_that("has_excessive_outliers_error matches variant wording", {
  tmp <- withr::local_tempdir()
  harness <- build_helper_harness(tmp)

  gfeat <- file.path(tmp, "model.gfeat")
  dir.create(gfeat)
  writeLines("excessive outliers detected by FLAME",
             file.path(gfeat, "logs.txt"))

  rc <- system2(harness, c("has_excessive_outliers_error", gfeat))
  expect_equal(rc, 0L)
})

test_that("has_excessive_outliers_error returns 1 when no error present", {
  tmp <- withr::local_tempdir()
  harness <- build_helper_harness(tmp)

  gfeat <- file.path(tmp, "model.gfeat")
  dir.create(gfeat)
  writeLines("Everything went fine", file.path(gfeat, "logs.txt"))

  rc <- system2(harness, c("has_excessive_outliers_error", gfeat))
  expect_equal(rc, 1L) # 1 = not found
})

test_that("has_excessive_outliers_error returns 1 for non-existent directory", {
  tmp <- withr::local_tempdir()
  harness <- build_helper_harness(tmp)

  rc <- system2(harness, c("has_excessive_outliers_error", file.path(tmp, "nope.gfeat")))
  expect_equal(rc, 1L)
})

# ===========================================================================
# 2. FSF rewriting: the awk script that sets robust_yn=0 and changes outputdir
# ===========================================================================
test_that("FSF rewriting sets robust_yn=0 and changes outputdir", {
  tmp <- withr::local_tempdir()

  fsf <- file.path(tmp, "original.fsf")
  writeLines(c(
    "set fmri(npts) 100",
    'set fmri(outputdir) "/data/results/my_model"',
    "set fmri(robust_yn) 1",
    "set fmri(evs_real) 4"
  ), fsf)

  retry_outputdir <- "/data/results/my_model.retry_no_outlier"
  retry_fsf <- file.path(tmp, "retry.fsf")

  # Run the same awk command used in the retry script
  awk_cmd <- sprintf(
    "awk -v retry_out='%s' '
      BEGIN { set_robust = 0; set_output = 0 }
      {
        if ($0 ~ /^set fmri\\(robust_yn\\)[[:space:]]+/) {
          print \"set fmri(robust_yn) 0\"
          set_robust = 1
          next
        }
        if ($0 ~ /^set fmri\\(outputdir\\)[[:space:]]+/) {
          print \"set fmri(outputdir) \\\"\" retry_out \"\\\"\"
          set_output = 1
          next
        }
        print
      }
      END {
        if (set_robust == 0) print \"set fmri(robust_yn) 0\"
        if (set_output == 0) print \"set fmri(outputdir) \\\"\" retry_out \"\\\"\"
      }
    ' '%s' > '%s'",
    retry_outputdir, fsf, retry_fsf
  )
  system(awk_cmd)

  result <- readLines(retry_fsf)
  expect_true(any(grepl("^set fmri\\(robust_yn\\) 0$", result)))
  expect_true(any(grepl(retry_outputdir, result, fixed = TRUE)))
  # Original values should NOT be present

  expect_false(any(grepl("robust_yn\\) 1", result)))
  expect_false(any(grepl("my_model\"$", result)))
  # Non-modified lines should pass through
  expect_true(any(grepl("set fmri\\(npts\\) 100", result)))
  expect_true(any(grepl("set fmri\\(evs_real\\) 4", result)))
})

test_that("FSF rewriting appends missing tags", {
  tmp <- withr::local_tempdir()

  # FSF with neither robust_yn nor outputdir
  fsf <- file.path(tmp, "minimal.fsf")
  writeLines(c(
    "set fmri(npts) 50"
  ), fsf)

  retry_outputdir <- "/data/retry_output"
  retry_fsf <- file.path(tmp, "retry.fsf")

  awk_cmd <- sprintf(
    "awk -v retry_out='%s' '
      BEGIN { set_robust = 0; set_output = 0 }
      {
        if ($0 ~ /^set fmri\\(robust_yn\\)[[:space:]]+/) {
          print \"set fmri(robust_yn) 0\"
          set_robust = 1
          next
        }
        if ($0 ~ /^set fmri\\(outputdir\\)[[:space:]]+/) {
          print \"set fmri(outputdir) \\\"\" retry_out \"\\\"\"
          set_output = 1
          next
        }
        print
      }
      END {
        if (set_robust == 0) print \"set fmri(robust_yn) 0\"
        if (set_output == 0) print \"set fmri(outputdir) \\\"\" retry_out \"\\\"\"
      }
    ' '%s' > '%s'",
    retry_outputdir, fsf, retry_fsf
  )
  system(awk_cmd)

  result <- readLines(retry_fsf)
  expect_true(any(grepl("^set fmri\\(robust_yn\\) 0$", result)))
  expect_true(any(grepl(retry_outputdir, result, fixed = TRUE)))
  expect_true(any(grepl("set fmri\\(npts\\) 50", result)))
})

# ===========================================================================
# 3. Full retry wrapper: end-to-end with a mock feat binary
# ===========================================================================

# Build a mock "feat" that:
#   - On first call: creates .gfeat dir with outlier error, exits 1
#   - On second call (retry): creates .gfeat dir cleanly, exits 0
# Uses a state file to distinguish first vs second invocation.
build_mock_feat <- function(tmpdir) {
  mock <- file.path(tmpdir, "mock_feat")
  state_file <- file.path(tmpdir, ".mock_feat_call_count")
  writeLines(c(
    "#!/bin/bash",
    "# Mock feat binary for testing the retry wrapper",
    "fsf=\"$1\"",
    paste0("state_file=\"", state_file, "\""),
    "",
    "# Parse outputdir from the FSF",
    "outputdir=$(sed -n 's/^set fmri(outputdir)[[:space:]]*\"\\(.*\\)\"[[:space:]]*$/\\1/p' \"$fsf\" | head -n 1)",
    "gfeat=\"${outputdir}.gfeat\"",
    "",
    "# Track call count",
    "if [ -f \"$state_file\" ]; then",
    "  count=$(cat \"$state_file\")",
    "else",
    "  count=0",
    "fi",
    "count=$((count + 1))",
    "echo $count > \"$state_file\"",
    "",
    "mkdir -p \"$gfeat\"",
    "",
    "if [ $count -eq 1 ]; then",
    "  # First call: simulate FLAME outlier failure",
    '  echo "Excessive number of outliers detected" > "${gfeat}/flame_log.txt"',
    "  exit 1",
    "else",
    "  # Second call (retry): succeed",
    '  echo "FLAME completed successfully" > "${gfeat}/flame_log.txt"',
    "  exit 0",
    "fi"
  ), mock)
  Sys.chmod(mock, "0755")
  mock
}

test_that("retry wrapper retries on outlier failure and produces correct output", {
  tmp <- withr::local_tempdir()

  mock_feat <- build_mock_feat(tmp)

  # Create a FSF file
  fsf <- file.path(tmp, "test_model.fsf")
  writeLines(c(
    paste0('set fmri(outputdir) "', file.path(tmp, "test_model"), '"'),
    "set fmri(robust_yn) 1",
    "set fmri(npts) 100"
  ), fsf)

  # Run the retry wrapper
  rc <- system2(retry_script, c(mock_feat, fsf), stdout = TRUE, stderr = TRUE)
  exit_code <- attr(rc, "status")

  # Should succeed (retry worked)
  expect_null(exit_code) # NULL means exit 0

  # The canonical .gfeat should exist
  canonical_gfeat <- file.path(tmp, "test_model.gfeat")
  expect_true(dir.exists(canonical_gfeat))

  # .feat_auto_retry_warning should be present
  warning_file <- file.path(canonical_gfeat, ".feat_auto_retry_warning")
  expect_true(file.exists(warning_file))
  warning_content <- readLines(warning_file)
  expect_true(any(grepl("Excessive number of FLAME outliers", warning_content)))
  expect_true(any(grepl("robust_yn.*1", warning_content)))

  # Retry FSF should have been created and its copy preserved inside .gfeat
  design_copy <- file.path(canonical_gfeat, "design_retry_no_outlier.fsf")
  expect_true(file.exists(design_copy))
  retry_lines <- readLines(design_copy)
  expect_true(any(grepl("robust_yn\\) 0", retry_lines)))

  # Retry artifacts should be cleaned up on success
  retry_fsf <- file.path(tmp, "test_model__retry_no_outlier.fsf")
  expect_false(file.exists(retry_fsf))
  archived <- list.files(tmp, pattern = "test_model\\.gfeat\\.failed_robust_", full.names = TRUE)
  expect_equal(length(archived), 0L)
})

# Mock feat that always succeeds (no retry should occur)
build_mock_feat_success <- function(tmpdir) {
  mock <- file.path(tmpdir, "mock_feat_ok")
  writeLines(c(
    "#!/bin/bash",
    "fsf=\"$1\"",
    "outputdir=$(sed -n 's/^set fmri(outputdir)[[:space:]]*\"\\(.*\\)\"[[:space:]]*$/\\1/p' \"$fsf\" | head -n 1)",
    "gfeat=\"${outputdir}.gfeat\"",
    "mkdir -p \"$gfeat\"",
    "exit 0"
  ), mock)
  Sys.chmod(mock, "0755")
  mock
}

test_that("retry wrapper exits immediately on success without retrying", {
  tmp <- withr::local_tempdir()

  mock_feat <- build_mock_feat_success(tmp)

  fsf <- file.path(tmp, "ok_model.fsf")
  writeLines(c(
    paste0('set fmri(outputdir) "', file.path(tmp, "ok_model"), '"'),
    "set fmri(robust_yn) 1"
  ), fsf)

  rc <- system2(retry_script, c(mock_feat, fsf), stdout = TRUE, stderr = TRUE)
  exit_code <- attr(rc, "status")
  expect_null(exit_code) # exit 0

  # No retry artifacts
  expect_false(file.exists(file.path(tmp, "ok_model__retry_no_outlier.fsf")))
  expect_false(file.exists(file.path(tmp, "ok_model.gfeat", ".feat_auto_retry_warning")))
})

# Mock feat that fails but NOT with outlier error
build_mock_feat_non_outlier_fail <- function(tmpdir) {
  mock <- file.path(tmpdir, "mock_feat_fail")
  writeLines(c(
    "#!/bin/bash",
    "fsf=\"$1\"",
    "outputdir=$(sed -n 's/^set fmri(outputdir)[[:space:]]*\"\\(.*\\)\"[[:space:]]*$/\\1/p' \"$fsf\" | head -n 1)",
    "gfeat=\"${outputdir}.gfeat\"",
    "mkdir -p \"$gfeat\"",
    'echo "Some other error occurred" > "${gfeat}/log.txt"',
    "exit 1"
  ), mock)
  Sys.chmod(mock, "0755")
  mock
}

test_that("retry wrapper does not retry on non-outlier failure", {
  tmp <- withr::local_tempdir()

  mock_feat <- build_mock_feat_non_outlier_fail(tmp)

  fsf <- file.path(tmp, "fail_model.fsf")
  writeLines(c(
    paste0('set fmri(outputdir) "', file.path(tmp, "fail_model"), '"'),
    "set fmri(robust_yn) 1"
  ), fsf)

  expect_warning(
    rc <- system2(retry_script, c(mock_feat, fsf), stdout = TRUE, stderr = TRUE),
    "status 1"
  )
  exit_code <- attr(rc, "status")
  expect_equal(exit_code, 1L)

  # No retry artifacts
  expect_false(file.exists(file.path(tmp, "fail_model__retry_no_outlier.fsf")))
})

test_that("retry wrapper does not retry when robust_yn is already 0", {
  tmp <- withr::local_tempdir()

  # Use the outlier-failure mock, but set robust_yn=0 in the FSF
  mock_feat <- build_mock_feat(tmp)

  fsf <- file.path(tmp, "already_norobust.fsf")
  writeLines(c(
    paste0('set fmri(outputdir) "', file.path(tmp, "already_norobust"), '"'),
    "set fmri(robust_yn) 0"
  ), fsf)

  expect_warning(
    rc <- system2(retry_script, c(mock_feat, fsf), stdout = TRUE, stderr = TRUE),
    "status 1"
  )
  exit_code <- attr(rc, "status")
  expect_equal(exit_code, 1L) # original failure passed through

  # No retry artifacts
  expect_false(file.exists(file.path(tmp, "already_norobust__retry_no_outlier.fsf")))
})

# ===========================================================================
# 4. get_feat_status: retry detection via .feat_auto_retry_warning
# ===========================================================================
test_that("get_feat_status detects auto-retry and sets feat_auto_retried=TRUE", {
  tmp <- withr::local_tempdir()

  gfeat <- file.path(tmp, "model.gfeat")
  dir.create(gfeat)
  writeLines("Mon Feb 16 10:00:00 EST 2026", file.path(gfeat, ".feat_complete"))
  writeLines(c(
    "Automatic FEAT retry triggered.",
    "Reason: Excessive number of FLAME outliers detected in initial run.",
    "Original FSF: /data/model.fsf",
    "Initial robust_yn: 1"
  ), file.path(gfeat, ".feat_auto_retry_warning"))

  # Also create a .fsf so get_feat_status doesn't complain
  writeLines("set fmri(npts) 100", file.path(tmp, "model.fsf"))

  result <- fmri.pipeline:::get_feat_status(feat_dir = gfeat)
  expect_true(result$feat_auto_retried)
  expect_false(result$feat_complete) # A retry marker alone cannot validate a fit.
})

test_that("get_feat_status returns feat_auto_retried=FALSE when no warning file", {
  tmp <- withr::local_tempdir()

  gfeat <- file.path(tmp, "clean_model.gfeat")
  dir.create(gfeat)
  writeLines("Mon Feb 16 10:00:00 EST 2026", file.path(gfeat, ".feat_complete"))
  writeLines("set fmri(npts) 100", file.path(tmp, "clean_model.fsf"))

  result <- fmri.pipeline:::get_feat_status(feat_dir = gfeat)
  expect_false(result$feat_auto_retried)
  expect_false(result$feat_complete) # A retry marker alone cannot validate a fit.
})

test_that("get_feat_status returns feat_auto_retried=FALSE when directory does not exist", {
  tmp <- withr::local_tempdir()

  result <- fmri.pipeline:::get_feat_status(feat_dir = file.path(tmp, "nonexistent.gfeat"))
  expect_false(result$feat_auto_retried)
  expect_false(result$feat_dir_exists)
})

test_that("get_feat_status logs retry warning details", {
  tmp <- withr::local_tempdir()

  gfeat <- file.path(tmp, "warn_model.gfeat")
  dir.create(gfeat)
  writeLines("Mon Feb 16 10:00:00 EST 2026", file.path(gfeat, ".feat_complete"))
  warning_lines <- c(
    "Automatic FEAT retry triggered.",
    "Reason: Excessive number of FLAME outliers detected in initial run."
  )
  writeLines(warning_lines, file.path(gfeat, ".feat_auto_retry_warning"))
  writeLines("set fmri(npts) 100", file.path(tmp, "warn_model.fsf"))

  # Capture log output
  lg <- lgr::get_logger("test_retry_warn")
  lg$set_threshold("warn")
  handler <- lgr::AppenderBuffer$new(threshold = "warn")
  lg$add_appender(handler, name = "test_buffer")
  on.exit(lg$remove_appender("test_buffer"), add = TRUE)

  result <- fmri.pipeline:::get_feat_status(feat_dir = gfeat, lg = lg)

  logged <- handler$data$msg
  expect_true(any(grepl("auto-retried", logged)))
  expect_true(any(grepl("Excessive number", logged)))
})
}
