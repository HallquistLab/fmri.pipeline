# Exclude the three development-pilot seed indices and append three fresh ones.
# Original pass: 40 datasets, seed base 9420000.
# Continuation: 3 datasets, seed base 9420040 (indices 41--43).
# Rscript --vanilla inst/validation/assemble_hrf_noise_validation.R first extra out
source(file.path("inst", "validation", "hrf_recovery.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
first <- readRDS(file.path(args[1], "settings.rds"))
extra <- readRDS(file.path(args[2], "settings.rds"))
stopifnot(first$repetitions == 40L, extra$repetitions == 3L,
  first$seed_base == 9420000L, extra$seed_base == 9420040L,
  identical(first$scenarios, extra$scenarios), identical(first$lambdas, extra$lambdas),
  identical(first$script_md5, extra$script_md5))
out <- args[3]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
join <- function(file) {
  a <- read.csv(file.path(args[1], file)); b <- read.csv(file.path(args[2], file))
  b$replicate <- b$replicate + 40L
  rbind(a[a$replicate > 3L, ], b)
}
metrics <- join("replicate_metrics.csv"); runs <- join("run_metrics.csv")
stopifnot(all(table(metrics$scenario, metrics$method, metrics$lambda) == 40L))
write.csv(metrics, file.path(out, "replicate_metrics.csv"), row.names = FALSE)
write.csv(runs, file.path(out, "run_metrics.csv"), row.names = FALSE)
tables <- hrfv_summary(metrics)
matched <- metrics[metrics$method %in% c("canonical_gls", "gls_heuristic", "gls_calibrated"), ]
matched$method[matched$method == "canonical_gls"] <- "canonical"
write.csv(tables$summary, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(tables$paired, file.path(out, "paired_vs_legacy_canonical.csv"), row.names = FALSE)
write.csv(hrfv_summary(matched)$paired, file.path(out, "paired_vs_matched_canonical.csv"), row.names = FALSE)
first$replicate_indices <- 4:43
first$excluded_development_pilot_indices <- 1:3
first$continuation_seed_base <- extra$seed_base
saveRDS(first, file.path(out, "settings.rds"))
writeLines(c("Final analysis uses 40 fresh datasets per scenario: seed indices 4--43.",
  "Development-pilot indices 1--3 are excluded from the reported analysis.",
  "The continuation's replicate 1--3 correspond exactly to indices 41--43.",
  "Noise/evidence/regularization settings were held fixed across both passes.",
  readLines(file.path(args[1], "environment.txt")),
  "Continuation environment:", readLines(file.path(args[2], "environment.txt"))),
  file.path(out, "environment.txt"))
