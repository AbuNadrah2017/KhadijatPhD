## Run automated checks for the KhadijatPhD Shiny app.
##
## Usage:
##   Rscript scripts/run_tests.R
##
## This repo is not an R package, so we run tests via testthat::test_dir().

suppressPackageStartupMessages({
  library(testthat)
})

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) < 1) return(NULL)
  sub("^--file=", "", file_arg[[1]])
}

find_root_simple <- function(start_dir = getwd(), max_up = 12L) {
  start_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  cur <- start_dir
  for (i in seq_len(max_up)) {
    has_agents <- file.exists(file.path(cur, "AGENTS.md"))
    has_legacy_data <- dir.exists(file.path(cur, "Data"))
    has_r_app <- file.exists(file.path(cur, "R", "app.R"))
    if (has_agents && (has_legacy_data || has_r_app)) return(cur)

    parent <- normalizePath(file.path(cur, ".."), winslash = "/", mustWork = FALSE)
    if (identical(parent, cur)) break
    cur <- parent
  }
  start_dir
}

script_path <- get_script_path()
proj_guess <- if (!is.null(script_path)) {
  script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
  file.path(script_dir, "..")
} else {
  getwd()
}

proj <- find_root_simple(proj_guess)
setwd(proj)
source(file.path(proj, "R", "utils.R"))

proj <- find_project_root(getwd())
setwd(proj)
options(khadijat_phd_root = proj)

cat("Project root:", proj, "\n")
cat("Running tests...\n")

res <- testthat::test_dir("tests/testthat", reporter = "summary")

# Summarize expectations (success/failure/error/warning) in a stable way.
all_results <- unlist(lapply(res, function(x) x$results), recursive = FALSE, use.names = FALSE)
is_fail <- vapply(all_results, inherits, logical(1), "expectation_failure")
is_err <- vapply(all_results, inherits, logical(1), "expectation_error")
is_warn <- vapply(all_results, inherits, logical(1), "expectation_warning")

cat("\n")
cat("Done.\n")
cat("Failures:", sum(is_fail), "Errors:", sum(is_err), "Warnings:", sum(is_warn), "\n")
