## Utilities shared across the Shiny app.
##
## Keep these functions dependency-light so they can be sourced early.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

find_project_root <- function(start_dir = getwd(), max_up = 8L) {
  # Walk up until we find the repo root markers.
  start_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  cur <- start_dir
  for (i in seq_len(max_up)) {
    has_agents <- file.exists(file.path(cur, "AGENTS.md"))
    has_legacy_data <- dir.exists(file.path(cur, "Data"))
    has_r_app <- file.exists(file.path(cur, "R", "app.R"))
    if (has_agents && (has_legacy_data || has_r_app)) {
      return(cur)
    }
    parent <- normalizePath(file.path(cur, ".."), winslash = "/", mustWork = FALSE)
    if (identical(parent, cur)) break
    cur <- parent
  }
  start_dir
}

copy_to_tempfile <- function(path) {
  stopifnot(is.character(path), length(path) == 1L)
  ext <- tools::file_ext(path)
  tmp <- tempfile(fileext = if (nzchar(ext)) paste0(".", ext) else "")

  ok <- file.copy(path, tmp, overwrite = TRUE)
  if (!ok && tolower(Sys.info()[["sysname"]] %||% "") == "windows") {
    # Work around some Windows sharing/locking edge-cases by delegating to PowerShell.
    src <- normalizePath(path, winslash = "\\", mustWork = FALSE)
    dst <- normalizePath(tmp, winslash = "\\", mustWork = FALSE)
    cmd <- paste0("Copy-Item -LiteralPath ", shQuote(src), " -Destination ", shQuote(dst), " -Force")
    status <- suppressWarnings(system2("powershell", args = c("-NoProfile", "-Command", cmd), stdout = TRUE, stderr = TRUE))
    ok <- file.exists(tmp)
  }

  if (!ok) return(NULL)
  tmp
}

excel_sheets_safe <- function(path) {
  # Work around Windows file locks (Excel open) by copying to temp on failure.
  stopifnot(is.character(path), length(path) == 1L)
  tryCatch(
    suppressMessages(readxl::excel_sheets(path)),
    error = function(e1) {
      tmp <- copy_to_tempfile(path)
      if (is.null(tmp)) stop(e1)
      suppressMessages(readxl::excel_sheets(tmp))
    }
  )
}

read_excel_safe <- function(path, ...) {
  stopifnot(is.character(path), length(path) == 1L)
  tryCatch(
    suppressMessages(readxl::read_excel(path, ...)),
    error = function(e1) {
      tmp <- copy_to_tempfile(path)
      if (is.null(tmp)) stop(e1)
      suppressMessages(readxl::read_excel(tmp, ...))
    }
  )
}

normalize_hyphens <- function(x) {
  # Replace common Unicode dash/hyphen variants with ASCII hyphen-minus.
  if (is.null(x)) return(x)
  x <- as.character(x)
  x <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", x, perl = TRUE)
  x
}

standardize_key <- function(x) {
  # Canonical text key for matching (case-insensitive, whitespace/hyphen normalized).
  x <- normalize_hyphens(x)
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  stringr::str_to_lower(x)
}

canonical_strain_id <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x <- stringr::str_to_upper(x)
  x[x == "6B"] <- "Y6B"
  x
}

parse_strain_time <- function(sample_h) {
  # Parse strings like "Y1-12" or "6B-0" into strain/time_h.
  sample_h <- stringr::str_trim(as.character(sample_h))
  sample_h <- normalize_hyphens(sample_h)
  m <- stringr::str_match(sample_h, "^(.+?)-([0-9]+)$")
  tibble::tibble(
    strain_raw = m[, 2],
    time_h = suppressWarnings(as.integer(m[, 3]))
  ) |>
    dplyr::mutate(strain = canonical_strain_id(.data$strain_raw))
}

is_non_detect <- function(x) {
  # v1 rule: numeric 0 means "not detected" (non-detect) for concentration-like measures.
  # Do NOT apply this to ordinal scores (agar-score) or categorical variables.
  suppressWarnings(!is.na(as.numeric(x)) & as.numeric(x) == 0)
}

as_numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}

fmt_pct <- function(x, digits = 0) {
  scales::percent(x, accuracy = 10^(-digits))
}

dt_num_fixed_js <- function(digits = 4L) {
  # JavaScript renderer for DataTables numeric columns (fixed decimals).
  # - shows exactly `digits` decimals
  # - avoids "-0.0000" artifacts by snapping near-zero to 0
  digits <- as.integer(digits)
  stopifnot(length(digits) == 1L, is.finite(digits), digits >= 0L)

  paste0(
    "function(data, type, row, meta) {",
    "  if (data === null || data === undefined || data === '') return data;",
    "  var num = parseFloat(data);",
    "  if (!isFinite(num)) return data;",
    "  if (type === 'display' || type === 'filter') {",
    "    var D = ", digits, ";",
    "    if (Math.abs(num) < 0.5 * Math.pow(10, -D)) num = 0;",
    "    return num.toFixed(D);",
    "  }",
    "  return num;",
    "}"
  )
}

dt_pval_fixed_js <- function(digits = 4L) {
  # JavaScript renderer for DataTables p-value columns.
  # - shows p-values with fixed decimals
  # - prints very small non-zero p-values as "<0.0001" (threshold depends on digits)
  digits <- as.integer(digits)
  stopifnot(length(digits) == 1L, is.finite(digits), digits >= 0L)

  paste0(
    "function(data, type, row, meta) {",
    "  if (data === null || data === undefined || data === '') return data;",
    "  var num = parseFloat(data);",
    "  if (!isFinite(num)) return data;",
    "  if (type === 'display' || type === 'filter') {",
    "    var D = ", digits, ";",
    "    var thr = Math.pow(10, -D);",
    "    if (num > 0 && num < thr) return '<' + thr.toFixed(D);",
    "    return num.toFixed(D);",
    "  }",
    "  return num;",
    "}"
  )
}
