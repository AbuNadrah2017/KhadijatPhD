kh_find_root <- function(start_dir = getwd(), max_up = 12L) {
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

kh_root <- function() {
  root <- getOption("khadijat_phd_root")
  if (is.character(root) && length(root) == 1L && file.exists(file.path(root, "AGENTS.md")) &&
    (dir.exists(file.path(root, "Data")) || file.exists(file.path(root, "R", "app.R")))) {
    return(normalizePath(root, winslash = "/", mustWork = FALSE))
  }
  kh_find_root(getwd())
}

kh_path <- function(...) file.path(kh_root(), ...)

kh_app_dir <- function() {
  root <- kh_root()
  # Current layout: app + data live under `R/`.
  if (file.exists(file.path(root, "R", "app.R"))) return(file.path(root, "R"))
  # Legacy layout: app at repo root.
  root
}

kh_data_path <- function(filename) {
  # Support both the legacy `Data/` layout and the current layout where xlsx files live under `R/`.
  root <- kh_root()
  candidates <- c(
    file.path(root, "Data", filename),
    file.path(root, "R", filename),
    file.path(root, filename)
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) return(existing[[1]])
  candidates[[2]]
}

kh_source <- function(rel_path) {
  source(kh_path(rel_path), local = parent.frame())
}
