## GC volatiles parsing and standardization

gc_extract_time_from_header <- function(header) {
  header <- as.character(header %||% "")
  header <- normalize_hyphens(header)
  # Try "-72-RA" or "-72-Names" patterns.
  m <- stringr::str_match(header, "-([0-9]{1,3})-(?:ra|names)\\b")
  if (!is.na(m[, 2])) return(as.integer(m[, 2]))
  # Try any "-72-" pattern.
  m2 <- stringr::str_match(header, "-([0-9]{1,3})-")
  if (!is.na(m2[, 2])) return(as.integer(m2[, 2]))
  # Try trailing "-72" pattern.
  m3 <- stringr::str_match(header, "-([0-9]{1,3})\\b")
  if (!is.na(m3[, 2])) return(as.integer(m3[, 2]))
  NA_integer_
}

gc_is_ra_col <- function(name) {
  stringr::str_detect(stringr::str_to_lower(as.character(name)), "\\bra\\b|ra\\s*\\(|ra_")
}

gc_is_names_col <- function(name) {
  stringr::str_detect(stringr::str_to_lower(as.character(name)), "name")
}

standardize_compound_name <- function(x) {
  # Controlled standardization for matching + de-duplication across sheets.
  # Keep this conservative; do not invent new compound identities.
  x <- as.character(x)
  x <- normalize_hyphens(x)
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)

  # Normalize common typos/synonyms observed in the GC sheets.
  # Keys must use standardize_key() so matching is consistent.
  key <- standardize_key(x)
  fixes <- c(
    "2-ocatanone" = "2-Octanone",
    "isoamylalcohol" = "Isoamyl alcohol",
    "ethyl acetate" = "Ethyl acetate",
    "isoamyl alc..." = "Isoamyl alcohol" # safety for truncated display; rarely used
  )
  if (key %in% names(fixes)) return(fixes[[key]])

  x
}

gc_tidy_from_sheet <- function(path, sheet, strain_override = NULL) {
  raw <- read_excel_safe(path, sheet = sheet, col_names = TRUE)
  if (nrow(raw) == 0) return(tibble::tibble())

  strain <- strain_override %||% sheet
  strain <- stringr::str_replace(strain, "-[Vv]olatiles$", "")
  strain <- canonical_strain_id(strain)

  cn <- names(raw)
  blocks <- list()

  i <- 1L
  while (i <= length(cn) - 1L) {
    ra_ok <- gc_is_ra_col(cn[[i]])
    names_ok <- gc_is_names_col(cn[[i + 1L]])
    if (!ra_ok || !names_ok) {
      i <- i + 1L
      next
    }

    time_h <- gc_extract_time_from_header(cn[[i]])
    if (is.na(time_h)) {
      # Some RA headers miss the time; try the paired Names header.
      time_h <- gc_extract_time_from_header(cn[[i + 1L]])
    }

    ra <- raw[[i]]
    nm <- raw[[i + 1L]]

    block <- tibble::tibble(
      strain = strain,
      time_h = time_h,
      compound_raw = as.character(nm),
      ra_ug_L_raw = as_numeric_or_na(ra)
    ) |>
      dplyr::filter(!is.na(.data$compound_raw) & .data$compound_raw != "") |>
      dplyr::mutate(
        compound_std = vapply(.data$compound_raw, standardize_compound_name, character(1)),
        compound_key = standardize_key(.data$compound_std),
        detected = !is.na(.data$ra_ug_L_raw) & .data$ra_ug_L_raw > 0,
        ra_ug_L = dplyr::if_else(.data$ra_ug_L_raw == 0, NA_real_, .data$ra_ug_L_raw)
      )

    blocks[[length(blocks) + 1L]] <- block
    i <- i + 2L
  }

  dplyr::bind_rows(blocks) |>
    dplyr::mutate(
      time_h = as.integer(.data$time_h),
      sheet = sheet
    ) |>
    dplyr::filter(!is.na(.data$time_h))
}

  gc_apply_duplicate_policy <- function(df, policy = c("sum", "mean", "keep_distinct")) {
    policy <- match.arg(policy)
    if (nrow(df) == 0) return(df)

  df <- df |>
    dplyr::arrange(.data$strain, .data$time_h, .data$compound_key, dplyr::desc(.data$detected))

  if (policy == "keep_distinct") {
    df |>
      dplyr::group_by(.data$strain, .data$time_h, .data$compound_key) |>
      dplyr::mutate(
        dup_index = dplyr::row_number(),
        compound_key2 = if (dplyr::n() > 1L) paste0(.data$compound_key, "#", .data$dup_index) else .data$compound_key,
        compound_std2 = if (dplyr::n() > 1L) paste0(.data$compound_std, " #", .data$dup_index) else .data$compound_std
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(compound_key = .data$compound_key2, compound_std = .data$compound_std2) |>
      dplyr::select(-dplyr::any_of(c("dup_index", "compound_key2", "compound_std2")))
    } else {
      fun <- switch(policy, sum = sum, mean = mean)
      df |>
        # Group by compound_key only so case variants (e.g., Ethanol vs ethanol) collapse.
        dplyr::group_by(.data$strain, .data$time_h, .data$compound_key) |>
        dplyr::summarise(
          # Avoid NaN from mean(numeric(0)) when all values are NA (non-detects).
          ra_ug_L = if (all(is.na(.data$ra_ug_L))) NA_real_ else fun(.data$ra_ug_L, na.rm = TRUE),
          ra_ug_L_raw = if (all(is.na(.data$ra_ug_L_raw))) NA_real_ else fun(.data$ra_ug_L_raw, na.rm = TRUE),
          detected = any(.data$detected, na.rm = TRUE),
          compound_std = dplyr::first(.data$compound_std),
          compound_raw = dplyr::first(.data$compound_raw),
          sheet = dplyr::first(.data$sheet),
          .groups = "drop"
        ) |>
      dplyr::mutate(
        ra_ug_L = dplyr::if_else(is.infinite(.data$ra_ug_L), NA_real_, .data$ra_ug_L),
        ra_ug_L_raw = dplyr::if_else(is.infinite(.data$ra_ug_L_raw), NA_real_, .data$ra_ug_L_raw)
      )
  }
}
