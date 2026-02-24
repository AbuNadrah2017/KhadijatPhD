## Data loading + tidying for the KhadijatPhD datasets.

load_metabolites <- function(path) {
  raw <- read_excel_safe(path, sheet = "Yeast metabolites", col_names = TRUE)
  stopifnot(all(c("Sample-Name", "Type-unit", "Analyte-Acid", "sample-h", "S-1", "S-dup") %in% names(raw)))

  df0 <- raw |>
    dplyr::mutate(
      strain = canonical_strain_id(.data[["Sample-Name"]]),
      analyte = as.character(.data[["Analyte-Acid"]]),
      type_unit = as.character(.data[["Type-unit"]]),
      sample_h = as.character(.data[["sample-h"]])
    )

  parsed <- parse_strain_time(df0$sample_h) |>
    dplyr::select(time_h, strain_from_sample_h = strain)

  df <- dplyr::bind_cols(df0, parsed) |>
    dplyr::mutate(
      strain = dplyr::coalesce(.data$strain, .data$strain_from_sample_h),
      strain = canonical_strain_id(.data$strain),
      analyte = stringr::str_trim(.data$analyte),
      analyte_key = standardize_key(.data$analyte)
    ) |>
    tidyr::pivot_longer(
      cols = c("S-1", "S-dup"),
      names_to = "tech_rep",
      values_to = "value_raw"
    ) |>
    dplyr::mutate(
      value_raw = as_numeric_or_na(.data$value_raw),
      detected = !is.na(.data$value_raw) & .data$value_raw > 0,
      value = dplyr::if_else(.data$value_raw == 0, NA_real_, .data$value_raw)
    )

  df
}

load_reducing_sugar <- function(path) {
  raw <- read_excel_safe(path, sheet = "Reducing sugar-24h", col_names = TRUE)
  stopifnot(all(c("Analyte", "sample", "S-1", "S-dup") %in% names(raw)))

  df <- raw |>
    dplyr::mutate(
      strain = canonical_strain_id(.data[["sample"]]),
      analyte = as.character(.data[["Analyte"]]),
      analyte_key = standardize_key(.data$analyte)
    ) |>
    tidyr::pivot_longer(cols = c("S-1", "S-dup"), names_to = "tech_rep", values_to = "value_raw") |>
    dplyr::mutate(
      value_raw = as_numeric_or_na(.data$value_raw),
      detected = !is.na(.data$value_raw) & .data$value_raw > 0,
      value = dplyr::if_else(.data$value_raw == 0, NA_real_, .data$value_raw),
      time_h = 24L
    )

  df
}

load_tolerance_agar <- function(path) {
  raw <- read_excel_safe(path, sheet = "tolerance test-agar", col_names = TRUE)
  stopifnot(all(c("sample", "analyte", "concentration", "agar-score") %in% names(raw)))

  raw |>
    dplyr::mutate(
      strain = canonical_strain_id(.data[["sample"]]),
      analyte = stringr::str_trim(as.character(.data[["analyte"]])),
      concentration = stringr::str_trim(as.character(.data[["concentration"]])),
      agar_score = as_numeric_or_na(.data[["agar-score"]]),
      analyte_key = standardize_key(.data$analyte),
      condition_key = standardize_key(paste(.data$analyte, .data$concentration, sep = " | "))
    )
}

load_tolerance_broth <- function(path) {
  raw <- read_excel_safe(path, sheet = "tolerance test-broth", col_names = TRUE)
  stopifnot(all(c("sample", "analyte", "concentration", "S-1", "S-dup") %in% names(raw)))

  raw |>
    dplyr::mutate(
      strain = canonical_strain_id(.data[["sample"]]),
      analyte = stringr::str_trim(as.character(.data[["analyte"]])),
      concentration = stringr::str_trim(as.character(.data[["concentration"]])),
      analyte_key = standardize_key(.data$analyte),
      condition_key = standardize_key(paste(.data$analyte, .data$concentration, sep = " | "))
    ) |>
    tidyr::pivot_longer(cols = c("S-1", "S-dup"), names_to = "tech_rep", values_to = "value_raw") |>
    dplyr::mutate(
      value_raw = as_numeric_or_na(.data$value_raw),
      # For broth tolerance, 0 means "no growth" (real 0), not a non-detect.
      value = .data$value_raw
    )
}

load_autolytic <- function(path) {
  raw <- read_excel_safe(path, sheet = "Autolytic-full", col_names = TRUE)
  stopifnot(all(c("Sample-Name", "Analyte-Acid", "sample-h", "S-1", "S-dup") %in% names(raw)))

  df0 <- raw |>
    dplyr::mutate(
      strain = canonical_strain_id(.data[["Sample-Name"]]),
      analyte = stringr::str_trim(as.character(.data[["Analyte-Acid"]])),
      sample_h = as.character(.data[["sample-h"]]),
      analyte_key = standardize_key(.data$analyte)
    )

  parsed <- parse_strain_time(df0$sample_h) |>
    dplyr::select(time_h)

  dplyr::bind_cols(df0, parsed) |>
    tidyr::pivot_longer(cols = c("S-1", "S-dup"), names_to = "tech_rep", values_to = "value_raw") |>
    dplyr::mutate(
      value_raw = as_numeric_or_na(.data$value_raw),
      detected = !is.na(.data$value_raw) & .data$value_raw > 0,
      value = dplyr::if_else(.data$value_raw == 0, NA_real_, .data$value_raw)
    )
}

load_enzymes <- function(path) {
  raw <- read_excel_safe(path, sheet = "Enzyme production", col_names = TRUE)
  stopifnot(all(c("sample", "analyte-Log(CFU/g )", "concentration", "S-1", "S-dup") %in% names(raw)))

  raw |>
    dplyr::mutate(
      strain = canonical_strain_id(.data[["sample"]]),
      category = stringr::str_trim(as.character(.data[["analyte-Log(CFU/g )"]])),
      measure = stringr::str_trim(as.character(.data[["concentration"]])),
      category_key = standardize_key(.data$category),
      measure_key = standardize_key(.data$measure)
    ) |>
    tidyr::pivot_longer(cols = c("S-1", "S-dup"), names_to = "tech_rep", values_to = "value_raw") |>
    dplyr::mutate(
      value_raw = as_numeric_or_na(.data$value_raw),
      detected = !is.na(.data$value_raw) & .data$value_raw > 0,
      value = dplyr::if_else(.data$value_raw == 0, NA_real_, .data$value_raw)
    )
}

load_compound_map <- function(path) {
  if (!file.exists(path)) return(tibble::tibble())
  raw <- read_excel_safe(path, sheet = "compound_map", col_names = TRUE)
  if (!all(c("Compound class", "Common (sensory) name") %in% names(raw))) return(tibble::tibble())

  raw |>
    dplyr::transmute(
      compound_class = stringr::str_trim(as.character(.data[["Compound class"]])),
      common_name = stringr::str_trim(as.character(.data[["Common (sensory) name"]])),
      sensory_note = stringr::str_trim(as.character(.data[["Sensory note / comment"]] %||% NA_character_)),
      compound_std = vapply(.data$common_name, standardize_compound_name, character(1)),
      compound_key = standardize_key(.data$compound_std)
    ) |>
    dplyr::filter(!is.na(.data$compound_key) & .data$compound_key != "")
}

load_gc_volatiles <- function(path) {
  sheets <- excel_sheets_safe(path)
  if (length(sheets) == 0) return(tibble::tibble())

  # Prefer strain sheets; keep QC separately for optional background inspection.
  strain_sheets <- sheets[!stringr::str_detect(stringr::str_to_lower(sheets), "^qc")]

  dfs <- lapply(strain_sheets, function(sh) gc_tidy_from_sheet(path, sheet = sh))
  dplyr::bind_rows(dfs) |>
    dplyr::mutate(strain = canonical_strain_id(.data$strain))
}
