testthat::test_that("All core datasets load and canonicalization rules hold", {
  kh_source("R/utils.R")
  kh_source("R/gc_parse.R")
  kh_source("R/load_data.R")

  met <- load_metabolites(kh_data_path("Metabolites.xlsx"))
  red <- load_reducing_sugar(kh_data_path("Metabolites.xlsx"))
  tol_a <- load_tolerance_agar(kh_data_path("Yeast tolerance test.xlsx"))
  tol_b <- load_tolerance_broth(kh_data_path("Yeast tolerance test.xlsx"))
  aut <- load_autolytic(kh_data_path("Yeast tolerance test.xlsx"))
  enz <- load_enzymes(kh_data_path("Yeast tolerance test.xlsx"))
  cmap <- load_compound_map(kh_data_path("exampleMapData.xlsx"))

  testthat::expect_true(nrow(met) > 0)
  testthat::expect_true(nrow(red) > 0)
  testthat::expect_true(nrow(tol_a) > 0)
  testthat::expect_true(nrow(tol_b) > 0)
  testthat::expect_true(nrow(aut) > 0)
  testthat::expect_true(nrow(enz) > 0)
  testthat::expect_true(nrow(cmap) > 0)

  # Canonical strain IDs: 6B should be normalized to Y6B everywhere.
  testthat::expect_false(any(met$strain == "6B", na.rm = TRUE))
  testthat::expect_false(any(tol_a$strain == "6B", na.rm = TRUE))
  testthat::expect_false(any(tol_b$strain == "6B", na.rm = TRUE))
  testthat::expect_false(any(aut$strain == "6B", na.rm = TRUE))
  testthat::expect_false(any(enz$strain == "6B", na.rm = TRUE))

  # Metabolites.xlsx contains legacy "6B" in the raw file; ensure it appears as Y6B after load.
  testthat::expect_true(any(met$strain == "Y6B", na.rm = TRUE))

  # Non-detect rule: numeric 0 must map to NA in value (but remain in value_raw).
  testthat::expect_true(any(met$value_raw == 0, na.rm = TRUE))
  testthat::expect_true(any(is.na(met$value) & met$value_raw == 0, na.rm = TRUE))

  # Agar score is ordinal; zeros are allowed and must not be coerced to NA.
  testthat::expect_true(any(tol_a$agar_score == 0, na.rm = TRUE))
  testthat::expect_false(any(is.na(tol_a$agar_score) & tol_a$agar_score == 0, na.rm = TRUE))

  # Broth tolerance: 0 means "no growth" (real 0), not missing.
  testthat::expect_true(any(tol_b$value_raw == 0, na.rm = TRUE))
  testthat::expect_true(any(tol_b$value == 0, na.rm = TRUE))
  testthat::expect_false(any(is.na(tol_b$value) & tol_b$value_raw == 0, na.rm = TRUE))
})

testthat::test_that("Hyphen/case standardization behaves as specified", {
  kh_source("R/utils.R")

  testthat::expect_identical(standardize_key("Ethanol"), standardize_key("ethanol"))
  testthat::expect_identical(standardize_key("1\u2011Heptanol"), standardize_key("1-Heptanol")) # non-breaking hyphen -> '-'
  testthat::expect_identical(canonical_strain_id("6B"), "Y6B")
  testthat::expect_identical(canonical_strain_id("Y6B"), "Y6B")
})
