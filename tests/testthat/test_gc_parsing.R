testthat::test_that("GC volatiles parse and duplicate policy behaves", {
  kh_source("R/utils.R")
  kh_source("R/gc_parse.R")
  kh_source("R/load_data.R")

  gc <- load_gc_volatiles(kh_data_path("GC_VolatilesYeast.xlsx"))
  testthat::expect_true(nrow(gc) > 0)
  testthat::expect_true(all(c("strain", "time_h", "compound_key", "ra_ug_L") %in% names(gc)))

  # Canonical strain IDs: no legacy 6B after load.
  testthat::expect_false(any(gc$strain == "6B", na.rm = TRUE))
  testthat::expect_true(any(gc$strain == "Y6B", na.rm = TRUE))

  # Duplicate handling: sum/mean should collapse to unique compound per strain/time/compound_key.
  g_sum <- gc_apply_duplicate_policy(gc, "sum")
  cnt_sum <- dplyr::count(g_sum, strain, time_h, compound_key, name = "n")
  testthat::expect_true(all(cnt_sum$n == 1L))

  g_mean <- gc_apply_duplicate_policy(gc, "mean")
  cnt_mean <- dplyr::count(g_mean, strain, time_h, compound_key, name = "n")
  testthat::expect_true(all(cnt_mean$n == 1L))

  # keep_distinct should retain duplicates by adding deterministic suffix when needed.
  g_keep <- gc_apply_duplicate_policy(gc, "keep_distinct")
  # If no duplicates exist, this would be 0; in practice the file contains duplicates.
  testthat::expect_true(any(stringr::str_detect(g_keep$compound_key, "#"), na.rm = TRUE))
})
