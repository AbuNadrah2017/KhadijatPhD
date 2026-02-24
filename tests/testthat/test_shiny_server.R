testthat::test_that("Shiny server reactives run for all tabs without error", {
  kh_source("R/utils.R")
  old_wd <- getwd()
  setwd(kh_app_dir())
  on.exit(setwd(old_wd), add = TRUE)

  # Load app into a dedicated environment so we can access ui/server objects.
  app_env <- new.env(parent = globalenv())
  sys.source("app.R", envir = app_env)

  testthat::expect_true(exists("server", envir = app_env))
  testthat::expect_true(exists("ui", envir = app_env))

  shiny::testServer(app_env$server, {
    # Trigger data load once.
    d <- data_all()
    testthat::expect_true(is.list(d))
    testthat::expect_true(nrow(d$metabolites) > 0)
    testthat::expect_true(nrow(d$tol_agar) > 0)
    testthat::expect_true(nrow(d$tol_broth) > 0)
    testthat::expect_true(nrow(d$enzymes) > 0)

    # Objective 1
    session$setInputs(growth_endpoint = "met_growth", growth_show_points = TRUE)
    session$flushReact()
    testthat::expect_s3_class(growth_plot_obj(), "ggplot")
    testthat::expect_true(nrow(growth_table_df()) > 0)
    testthat::expect_s3_class(growth_bar_plot_obj(), "ggplot")

    # Snapshot + delta plots require explicit time selections.
    tp <- sort(unique(d$metabolites$time_h))
    session$setInputs(growth_snapshot_time = tp[[1]])
    session$setInputs(growth_delta_t1 = tp[[1]], growth_delta_t2 = tp[[length(tp)]])
    session$flushReact()
    testthat::expect_s3_class(growth_snapshot_plot_obj(), "ggplot")
    testthat::expect_s3_class(growth_delta_plot_obj(), "ggplot")
    testthat::expect_s3_class(growth_heatmap_plot_obj(), "ggplot")

    # Objective 1 (growth models; descriptive/exploratory)
    session$setInputs(
      growth_model_models = c("logistic", "gompertz"),
      growth_model_display = "best",
      growth_model_ncol = 2
    )
    session$flushReact()
    testthat::expect_s3_class(growth_model_plot_obj(), "ggplot")
    testthat::expect_s3_class(growth_model_resid_plot_obj(), "ggplot")
    testthat::expect_true(nrow(growth_model_table_df()) > 0)

    # Optional: p-values are available (technical duplicates; exploratory).
    session$setInputs(growth_enable_pvals = TRUE, growth_pval_method = "aov_tukey", growth_pval_view = "by_time")
    session$flushReact()
    testthat::expect_true(is.data.frame(growth_pvals_filtered()))

    session$setInputs(growth_endpoint = "aut_survival")
    session$flushReact()
    testthat::expect_s3_class(growth_plot_obj(), "ggplot")

    session$setInputs(growth_endpoint = "aut_od")
    session$flushReact()
    testthat::expect_s3_class(growth_plot_obj(), "ggplot")

    # Objective 2
    session$setInputs(tol_view = "agar")
    session$flushReact()
    a_choice <- sort(unique(d$tol_agar$analyte))[1]
    session$setInputs(tol_analyte = a_choice)
    session$flushReact()
    testthat::expect_s3_class(tol_plot_obj(), "ggplot")
    testthat::expect_true(nrow(tol_table_df()) > 0)

    session$setInputs(tol_view = "broth")
    session$flushReact()
    b_choice <- sort(unique(d$tol_broth$analyte))[1]
    session$setInputs(tol_analyte = b_choice)
    session$flushReact()
    testthat::expect_s3_class(tol_plot_obj(), "ggplot")
    testthat::expect_true(nrow(tol_table_df()) > 0)
    testthat::expect_s3_class(tol_dose_plot_obj(), "ggplot")

    # Objective 2 (composite tolerance index; broth-only; exploratory)
    session$setInputs(tol_index_preset = "balanced")
    session$flushReact()
    testthat::expect_s3_class(tol_index_plot_obj(), "ggplot")
    testthat::expect_true(nrow(tol_index_rank_df()) > 0)
    testthat::expect_true(nrow(tol_index_weights_tbl()) > 0)
    testthat::expect_true(nrow(tol_index_sensitivity_df()) > 0)

    # Optional: broth p-values (technical duplicates; exploratory).
    session$setInputs(tol_enable_pvals = TRUE, tol_pval_method = "aov_tukey")
    session$flushReact()
    testthat::expect_true(is.data.frame(tol_pvals_filtered()))

    # Objective 3
    m_choice <- d$enzymes |>
      dplyr::filter(category == "Enzyme production rate") |>
      dplyr::distinct(measure) |>
      dplyr::arrange(measure) |>
      dplyr::slice(1) |>
      dplyr::pull(measure)
    session$setInputs(enz_measure = m_choice, enz_rep_handling = "mean", enz_heatmap_scale = "z")
    session$flushReact()
    testthat::expect_s3_class(enz_plot_obj(), "ggplot")
    testthat::expect_s3_class(enz_heatmap_plot_obj(), "ggplot")
    testthat::expect_true(nrow(enz_table_df()) > 0)

    session$setInputs(enz_enable_pvals = TRUE, enz_pval_method = "aov_tukey", enz_pval_scope = "selected")
    session$flushReact()
    testthat::expect_true(is.data.frame(enz_pvals_filtered()))

    # Objective 4 (GC)
    if (nrow(d$gc) > 0) {
      session$setInputs(gc_dup_policy = "sum", gc_top_n = 20)
      session$flushReact()
      testthat::expect_s3_class(gc_heatmap_plot(), "ggplot")
      testthat::expect_s3_class(gc_pca_plot(), "ggplot")
      testthat::expect_s3_class(gc_compare_plot_obj(), "ggplot")
      testthat::expect_s3_class(gc_snapshot_plot_obj(), "ggplot")
      testthat::expect_s3_class(gc_total_plot_obj(), "ggplot")
      testthat::expect_s3_class(gc_richness_plot_obj(), "ggplot")
      testthat::expect_true(nrow(gc_table_df()) > 0)
      testthat::expect_true(nrow(gc_pca_loadings_df()) > 0)

      session$setInputs(gc_dup_policy = "mean")
      session$flushReact()
      testthat::expect_s3_class(gc_heatmap_plot(), "ggplot")
      testthat::expect_s3_class(gc_pca_plot(), "ggplot")

      session$setInputs(gc_dup_policy = "keep_distinct")
      session$flushReact()
      testthat::expect_s3_class(gc_heatmap_plot(), "ggplot")
      testthat::expect_s3_class(gc_pca_plot(), "ggplot")
    }

    # Objective 5
    sel <- d$metabolites |>
      dplyr::distinct(type_unit, analyte) |>
      dplyr::arrange(type_unit, analyte)
    met_choice <- paste0(sel$type_unit[1], " | ", sel$analyte[1])
    session$setInputs(met_analyte = met_choice, met_rep_handling = "mean", met_show_points = FALSE)
    session$flushReact()
    testthat::expect_s3_class(met_plot_obj(), "ggplot")
    testthat::expect_true(nrow(met_table_df()) > 0)

    tp_met <- sort(unique(d$metabolites$time_h))
    session$setInputs(met_snapshot_time = tp_met[[1]], met_delta_t1 = tp_met[[1]], met_delta_t2 = tp_met[[length(tp_met)]], met_compare_geom = "dot")
    session$flushReact()
    testthat::expect_s3_class(met_snapshot_plot_obj(), "ggplot")
    testthat::expect_s3_class(met_delta_plot_obj(), "ggplot")
    testthat::expect_s3_class(met_heatmap_plot_obj(), "ggplot")

    session$setInputs(met_enable_pvals = TRUE, met_pval_method = "aov_tukey", met_pval_view = "by_time")
    session$flushReact()
    testthat::expect_true(is.data.frame(met_pvals_filtered()))

    # Objective 5 (PCA; exploratory)
    session$setInputs(
      met_pca_transform = "log1p",
      met_pca_features = "top",
      met_pca_top_n = 15,
      met_pca_scale = TRUE,
      met_pca_color_by = "strain",
      met_pca_zero_fill = TRUE
    )
    session$flushReact()
    testthat::expect_s3_class(met_pca_plot_obj(), "ggplot")
    testthat::expect_true(nrow(met_pca_loadings_df()) > 0)

    # Objective 5 (yield/efficiency; descriptive)
    feats <- d$metabolites |>
      dplyr::distinct(type_unit, analyte) |>
      dplyr::arrange(type_unit, analyte)
    choices <- paste0(feats$type_unit, " | ", feats$analyte)
    idx_sugar <- which(grepl("sugar", feats$type_unit, ignore.case = TRUE))
    idx_acid <- which(grepl("acid", feats$type_unit, ignore.case = TRUE))
    sub_choice <- if (length(idx_sugar) >= 1) choices[[idx_sugar[1]]] else choices[[1]]
    prod_choice <- if (length(idx_acid) >= 1) choices[[idx_acid[1]]] else choices[[min(2, length(choices))]]
    session$setInputs(
      yield_substrate = sub_choice,
      yield_product = prod_choice,
      yield_t1 = tp_met[[1]],
      yield_t2 = tp_met[[length(tp_met)]]
    )
    session$flushReact()
    testthat::expect_true(is.data.frame(met_yield_df()))
    p_yield <- met_yield_plot_obj()
    testthat::expect_true(inherits(p_yield, "patchwork") || inherits(p_yield, "ggplot"))

    # Objective 5 (alternate sources)
    if (nrow(d$autolytic) > 0) {
      a_choice2 <- sort(unique(d$autolytic$analyte))[1]
      session$setInputs(met_source = "aut", met_analyte = a_choice2)
      session$flushReact()
      testthat::expect_s3_class(met_plot_obj(), "ggplot")
      testthat::expect_s3_class(met_heatmap_plot_obj(), "ggplot")
      testthat::expect_s3_class(met_pca_plot_obj(), "ggplot")
    }
    if (nrow(d$reducing) > 0) {
      r_choice <- sort(unique(d$reducing$analyte))[1]
      session$setInputs(met_source = "red24", met_analyte = r_choice)
      session$flushReact()
      testthat::expect_s3_class(met_plot_obj(), "ggplot")
      testthat::expect_s3_class(met_heatmap_plot_obj(), "ggplot")
      # PCA should return a "not available" plot (still a ggplot object).
      testthat::expect_s3_class(met_pca_plot_obj(), "ggplot")
    }

    # Objective 6
    session$setInputs(w_growth = 0.4, w_tolerance = 0.2, w_enz = 0.2, w_vol = 0.2)
    session$flushReact()
    testthat::expect_true(nrow(rank_df()) > 0)
    testthat::expect_true(nrow(rank_table_df()) > 0)
    testthat::expect_true(nrow(rank_components_df()) > 0)
    testthat::expect_true(nrow(rank_weights_tbl()) > 0)
    testthat::expect_s3_class(rank_plot_obj(), "ggplot")
    testthat::expect_s3_class(rank_contrib_plot_obj(), "ggplot")
    testthat::expect_true(nrow(rank_sensitivity_df()) > 0)
    testthat::expect_s3_class(rank_sensitivity_plot_obj(), "ggplot")

    session$setInputs(rank_trade_x = "growth_max", rank_trade_y = "vol_total_mean", rank_trade_scaled = TRUE)
    session$flushReact()
    testthat::expect_s3_class(rank_trade_plot_obj(), "ggplot")

    # QC
    session$setInputs(qc_dataset = "met")
    session$flushReact()
    qc_choice <- sort(unique(d$metabolites$analyte))[1]
    session$setInputs(qc_analyte = qc_choice)
    session$flushReact()
    testthat::expect_s3_class(qc_plot_obj(), "ggplot")
  })
})
