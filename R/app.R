#KhadijatPhD
## KhadijatPhD: Internal-data Shiny app for peer-review-safe exploratory analyses.
##
## Key constraints (v1):
## - Mono-culture only (no co-culture IDs in current dataset).
## - S-1 / S-dup are technical duplicates of the same sample/run (not biological replicates).
## - 0 means "not detected" for concentration-like measures; do not impute; report detection frequency.
## - No inferential p-values/CIs presented as evidence of between-strain differences for this dataset.
## - GC volatiles have no replicates; descriptive/exploratory only.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(DT)
  library(scales)
  library(patchwork)
  library(jsonlite)
  library(writexl)
})

source("utils.R", local = TRUE)
source("gc_parse.R", local = TRUE)
source("load_data.R", local = TRUE)
source("stats.R", local = TRUE)
source("plots_tables.R", local = TRUE)

dt_coldefs_numeric <- function(df, digits = 4L) {
  # Apply publication-friendly numeric formatting to DT tables:
  # - fixed decimals for continuous numeric columns
  # - p-values get "<0.0001" style display
  # - integer columns (timepoints, counts) are left as-is
  if (!is.data.frame(df) || ncol(df) == 0) return(list())

  is_num <- vapply(df, is.numeric, logical(1))
  is_int <- vapply(df, is.integer, logical(1))
  is_fmt <- is_num & !is_int
  if (!any(is_fmt)) return(list())

  p_cols <- grepl("(^p($|_)|p_adj|global_p)", names(df), ignore.case = TRUE)
  p_targets <- unname(which(is_fmt & p_cols)) - 1L
  other_targets <- unname(which(is_fmt & !p_cols)) - 1L

  defs <- list()
  if (length(other_targets) > 0) defs[[length(defs) + 1L]] <- list(targets = other_targets, render = JS(dt_num_fixed_js(digits)))
  if (length(p_targets) > 0) defs[[length(defs) + 1L]] <- list(targets = p_targets, render = JS(dt_pval_fixed_js(digits)))
  defs
}

project_root <- getwd()  #find_project_root(getwd())
paths <- list(
  metabolites = file.path(project_root,  "Metabolites.xlsx"),
  tolerance = file.path(project_root,  "Yeast tolerance test.xlsx"),
  gc = file.path(project_root, "GC_VolatilesYeast.xlsx"),
  map = file.path(project_root, "exampleMapData.xlsx") #,
  #agents = file.path(project_root, "AGENTS.md")
)

objective_coverage <- function(d) {
  out <- tibble::tibble(
    objective = c(
      "Objective 1: Growth (mono-culture)",
      "Objective 2: Tolerance (conditions present)",
      "Objective 3: Enzymes (assays present)",
      "Objective 4: GC volatiles (mono-culture)",
      "Objective 5: Metabolic dynamics/associations (mono-culture)",
      "Objective 6: Best yeast (multi-criteria ranking)"
    ),
    status = c("Available", "Available", "Available", "Available", "Available", "Available"),
    notes = c(
      "Co-culture comparisons are deferred (no co-culture IDs in dataset).",
      "Only conditions present in tolerance sheets are analysed; missing items show as Not available.",
      "Only assays present in Enzyme production sheet are analysed; missing items show as Not available.",
      "No replicate-based inference; descriptive/exploratory only.",
      "No synergy/antagonism claims; no co-culture analyses in v1.",
      "Ranking is transparent and sensitivity-checked; still exploratory without biological replication."
    )
  )

  tol_text <- paste(c(d$tol_agar$concentration, d$tol_broth$concentration), collapse = " | ")
  tol_key <- standardize_key(tol_text)
  tol_items <- tibble::tibble(
    item = c("maltose", "citric acid", "fructose"),
    present = c(
      str_detect(tol_key, "maltose"),
      str_detect(tol_key, "citric"),
      str_detect(tol_key, "fructose")
    )
  )

  enz_text <- paste(d$enzymes$measure, collapse = " | ")
  enz_key <- standardize_key(enz_text)
  enz_items <- tibble::tibble(
    item = c("invertase", "cellulase"),
    present = c(
      str_detect(enz_key, "invertase"),
      str_detect(enz_key, "cellulase")
    )
  )

  list(summary = out, tol_items = tol_items, enz_items = enz_items)
}

ui <- page_navbar(
  title = "KhadijatPhD: Yeast Fermentation Analyses (v1)",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    bg = "#ffffff",
    fg = "#1d2733",
    primary = "#0b5cab"
  ),
    header = tagList(
      tags$style(HTML("
        .app-note { font-size: 0.95rem; line-height: 1.25rem; }
        .pill { display:inline-block; padding:0.2rem 0.55rem; border-radius: 999px; background:#eef3fb; border:1px solid #d8e6fb; margin-right:0.35rem; }
        /* DataTables uses floats for controls; ensure containers clear so headings don't overlap. */
        .dataTables_wrapper::after { content: ''; display: block; clear: both; }
        .dataTables_wrapper { width: 100%; overflow-x: auto; }
        /* Additional hardening for the Overview page: force DT blocks to participate in normal flow. */
        .dt-block { display: block; width: 100%; clear: both; }
      "))
    ),

  nav_panel(
    "Overview",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        h4("App rules (peer-review constraints)"),
        div(
          class = "app-note",
          div(class = "pill", "Mono-culture only"),
          #div(class = "pill", "No p-values/CIs as evidence"),
          div(class = "pill", "0 = non-detect (no imputation)"),
          div(class = "pill", "Volatiles: exploratory only"),
          tags$p("The v1 dataset has technical duplicates (S-1/S-dup) but no encoded independent fermentation replicates.")
        ),
        hr(),
        h4("Download"),
        downloadButton("dl_session", "Download session info"),
        downloadButton("dl_params", "Download analysis parameters (JSON)")
      ),
        card(
          card_header("Data status and objective coverage"),
          uiOutput("data_status"),
          hr(),
          h4("Objective coverage (current dataset)"),
          div(class = "dt-block", DTOutput("tbl_coverage")),
          hr(),
          h4("Objective gaps (explicit items)"),
          fluidRow(
            column(6, h5("Tolerance items"), div(class = "dt-block", DTOutput("tbl_tol_items"))),
            column(6, h5("Enzyme items"), div(class = "dt-block", DTOutput("tbl_enz_items")))
          ),
          hr(),
          h4("Notes"),
          tags$ul(
          tags$li("No co-culture comparisons or synergy/antagonism scoring are implemented in v1."),
          tags$li("All plots/tables are computed from the loaded spreadsheets."),
          tags$li("If a condition/assay is absent, the app reports it as Not available in current dataset.")
        )
      )
    )
  ),

  nav_panel(
    "Objective 1: Growth",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        h4("Endpoint"),
        selectInput(
          "growth_endpoint",
          "Choose endpoint",
          choices = c(
            "Microbial growth (log CFU/mL) [Metabolites]" = "met_growth",
            "Survival (log CFU/mL) [Autolytic]" = "aut_survival",
            "OD (AU) [Autolytic]" = "aut_od"
          )
        ),
        h4("Technical duplicates"),
        selectInput(
          "growth_rep_handling",
          "Replicate handling",
          choices = c(
            "Mean(S-1, S-dup) (Recommended)" = "mean",
            "S-1 only" = "s1",
            "S-dup only" = "sdup",
            "Difference (S-1 − S-dup) (QC)" = "diff"
          ),
          selected = "mean"
        ),
        checkboxInput("growth_show_points", "Show technical duplicates (S-1/S-dup)", value = TRUE),
        hr(),
        h4("Downloads"),
        downloadButton("dl_growth_plot", "Download plot (PNG)"),
        downloadButton("dl_growth_snapshot_plot", "Download snapshot (PNG)"),
        downloadButton("dl_growth_model_plot", "Download model fit (PNG)"),
        downloadButton("dl_growth_model_table", "Download model table (XLSX)"),
        downloadButton("dl_growth_pvals", "Download p-values (CSV)"),
        downloadButton("dl_growth_table", "Download table (XLSX)")
      ),
      card(
        card_header("Growth dynamics (mono-culture; descriptive)"),
        div(
          class = "app-note",
          tags$p("Shows technical duplicate measurements (S-1/S-dup) and a user-selected summary (e.g., mean). P-values can be computed from technical duplicates (exploratory only)."),
          tags$p(tags$b("Data sources (file / sheet):")),
          tags$ul(
            tags$li(tags$code("Data/Metabolites.xlsx"), " / ", tags$code("Yeast metabolites"), " (microbial growth log CFU/mL)."),
            tags$li(tags$code("Data/Yeast tolerance test.xlsx"), " / ", tags$code("Autolytic-full"), " (survival log CFU/mL and OD).")
          )
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Time series",
            plotOutput("growth_plot", height = 420),
            hr(),
            h5("Bar chart (optional)"),
            fluidRow(
              column(
                4,
                selectInput(
                  "growth_bar_position",
                  "Bar layout",
                  choices = c("Dodged (recommended)" = "dodge", "Stacked (visual only)" = "stack"),
                  selected = "dodge"
                )
              ),
              column(8, div(class = "app-note", tags$p("Bar chart shows the selected replicate handling by timepoint and strain.")))
            ),
            plotOutput("growth_bar_plot", height = 420)
          ),
          tabPanel(
            "Timepoint compare",
            fluidRow(
              column(4, uiOutput("growth_timepoint_ui")),
              column(4, uiOutput("growth_delta_ui")),
              column(4, selectInput("growth_compare_geom", "Geometry", choices = c("Dot" = "dot", "Bar" = "bar"), selected = "dot"))
            ),
            plotOutput("growth_snapshot_plot", height = 420),
            hr(),
            plotOutput("growth_delta_plot", height = 360)
          ),
          tabPanel(
            "Heatmap",
            plotOutput("growth_heatmap_plot", height = 520)
          ),
          tabPanel(
            "Stats (p-values)",
            div(
              class = "app-note",
              tags$p(tags$b("Important:"), " p-values in v1 are based on technical duplicates (repeat measurements of the same sample/run).")
            ),
            checkboxInput("growth_enable_pvals", "Enable p-values (technical duplicates; exploratory)", value = FALSE),
            conditionalPanel(
              condition = "input.growth_enable_pvals === true",
              fluidRow(
                column(
                  4,
                  selectInput(
                    "growth_pval_method",
                    "Method",
                    choices = c("ANOVA + Tukey" = "aov_tukey", "Kruskal + Wilcoxon" = "kw_wilcox"),
                    selected = "aov_tukey"
                  )
                ),
                column(
                  4,
                  selectInput(
                    "growth_pval_view",
                    "Compare",
                    choices = c("Strains at each timepoint" = "by_time", "Timepoints within each strain" = "by_strain"),
                    selected = "by_time"
                  )
                ),
                column(
                  4,
                  selectInput(
                    "growth_padj_overall",
                    "Overall adjustment (across all pairwise tests shown)",
                    choices = c("BH (FDR)" = "BH", "Holm" = "holm", "Bonferroni" = "bonferroni"),
                    selected = "BH"
                  )
                )
              ),
              uiOutput("growth_pval_filter_ui"),
              DTOutput("growth_pvals_table")
            )
          ),
          tabPanel(
            "Model fit (exploratory)",
            div(
              class = "app-note",
              tags$p("Fits logistic/Gompertz/Baranyi growth curves per strain for descriptive comparison. With the current dataset (technical duplicates only), do not interpret fits/parameters as biological inference.")
            ),
            fluidRow(
              column(
                4,
                checkboxGroupInput(
                  "growth_model_models",
                  "Models",
                  choices = c("Logistic" = "logistic", "Gompertz" = "gompertz", "Baranyi" = "baranyi"),
                  selected = c("logistic", "gompertz", "baranyi")
                )
              ),
              column(
                4,
                selectInput(
                  "growth_model_display",
                  "Display",
                  choices = c("Best model per strain (AIC)" = "best", "Show all selected models" = "all"),
                  selected = "best"
                )
              ),
              column(
                4,
                sliderInput("growth_model_ncol", "Facet columns", min = 1, max = 4, value = 2, step = 1)
              )
            ),
            conditionalPanel(
              condition = "input.growth_rep_handling === 'diff'",
              div(class = "app-note", tags$p(tags$b("Note:"), " model fitting is disabled when replicate handling is set to Difference (QC)."))
            ),
            plotOutput("growth_model_plot", height = 520),
            hr(),
            plotOutput("growth_model_resid_plot", height = 380),
            hr(),
            DTOutput("growth_model_table")
          ),
          tabPanel(
            "Summary table",
            DTOutput("growth_table")
          )
        )
      )
    )
  ),

  nav_panel(
    "Objective 2: Tolerance",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        h4("View"),
        radioButtons("tol_view", "Dataset", choices = c("Agar (ordinal score)" = "agar", "Broth (quantitative)" = "broth")),
        hr(),
        h4("Filter"),
        uiOutput("tol_filter_ui"),
        uiOutput("tol_rep_ui"),
        hr(),
        h4("Downloads"),
        downloadButton("dl_tol_plot", "Download heatmap (PNG)"),
        downloadButton("dl_tol_dose_plot", "Download dose-response (PNG)"),
        downloadButton("dl_tol_index_plot", "Download tolerance index (PNG)"),
        downloadButton("dl_tol_index_table", "Download tolerance index (XLSX)"),
        downloadButton("dl_tol_pvals", "Download p-values (CSV)"),
        downloadButton("dl_tol_table", "Download table (XLSX)")
      ),
      card(
        card_header("Stress tolerance profiles (descriptive)"),
        div(
          class = "app-note",
          tags$p("Heatmaps summarize tolerance by strain and condition. For broth, 0 means no growth (real 0). Technical duplicates (S-1/S-dup) can be summarised (mean) or inspected individually."),
          tags$p(tags$b("Data sources (file / sheet):")),
          tags$ul(
            tags$li(tags$code("Data/Yeast tolerance test.xlsx"), " / ", tags$code("tolerance test-agar"), " (agar-score)."),
            tags$li(tags$code("Data/Yeast tolerance test.xlsx"), " / ", tags$code("tolerance test-broth"), " (quantitative; S-1/S-dup).")
          )
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel("Heatmap", plotOutput("tol_plot", height = 520)),
          tabPanel(
            "Dose-response",
            fluidRow(
              column(
                4,
                selectInput(
                  "tol_dose_geom",
                  "Geometry",
                  choices = c("Line + points" = "line", "Bars" = "bar"),
                  selected = "line"
                )
              ),
              column(
                4,
                selectInput(
                  "tol_bar_position",
                  "Bar layout",
                  choices = c("Dodged (recommended)" = "dodge", "Stacked (visual only)" = "stack"),
                  selected = "dodge"
                )
              ),
              column(4, div(class = "app-note", tags$p("Concentration is parsed to numeric where possible (%, pH, temperature).")))
            ),
            plotOutput("tol_dose_plot", height = 520)
          ),
          tabPanel(
            "Stats (p-values)",
            div(
              class = "app-note",
              tags$p(tags$b("Important:"), " broth p-values are based on technical duplicates (repeat measurements) in v1.")
            ),
            checkboxInput("tol_enable_pvals", "Enable p-values for broth (technical duplicates; exploratory)", value = FALSE),
            conditionalPanel(
              condition = "input.tol_enable_pvals === true && input.tol_view === 'broth'",
              fluidRow(
                column(
                  4,
                  selectInput(
                    "tol_pval_method",
                    "Method",
                    choices = c("ANOVA + Tukey" = "aov_tukey", "Kruskal + Wilcoxon" = "kw_wilcox"),
                    selected = "aov_tukey"
                  )
                ),
                column(
                  4,
                  selectInput(
                    "tol_padj_overall",
                    "Overall adjustment (across all pairwise tests shown)",
                    choices = c("BH (FDR)" = "BH", "Holm" = "holm", "Bonferroni" = "bonferroni"),
                    selected = "BH"
                  )
                ),
                column(4, uiOutput("tol_pval_conc_ui"))
              ),
              DTOutput("tol_pvals_table")
            ),
            conditionalPanel(
              condition = "input.tol_enable_pvals === true && input.tol_view === 'agar'",
              div(class = "app-note", tags$p("Agar tolerance has one score per strain/condition (no replicates), so inferential p-values are not available."))
            )
          ),
          tabPanel("Table", DTOutput("tol_table")),
          tabPanel(
            "Tolerance index (broth; exploratory)",
            div(
              class = "app-note",
              tags$p("Computes a transparent composite tolerance score from broth responses for the selected analyte group. This is descriptive/exploratory (technical duplicates only); do not interpret rankings as biologically replicated inference.")
            ),
            conditionalPanel(
              condition = "input.tol_view === 'agar'",
              div(class = "app-note", tags$p(tags$b("Not available:"), " composite tolerance scoring is implemented for broth only (quantitative)."))
            ),
            conditionalPanel(
              condition = "input.tol_view === 'broth'",
              fluidRow(
                column(
                  3,
                  selectInput(
                    "tol_index_preset",
                    "Preset weights",
                    choices = c(
                      "Balanced (Recommended)" = "balanced",
                      "Equal weights" = "equal",
                      "Mean-focused" = "mean",
                      "Worst-case (min) focused" = "min",
                      "AUC-focused (numeric concentrations)" = "auc",
                      "IC50-focused (numeric concentrations)" = "ic50",
                      "Custom (use sliders)" = "custom"
                    ),
                    selected = "balanced"
                  )
                ),
                column(3, sliderInput("w_tol_mean", "Mean weight", min = 0, max = 1, value = 0.40, step = 0.05)),
                column(3, sliderInput("w_tol_min", "Min weight", min = 0, max = 1, value = 0.30, step = 0.05)),
                column(3, sliderInput("w_tol_auc", "AUC weight", min = 0, max = 1, value = 0.30, step = 0.05))
              ),
              fluidRow(
                column(3, sliderInput("w_tol_ic50", "IC50 weight", min = 0, max = 1, value = 0.00, step = 0.05)),
                column(9, div(class = "app-note", tags$p("AUC/IC50 require numeric concentration parsing within the selected analyte group. If not available, these metrics are excluded and weights are renormalized.")))
              ),
              hr(),
              fluidRow(
                column(6, plotOutput("tol_index_plot", height = 420)),
                column(6, DTOutput("tol_index_table"))
              ),
              hr(),
              h4("Weights used (after normalization)"),
              DTOutput("tol_index_weights"),
              hr(),
              h4("Sensitivity across presets"),
              plotOutput("tol_index_sensitivity_plot", height = 420),
              hr(),
              DTOutput("tol_index_sensitivity_table"),
              hr(),
              h4("Underlying metrics (selected analyte group)"),
              DTOutput("tol_index_metrics_table")
            )
          )
        )
      )
    )
  ),

    nav_panel(
      "Objective 3: Enzymes",
      layout_sidebar(
        sidebar = sidebar(
          width = 360,
          h4("Filter"),
          uiOutput("enz_select_ui"),
          hr(),
          h4("Downloads"),
          downloadButton("dl_enz_plot", "Download plot (PNG)"),
          downloadButton("dl_enz_heatmap", "Download heatmap (PNG)"),
          downloadButton("dl_enz_pvals", "Download p-values (CSV)"),
          downloadButton("dl_enz_table", "Download table (XLSX)")
        ),
        card(
          card_header("Enzymatic activity: production rate"),
          div(
            class = "app-note",
            tags$p("This tab is restricted to the ", tags$b("Enzyme production rate"), " category only."),
            tags$p("Values are summarized from technical duplicates (S-1/S-dup). P-values (if enabled) are technical-duplicate / exploratory."),
            tags$p(tags$b("Data sources (file / sheet):")),
            tags$ul(
              tags$li(tags$code("Data/Yeast tolerance test.xlsx"), " / ", tags$code("Enzyme production"), " (category: Enzyme production rate).")
            )
          ),
          tabsetPanel(
            type = "tabs",
            tabPanel("By measure", plotOutput("enz_plot", height = 520)),
            tabPanel(
              "Heatmap",
              fluidRow(
                column(
                  4,
                  selectInput(
                    "enz_heatmap_scale",
                    "Heatmap scale",
                    choices = c("Raw values" = "raw", "Z-score per measure (recommended)" = "z"),
                    selected = "z"
                  )
                ),
                column(8, div(class = "app-note", tags$p("Z-score scaling standardizes each measure across strains (useful because measures have different units/magnitudes).")))
              ),
              plotOutput("enz_heatmap", height = 560)
            ),
            tabPanel(
              "Stats (p-values)",
              div(
                class = "app-note",
                tags$p(tags$b("Important:"), " p-values are based on technical duplicates (repeat measurements). They are exploratory."),
                tags$p("Time comparisons are not available because the Enzyme production sheet has no time variable; to enable time-based comparisons, add a time column/ID to the dataset.")
              ),
              checkboxInput("enz_enable_pvals", "Enable p-values (technical duplicates; exploratory)", value = FALSE),
              conditionalPanel(
                condition = "input.enz_enable_pvals === true",
                fluidRow(
                  column(
                    4,
                    selectInput(
                      "enz_pval_method",
                      "Method",
                      choices = c("ANOVA + Tukey" = "aov_tukey", "Kruskal + Wilcoxon" = "kw_wilcox"),
                      selected = "aov_tukey"
                    )
                  ),
                  column(
                    4,
                    selectInput(
                      "enz_padj_overall",
                      "Overall adjustment (across all pairwise tests shown)",
                      choices = c("BH (FDR)" = "BH", "Holm" = "holm", "Bonferroni" = "bonferroni"),
                      selected = "BH"
                    )
                  ),
                  column(
                    4,
                    radioButtons(
                      "enz_pval_scope",
                      "Scope",
                      choices = c("Selected measure" = "selected", "All measures" = "all"),
                      selected = "selected",
                      inline = TRUE
                    )
                  )
                ),
                uiOutput("enz_pval_filter_ui"),
                DTOutput("enz_pvals_table")
              )
            ),
            tabPanel("Table", DTOutput("enz_table"))
          )
        )
      )
    ),

    nav_panel(
      "Objective 4: Volatiles (GC)",
      layout_sidebar(
        sidebar = sidebar(
          width = 360,
          h4("GC options"),
          selectInput(
            "gc_dup_policy",
            "Duplicate compound names within a timepoint",
            choices = c(
              "Sum duplicates" = "sum",
              "Mean duplicates" = "mean",
              "Keep distinct (suffix #1/#2)" = "keep_distinct"
            ),
            selected = "sum"
          ),
          sliderInput("gc_top_n", "Top N compounds per timepoint (by mean RA)", min = 1, max = 100, value = 20, step = 1),
          uiOutput("gc_heatmap_strain_ui"),
          selectInput(
            "gc_heatmap_ncol",
            "Heatmap facet columns (All strains)",
            choices = c("1" = 1, "2 (Recommended)" = 2, "3" = 3),
            selected = 2
          ),
          hr(),
          h4("Downloads"),
          downloadButton("dl_gc_heatmap", "Download heatmap (PNG)"),
          downloadButton("dl_gc_compare", "Download compare plot (PNG)"),
          downloadButton("dl_gc_pca", "Download PCA (PNG)"),
          downloadButton("dl_gc_loadings", "Download PCA loadings (CSV)"),
          downloadButton("dl_gc_table", "Download table (XLSX)")
        ),
        card(
          card_header("GC volatile profiles (exploratory; no replicates)"),
          div(
            class = "app-note",
            tags$p("GC data contains one measurement per strain/timepoint block in the current file. The app provides descriptive/exploratory visualizations only; no p-values/CIs are produced for compound-level comparisons."),
            tags$p(tags$b("Data sources (file / sheets):")),
            tags$ul(
              tags$li(tags$code("Data/GC_VolatilesYeast.xlsx"), " / ", tags$code("QC-All"), " (present; not used in v1)."),
              tags$li(
                tags$code("Data/GC_VolatilesYeast.xlsx"),
                " / strain sheets: ",
                tags$code("SC-Volatiles"), ", ",
                tags$code("Y1-Volatiles"), ", ",
                tags$code("Y1A-Volatiles"), ", ",
                tags$code("Y3A-Volatiles"), ", ",
                tags$code("Y4-Volatiles"), ", ",
                tags$code("Y5-volatiles"), ", ",
                tags$code("Y6B-Volatiles"),
                "."
              ),
              tags$li(tags$code("exampleMapData.xlsx"), " / ", tags$code("compound_map"), " (compound class / sensory mapping coverage).")
            )
          ),
          tabsetPanel(
            type = "tabs",
            tabPanel(
              "Heatmap",
              div(class = "app-note", tags$p("Use the strain selector to focus on a single strain (cleaner) or show All strains (faceted).")),
              plotOutput("gc_heatmap", height = 560)
            ),
            tabPanel(
              "Compare",
              uiOutput("gc_compare_ui"),
              plotOutput("gc_compare_plot", height = 520),
              hr(),
              plotOutput("gc_snapshot_plot", height = 380)
            ),
            tabPanel(
              "PCA",
              uiOutput("gc_pca_ui"),
              plotOutput("gc_pca", height = 520),
              hr(),
              h4("Top loadings (helps interpret PCs)"),
              DTOutput("gc_pca_loadings")
            ),
            tabPanel(
              "Summary",
              plotOutput("gc_total_plot", height = 360),
              hr(),
              plotOutput("gc_richness_plot", height = 360)
            ),
            tabPanel("Mapping coverage", DTOutput("gc_map_coverage")),
            tabPanel("Table", DTOutput("gc_table"))
          ),
          hr(),
          div(class = "app-note", tags$p("All GC outputs depend on the duplicate handling policy and are computed from the tidy internal dataset."))
        )
      )
    ),

		  nav_panel(
		    "Objective 5: Metabolic dynamics",
		    layout_sidebar(
		      sidebar = sidebar(
		        width = 360,
		        h4("Selection"),
            selectInput(
              "met_source",
              "Data source",
              choices = c(
                "Metabolites (Yeast metabolites)" = "met",
                "Autolytic-full" = "aut",
                "Reducing sugar-24h" = "red24"
              ),
              selected = "met"
            ),
		        uiOutput("met_select_ui"),
		        hr(),
		        h4("Technical duplicates"),
		        selectInput(
	          "met_rep_handling",
	          "Replicate handling",
	          choices = c(
	            "Mean(S-1, S-dup) (Recommended)" = "mean",
	            "S-1 only" = "s1",
	            "S-dup only" = "sdup",
	            "Difference (S-1 ? S-dup) (QC)" = "diff"
	          ),
	          selected = "mean"
	        ),
	        checkboxInput("met_show_points", "Show raw points (detected only)", value = FALSE),
		        hr(),
		        h4("Downloads"),
		        downloadButton("dl_met_plot", "Download dynamics (PNG)"),
		        downloadButton("dl_met_snapshot_plot", "Download compare (PNG)"),
		        downloadButton("dl_met_heatmap_plot", "Download heatmap (PNG)"),
            downloadButton("dl_met_pca_plot", "Download PCA (PNG)"),
            downloadButton("dl_met_pca_loadings", "Download PCA loadings (CSV)"),
            downloadButton("dl_met_yield_plot", "Download yield plot (PNG)"),
            downloadButton("dl_met_yield_table", "Download yield table (XLSX)"),
		        downloadButton("dl_met_pvals", "Download p-values (CSV)"),
		        downloadButton("dl_met_table", "Download table (XLSX)")
		      ),
		      card(
		        card_header("Metabolic dynamics (descriptive; mono-culture)"),
		        div(
		          class = "app-note",
			          tags$p("Time series are summarized from technical duplicates (S-1/S-dup). 0 is treated as non-detect (excluded from means; detection rates are reported)."),
			          tags$p(tags$b("Peer-review note:"), " p-values are optional and are based on technical duplicates only (exploratory), not independent fermentations."),
			          tags$p(tags$b("Data sources (file / sheet):")),
			          tags$ul(
			            tags$li(tags$code("Data/Metabolites.xlsx"), " / ", tags$code("Yeast metabolites"), " (metabolites time series)."),
			            tags$li(tags$code("Data/Yeast tolerance test.xlsx"), " / ", tags$code("Autolytic-full"), " (FAN, protein, DNS reducing sugar, survival, OD context)."),
			            tags$li(tags$code("Data/Metabolites.xlsx"), " / ", tags$code("Reducing sugar-24h"), " (24h reducing sugar summaries).")
			          )
			        ),
			        tabsetPanel(
			          type = "tabs",
			          tabPanel("Dynamics", plotOutput("met_plot", height = 520)),
		          tabPanel(
	            "Timepoint compare",
	            fluidRow(
	              column(4, uiOutput("met_timepoint_ui")),
	              column(4, uiOutput("met_delta_ui")),
	              column(4, selectInput("met_compare_geom", "Geometry", choices = c("Dot" = "dot", "Bar" = "bar"), selected = "dot"))
	            ),
	            plotOutput("met_snapshot_plot", height = 420),
	            hr(),
	            plotOutput("met_delta_plot", height = 360)
	          ),
		          tabPanel("Heatmap", plotOutput("met_heatmap_plot", height = 520)),
              tabPanel(
                "Multivariate (PCA)",
                uiOutput("met_pca_ui"),
                plotOutput("met_pca_plot", height = 520),
                hr(),
                h4("Top loadings"),
                DTOutput("met_pca_loadings")
              ),
              tabPanel(
                "Yield/Efficiency",
                div(class = "app-note", tags$p("Uses detected values only (no imputation). If required analytes or timepoints are missing/non-detect, results show as Not available.")),
                uiOutput("met_yield_ui"),
                plotOutput("met_yield_plot", height = 520),
                hr(),
                DTOutput("met_yield_table")
              ),
		          tabPanel(
		            "Stats (p-values)",
		            div(
	              class = "app-note",
	              tags$p(tags$b("Important:"), " p-values in v1 are based on technical duplicates (repeat measurements of the same sample/run).")
	            ),
	            checkboxInput("met_enable_pvals", "Enable p-values (technical duplicates; exploratory)", value = FALSE),
	            conditionalPanel(
	              condition = "input.met_enable_pvals === true",
	              fluidRow(
	                column(
	                  4,
	                  selectInput(
	                    "met_pval_method",
	                    "Method",
	                    choices = c("ANOVA + Tukey" = "aov_tukey", "Kruskal + Wilcoxon" = "kw_wilcox"),
	                    selected = "aov_tukey"
	                  )
	                ),
	                column(
	                  4,
	                  selectInput(
	                    "met_pval_view",
	                    "Compare",
	                    choices = c("Strains at each timepoint" = "by_time", "Timepoints within each strain" = "by_strain"),
	                    selected = "by_time"
	                  )
	                ),
	                column(
	                  4,
	                  selectInput(
	                    "met_padj_overall",
	                    "Overall adjustment (across all pairwise tests shown)",
	                    choices = c("BH (FDR)" = "BH", "Holm" = "holm", "Bonferroni" = "bonferroni"),
	                    selected = "BH"
	                  )
	                )
	              ),
	              uiOutput("met_pval_filter_ui"),
	              DTOutput("met_pvals_table")
	            )
	          ),
	          tabPanel("Summary table", DTOutput("met_table"))
	        )
	      )
	    )
	  ),

	  nav_panel(
	    "Objective 6: Best yeast",
	    layout_sidebar(
	      sidebar = sidebar(
	        width = 360,
	        h4("Ranking weights (transparent)"),
	        selectInput(
	          "rank_preset",
	          "Preset weights",
	          choices = c(
	            "Balanced (Recommended)" = "balanced",
	            "Equal weights" = "equal",
	            "Growth-focused" = "growth",
	            "Tolerance-focused" = "tolerance",
	            "Enzymes-focused" = "enzymes",
	            "Volatiles-focused" = "volatiles",
	            "Custom (use sliders)" = "custom"
	          ),
	          selected = "balanced"
	        ),
	        sliderInput("w_growth", "Growth weight", min = 0, max = 1, value = 0.35, step = 0.05),
	        sliderInput("w_tolerance", "Tolerance weight", min = 0, max = 1, value = 0.25, step = 0.05),
	        sliderInput("w_enz", "Enzymes weight", min = 0, max = 1, value = 0.20, step = 0.05),
	        sliderInput("w_vol", "Volatiles weight", min = 0, max = 1, value = 0.20, step = 0.05),
	        div(class = "app-note", tags$p("Weights are normalized internally to sum to 1. Ranking is exploratory given the current design (no biological replication).")),
	        hr(),
	        h4("Downloads"),
	        downloadButton("dl_rank_plot", "Download ranking plot (PNG)"),
	        downloadButton("dl_rank_contrib_plot", "Download contributions (PNG)"),
	        downloadButton("dl_rank_sensitivity_plot", "Download sensitivity (PNG)"),
	        downloadButton("dl_rank_table", "Download results (XLSX)")
	      ),
	      card(
	        card_header("Multi-criteria ranking (exploratory)"),
	        div(
	          class = "app-note",
	          tags$p("This tab computes a reproducible, transparent score using only available endpoints. Co-culture synergy/antagonism is explicitly excluded in v1. The app reports which metrics were used and which were not available."),
	          tags$p(tags$b("Data sources (file / sheets):")),
	          tags$ul(
	            tags$li(tags$code("Data/Metabolites.xlsx"), " / ", tags$code("Yeast metabolites"), " (growth metrics)."),
	            tags$li(tags$code("Data/Yeast tolerance test.xlsx"), " / ", tags$code("Autolytic-full"), " (survival metrics)."),
	            tags$li(tags$code("Data/Yeast tolerance test.xlsx"), " / ", tags$code("tolerance test-broth"), " and ", tags$code("tolerance test-agar"), " (tolerance metrics)."),
	            tags$li(tags$code("Data/Yeast tolerance test.xlsx"), " / ", tags$code("Enzyme production"), " (enzyme production rate metrics)."),
	            tags$li(tags$code("Data/GC_VolatilesYeast.xlsx"), " / strain sheets (volatiles metrics; depends on GC duplicate policy).")
	          )
	        ),
	        tabsetPanel(
	          type = "tabs",
	          tabPanel(
	            "Ranking",
	            fluidRow(
	              column(6, plotOutput("rank_plot", height = 420)),
	              column(6, DTOutput("rank_table"))
	            ),
	            hr(),
	            h4("Weights used (after normalization)"),
	            DTOutput("rank_weights")
	          ),
	          tabPanel(
	            "Components",
	            plotOutput("rank_contrib_plot", height = 420),
	            hr(),
	            DTOutput("rank_components_table")
	          ),
	          tabPanel(
	            "Sensitivity",
	            div(
	              class = "app-note",
	              tags$p("Sensitivity analysis shows how ranks change under pre-defined weight presets. This does not add biological replication; it only checks ranking robustness to weight choices.")
	            ),
	            plotOutput("rank_sensitivity_plot", height = 420),
	            hr(),
	            DTOutput("rank_sensitivity_table")
	          ),
	          tabPanel(
	            "Trade-offs",
	            uiOutput("rank_trade_ui"),
	            plotOutput("rank_trade_plot", height = 520)
	          )
	        )
	      )
	    )
	  ),

  nav_panel(
    "QC & Downloads",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        h4("Technical duplicate QC"),
        selectInput(
          "qc_dataset",
          "Dataset",
          choices = c(
            "Metabolites (Yeast metabolites)" = "met",
            "Autolytic" = "aut",
            "Tolerance broth" = "tol",
            "Enzymes" = "enz"
          )
        ),
        uiOutput("qc_analyte_ui"),
        hr(),
        h4("Downloads"),
        downloadButton("dl_qc_plot", "Download QC plot (PNG)")
      ),
      card(
        card_header("S-1 vs S-dup agreement (QC)"),
        div(
          class = "app-note",
          tags$p("These are technical duplicates (repeat measurements of the same sample/run). Large disagreement flags potential measurement issues; it does not represent biological variability.")
        ),
        plotOutput("qc_plot", height = 520)
      )
    )
  )
)

server <- function(input, output, session) {
  data_all <- reactive({
    req(file.exists(paths$metabolites), file.exists(paths$tolerance))
    list(
      metabolites = load_metabolites(paths$metabolites),
      reducing = load_reducing_sugar(paths$metabolites),
      tol_agar = load_tolerance_agar(paths$tolerance),
      tol_broth = load_tolerance_broth(paths$tolerance),
      autolytic = load_autolytic(paths$tolerance),
      enzymes = load_enzymes(paths$tolerance),
      gc = if (file.exists(paths$gc)) load_gc_volatiles(paths$gc) else tibble::tibble(),
      compound_map = load_compound_map(paths$map)
    )
  })

  output$data_status <- renderUI({
    ok <- c(
      Metabolites = file.exists(paths$metabolites),
      Tolerance = file.exists(paths$tolerance),
      `GC volatiles` = file.exists(paths$gc),
      `Mapping (example)` = file.exists(paths$map)
    )
    tags$div(
      tags$ul(
        lapply(names(ok), function(nm) {
          tags$li(tags$b(nm), ": ", if (ok[[nm]]) "found" else "missing")
        })
      ),
      tags$p(class = "app-note", tags$b("Project root:"), " ", project_root)
    )
  })

  coverage <- reactive({
    objective_coverage(data_all())
  })

    output$tbl_coverage <- renderDT({
      datatable(
        coverage()$summary,
        rownames = FALSE,
        # Small table; disable DT controls to avoid layout overlap on the Overview page.
        options = list(dom = "t", paging = FALSE, searching = FALSE, info = FALSE, ordering = FALSE)
      )
    })

    output$tbl_tol_items <- renderDT({
      datatable(
        coverage()$tol_items,
        rownames = FALSE,
        options = list(dom = "t", paging = FALSE, searching = FALSE, info = FALSE, ordering = FALSE)
      )
    })

    output$tbl_enz_items <- renderDT({
      datatable(
        coverage()$enz_items,
        rownames = FALSE,
        options = list(dom = "t", paging = FALSE, searching = FALSE, info = FALSE, ordering = FALSE)
      )
    })

  output$dl_session <- downloadHandler(
    filename = function() paste0("sessionInfo_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content = function(file) {
      con <- file(file, open = "wt", encoding = "UTF-8")
      on.exit(close(con), add = TRUE)
      writeLines(capture.output(sessionInfo()), con = con)
    }
  )

  output$dl_params <- downloadHandler(
    filename = function() paste0("analysis_params_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"),
    content = function(file) {
      params <- list(
        timestamp = as.character(Sys.time()),
        project_root = project_root,
        paths = paths,
        gc_duplicate_policy = input$gc_dup_policy %||% NA_character_,
        non_detect_policy = "0 means not detected; treated as NA for means/models; no imputation; report detection frequency",
        inference_policy = "No p-values/CIs as evidence of between-strain differences in v1"
      )
      writeLines(toJSON(params, pretty = TRUE, auto_unbox = TRUE), con = file)
    }
  )

  ## Objective 1: Growth
  growth_data <- reactive({
    d <- data_all()
    if (input$growth_endpoint == "met_growth") {
      df <- d$metabolites |>
        filter(standardize_key(type_unit) == standardize_key("microbial growth Log (CFU/mL)"))
      list(df = df, ylab = "Microbial growth (log CFU/mL)", title = "Microbial growth (Metabolites.xlsx)")
    } else if (input$growth_endpoint == "aut_survival") {
      df <- d$autolytic |>
        filter(analyte_key == standardize_key("microbial-survival- Log (CFU/mL)"))
      list(df = df, ylab = "Survival (log CFU/mL)", title = "Microbial survival (Autolytic)")
    } else {
      df <- d$autolytic |>
        filter(analyte_key == standardize_key("OD"))
      list(df = df, ylab = "OD (AU)", title = "Optical density (Autolytic)")
    }
  })

    growth_df_use <- reactive({
      gd <- growth_data()$df
      if (nrow(gd) == 0) return(tibble::tibble())
      apply_techrep_handling(
        gd,
        group_cols = c("strain", "time_h", "analyte"),
        handling = input$growth_rep_handling %||% "mean",
        value_col = "value"
      )
    })

  growth_timepoints <- reactive({
    gd <- growth_data()$df
    sort(unique(gd$time_h))
  })

  output$growth_timepoint_ui <- renderUI({
    tp <- growth_timepoints()
    selectInput("growth_snapshot_time", "Snapshot time (h)", choices = tp, selected = tp[[1]] %||% NULL)
  })

  output$growth_delta_ui <- renderUI({
    tp <- growth_timepoints()
    if (length(tp) < 2) return(NULL)
    tagList(
      selectInput("growth_delta_t1", "t1 (h)", choices = tp, selected = tp[[1]] %||% NULL),
      selectInput("growth_delta_t2", "t2 (h)", choices = tp, selected = tp[[length(tp)]] %||% NULL)
    )
  })

  growth_table_df <- reactive({
    gd <- growth_data()$df
    if (nrow(gd) == 0) return(tibble::tibble())

      wide <- gd |>
        dplyr::select(strain, time_h, analyte, tech_rep, value = .data$value) |>
        dplyr::group_by(.data$strain, .data$time_h, .data$analyte, .data$tech_rep) |>
        dplyr::summarise(
          value = if (all(is.na(.data$value))) NA_real_ else mean(.data$value, na.rm = TRUE),
          .groups = "drop"
        ) |>
        tidyr::pivot_wider(names_from = tech_rep, values_from = value)

      wide |>
        dplyr::rowwise() |>
        dplyr::mutate(
          mean_value = mean(dplyr::c_across(dplyr::any_of(c("S-1", "S-dup"))), na.rm = TRUE),
          sd_value = {
            x <- dplyr::c_across(dplyr::any_of(c("S-1", "S-dup")))
            if (sum(!is.na(x)) < 2) NA_real_ else stats::sd(x, na.rm = TRUE)
          },
          diff_s1_sdup = .data$`S-1` - .data$`S-dup`,
          cv_pct = dplyr::if_else(!is.na(.data$mean_value) & .data$mean_value != 0, 100 * .data$sd_value / .data$mean_value, NA_real_)
        ) |>
        dplyr::ungroup() |>
        dplyr::arrange(.data$strain, .data$time_h)
    })

    growth_plot_obj <- reactive({
      g <- growth_data()
      gd <- g$df
      if (nrow(gd) == 0) return(ggplot() + labs(title = "Not available in current dataset"))

      df_line <- growth_df_use()
      pal <- strain_palette(df_line$strain)
      rep_handling <- input$growth_rep_handling %||% "mean"
      ylab <- paste0(g$ylab, " (", rep_handling, ")")

    p_main <- df_line |>
      ggplot(aes(x = time_h, y = value_use, color = strain, group = strain)) +
      geom_line(linewidth = 0.8, alpha = 0.9) +
      geom_point(size = 2) +
      scale_color_manual(values = pal) +
      labs(title = g$title, x = "Time (h)", y = ylab, color = "Strain") +
      khadijat_theme()

      if (!isTRUE(input$growth_show_points)) return(p_main)

      # Show raw technical duplicates as faint points when plotting the mean.
      if (rep_handling != "mean") return(p_main)

    p_main +
      geom_point(
        data = gd,
        aes(x = time_h, y = value, color = strain, shape = tech_rep),
        size = 1.6,
        alpha = 0.55,
        inherit.aes = FALSE
      ) +
      scale_shape_manual(values = c("S-1" = 16, "S-dup" = 1)) +
      guides(shape = guide_legend(title = "Tech rep"))
  })

  output$growth_plot <- renderPlot({
    growth_plot_obj()
  })

    growth_bar_plot_obj <- reactive({
      g <- growth_data()
      df <- growth_df_use()
      if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))
      pal <- strain_palette(df$strain)
      rep_handling <- input$growth_rep_handling %||% "mean"

      pos <- if (identical(input$growth_bar_position, "stack")) "stack" else "dodge"

    ggplot(df, aes(x = factor(time_h), y = value_use, fill = strain)) +
      geom_col(position = pos, width = 0.85, color = "white") +
      scale_fill_manual(values = pal) +
        labs(
          title = paste0(g$title, " (bar chart)"),
          subtitle = paste0("Replicate handling: ", rep_handling),
          x = "Time (h)",
          y = g$ylab,
          fill = "Strain"
        ) +
	      khadijat_theme() +
	      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"))
  })

  output$growth_bar_plot <- renderPlot({ growth_bar_plot_obj() })

    growth_snapshot_plot_obj <- reactive({
      g <- growth_data()
      df <- growth_df_use()
      req(nrow(df) > 0)
      req(input$growth_snapshot_time)
      t0 <- as.integer(input$growth_snapshot_time)
      df_t <- df |>
        dplyr::filter(.data$time_h == t0)

      pal <- strain_palette(df_t$strain)
      rep_handling <- input$growth_rep_handling %||% "mean"

      if (identical(input$growth_compare_geom, "bar")) {
      ggplot(df_t, aes(x = reorder(strain, value_use), y = value_use, fill = strain)) +
        geom_col(width = 0.8, color = "white") +
        coord_flip() +
        scale_fill_manual(values = pal, guide = "none") +
          labs(
            title = paste0(g$title, " (snapshot at ", t0, " h)"),
            subtitle = paste0("Replicate handling: ", rep_handling),
            x = "Strain",
            y = g$ylab
          ) +
          khadijat_theme()
      } else {
      ggplot(df_t, aes(x = reorder(strain, value_use), y = value_use, color = strain)) +
        geom_point(size = 3, alpha = 0.9) +
        coord_flip() +
        scale_color_manual(values = pal, guide = "none") +
          labs(
            title = paste0(g$title, " (snapshot at ", t0, " h)"),
            subtitle = paste0("Replicate handling: ", rep_handling),
            x = "Strain",
            y = g$ylab
          ) +
          khadijat_theme()
      }
  })

  output$growth_snapshot_plot <- renderPlot({ growth_snapshot_plot_obj() })

    growth_delta_plot_obj <- reactive({
      g <- growth_data()
      df <- growth_df_use()
      req(nrow(df) > 0)
      req(input$growth_delta_t1, input$growth_delta_t2)
      t1 <- as.integer(input$growth_delta_t1)
      t2 <- as.integer(input$growth_delta_t2)
      if (is.na(t1) || is.na(t2) || t1 == t2) {
        return(ggplot() + labs(title = "Select two different timepoints"))
      }

    w <- df |>
      dplyr::select(strain, time_h, value_use) |>
      tidyr::pivot_wider(names_from = time_h, values_from = value_use)

    if (!(as.character(t1) %in% names(w)) || !(as.character(t2) %in% names(w))) {
      return(ggplot() + labs(title = "Selected timepoints not available"))
    }

      w <- w |>
        dplyr::mutate(delta = .data[[as.character(t2)]] - .data[[as.character(t1)]])

      pal <- strain_palette(w$strain)
      rep_handling <- input$growth_rep_handling %||% "mean"
      ggplot(w, aes(x = reorder(strain, delta), y = delta, fill = strain)) +
        geom_col(width = 0.8, color = "white") +
        coord_flip() +
        scale_fill_manual(values = pal, guide = "none") +
        labs(
          title = paste0(g$title, " (delta: ", t2, " h − ", t1, " h)"),
          subtitle = paste0("Replicate handling: ", rep_handling),
          x = "Strain",
          y = paste0("Delta (", g$ylab, ")")
        ) +
        khadijat_theme()
    })

  output$growth_delta_plot <- renderPlot({ growth_delta_plot_obj() })

    growth_heatmap_plot_obj <- reactive({
      g <- growth_data()
      df <- growth_df_use()
      if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))
      rep_handling <- input$growth_rep_handling %||% "mean"
      ggplot(df, aes(x = factor(time_h), y = strain, fill = value_use)) +
        geom_tile(color = "white", linewidth = 0.2) +
        scale_fill_viridis_c(option = "C", na.value = "#f1f1f1") +
        labs(
          title = paste0(g$title, " (heatmap)"),
          subtitle = paste0("Replicate handling: ", rep_handling),
          x = "Time (h)",
          y = "Strain",
          fill = g$ylab
        ) +
	        khadijat_theme() +
	        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"))
	    })

  output$growth_heatmap_plot <- renderPlot({ growth_heatmap_plot_obj() })

    output$growth_table <- renderDT({
      df <- growth_table_df()
      coldefs <- dt_coldefs_numeric(df, digits = 4L)
      datatable(
        df,
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          columnDefs = coldefs
        )
      )
    })

  output$dl_growth_plot <- downloadHandler(
    filename = function() paste0("growth_", input$growth_endpoint, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(growth_plot_obj(), file, width = 10, height = 5.5, dpi = 320)
  )

  output$dl_growth_snapshot_plot <- downloadHandler(
    filename = function() paste0("growth_snapshot_", input$growth_endpoint, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(growth_snapshot_plot_obj(), file, width = 10, height = 5.5, dpi = 320)
  )

  growth_model_fit_all <- reactive({
    # Model fitting is not meaningful for the QC "difference" view.
    rep_handling <- input$growth_rep_handling %||% "mean"
    if (identical(rep_handling, "diff")) {
      return(list(fits = tibble::tibble(), pred = tibble::tibble(), resid = tibble::tibble()))
    }

    df <- growth_df_use()
    if (nrow(df) == 0) return(list(fits = tibble::tibble(), pred = tibble::tibble(), resid = tibble::tibble()))

    models <- input$growth_model_models %||% c("logistic", "gompertz", "baranyi")
    fit_growth_models(df |>
      dplyr::transmute(strain = .data$strain, time_h = .data$time_h, y = .data$value_use), models = models)
  })

  growth_model_fit_use <- reactive({
    all <- growth_model_fit_all()
    fits <- all$fits
    if (nrow(fits) == 0) return(all)
    disp <- input$growth_model_display %||% "best"
    if (identical(disp, "all")) return(all)

    best <- fits |>
      mutate(
        aic_key = dplyr::if_else(is.na(.data$aic), Inf, .data$aic),
        rmse_key = dplyr::if_else(is.na(.data$rmse), Inf, .data$rmse)
      ) |>
      group_by(.data$strain) |>
      arrange(.data$aic_key, .data$rmse_key) |>
      slice(1) |>
      ungroup() |>
      select(strain, model)

    list(
      fits = fits |> semi_join(best, by = c("strain", "model")),
      pred = all$pred |> semi_join(best, by = c("strain", "model")),
      resid = all$resid |> semi_join(best, by = c("strain", "model"))
    )
  })

  growth_model_plot_obj <- reactive({
    g <- growth_data()
    df_obs <- growth_df_use()
    rep_handling <- input$growth_rep_handling %||% "mean"
    if (identical(rep_handling, "diff")) {
      return(ggplot() + labs(title = "Model fit not available", subtitle = "Disable Difference (QC) to fit growth models."))
    }

    fit_use <- growth_model_fit_use()
    pred <- fit_use$pred
    if (nrow(df_obs) == 0 || nrow(pred) == 0) return(ggplot() + labs(title = "Model fit not available"))

    ncol_f <- as.integer(input$growth_model_ncol %||% 2L)

    df_obs2 <- df_obs |>
      dplyr::transmute(strain = .data$strain, time_h = .data$time_h, y = .data$value_use) |>
      dplyr::filter(is.finite(.data$time_h) & is.finite(.data$y))

    ggplot() +
      geom_point(data = df_obs2, aes(x = time_h, y = y), size = 2, alpha = 0.85, color = "#1d2733") +
      geom_line(data = pred, aes(x = time_h, y = y_pred, color = model, group = model), linewidth = 0.9, alpha = 0.9) +
      facet_wrap(~strain, ncol = ncol_f, scales = "free_y") +
      scale_color_viridis_d(option = "C") +
      scale_x_continuous(breaks = sort(unique(df_obs2$time_h))) +
      labs(
        title = paste0(g$title, " (model fits; descriptive)"),
        subtitle = paste0(
          "Replicate handling: ", rep_handling,
          " | Display: ", (input$growth_model_display %||% "best")
        ),
        x = "Time (h)",
        y = g$ylab,
        color = "Model"
      ) +
      khadijat_theme()
  })

  output$growth_model_plot <- renderPlot({ growth_model_plot_obj() })

  growth_model_resid_plot_obj <- reactive({
    rep_handling <- input$growth_rep_handling %||% "mean"
    if (identical(rep_handling, "diff")) return(ggplot() + labs(title = "Not available for Difference (QC)"))

    fit_use <- growth_model_fit_use()
    r <- fit_use$resid
    if (nrow(r) == 0) return(ggplot() + labs(title = "Not available in current dataset"))

    ncol_f <- as.integer(input$growth_model_ncol %||% 2L)

    ggplot(r, aes(x = y_fit, y = resid, color = model)) +
      geom_hline(yintercept = 0, linewidth = 0.6, color = "#333333") +
      geom_point(size = 2, alpha = 0.8) +
      facet_wrap(~strain, ncol = ncol_f, scales = "free") +
      scale_color_viridis_d(option = "C") +
      labs(
        title = "Residual diagnostics (residual vs fitted)",
        x = "Fitted value",
        y = "Residual (obs - fitted)",
        color = "Model"
      ) +
      khadijat_theme()
  })

  output$growth_model_resid_plot <- renderPlot({ growth_model_resid_plot_obj() })

  growth_model_table_df <- reactive({
    rep_handling <- input$growth_rep_handling %||% "mean"
    if (identical(rep_handling, "diff")) return(tibble::tibble())
    f <- growth_model_fit_use()$fits
    if (nrow(f) == 0) return(f)
    f |>
      transmute(
        strain,
        model,
        n_points,
        aic,
        rmse,
        max_pred,
        t_inflect,
        mu_max,
        lag_tangent,
        t_95 = .data$t_frac,
        Asym,
        xmid,
        scal,
        b2,
        b3,
        y0,
        ymax,
        mu,
        lag
      ) |>
      arrange(.data$strain, .data$aic)
  })

  output$growth_model_table <- renderDT({
    df <- growth_model_table_df()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  output$dl_growth_model_plot <- downloadHandler(
    filename = function() paste0("growth_model_fit_", input$growth_endpoint, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(growth_model_plot_obj(), file, width = 12, height = 7.5, dpi = 320)
  )

  output$dl_growth_model_table <- downloadHandler(
    filename = function() paste0("growth_model_fit_", input$growth_endpoint, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) {
      all <- growth_model_fit_all()
      use <- growth_model_fit_use()
      writexl::write_xlsx(
        list(
          fits_all = all$fits,
          fits_displayed = use$fits,
          predictions_displayed = use$pred,
          residuals_displayed = use$resid
        ),
        path = file
      )
    }
  )

  growth_pvals_by_time <- reactive({
    req(isTRUE(input$growth_enable_pvals))
    gd <- growth_data()$df
    if (nrow(gd) == 0) return(tibble::tibble())

    method <- input$growth_pval_method %||% "aov_tukey"
    out <- list()
    k <- 0L
    for (t in sort(unique(gd$time_h))) {
      dat <- gd |>
        dplyr::filter(.data$time_h == t) |>
        dplyr::filter(!is.na(.data$value))
      if (nrow(dat) == 0) next

      res <- if (identical(method, "kw_wilcox")) safe_kw_wilcox(dat, response_col = "value", group_col = "strain") else safe_aov_tukey(dat, response_col = "value", group_col = "strain")

      if (nrow(res$pairwise) == 0) next
      k <- k + 1L
      out[[k]] <- res$pairwise |>
        dplyr::mutate(
          context = "by_time",
          time_h = t,
          strain = NA_character_,
          global_p = res$global_p,
          p_adj_within = .data$p_adj,
          method = method
        )
    }
    dplyr::bind_rows(out)
  })

  growth_pvals_by_strain <- reactive({
    req(isTRUE(input$growth_enable_pvals))
    gd <- growth_data()$df
    if (nrow(gd) == 0) return(tibble::tibble())

    method <- input$growth_pval_method %||% "aov_tukey"
    out <- list()
    k <- 0L
    for (s in sort(unique(gd$strain))) {
      dat <- gd |>
        dplyr::filter(.data$strain == s) |>
        dplyr::filter(!is.na(.data$value))
      if (nrow(dat) == 0) next

      # Compare timepoints within strain.
      res <- if (identical(method, "kw_wilcox")) safe_kw_wilcox(dat, response_col = "value", group_col = "time_h") else safe_aov_tukey(dat, response_col = "value", group_col = "time_h")

      if (nrow(res$pairwise) == 0) next
      k <- k + 1L
      out[[k]] <- res$pairwise |>
        dplyr::mutate(
          context = "by_strain",
          time_h = NA_integer_,
          strain = s,
          global_p = res$global_p,
          p_adj_within = .data$p_adj,
          method = method
        )
    }
    dplyr::bind_rows(out)
  })

  output$growth_pval_filter_ui <- renderUI({
    gd <- growth_data()$df
    if (nrow(gd) == 0) return(NULL)
    view <- input$growth_pval_view %||% "by_time"
    if (identical(view, "by_time")) {
      tp <- sort(unique(gd$time_h))
      selectInput("growth_pval_time", "Filter time (h)", choices = c("All" = "__all__", tp), selected = "__all__")
    } else {
      strains <- sort(unique(gd$strain))
      selectInput("growth_pval_strain", "Filter strain", choices = c("All" = "__all__", strains), selected = "__all__")
    }
  })

  growth_pvals_filtered <- reactive({
    req(isTRUE(input$growth_enable_pvals))
    view <- input$growth_pval_view %||% "by_time"
    padj_overall <- input$growth_padj_overall %||% "BH"

    df <- if (identical(view, "by_strain")) growth_pvals_by_strain() else growth_pvals_by_time()
    if (nrow(df) == 0) return(df)

    if (identical(view, "by_time") && !is.null(input$growth_pval_time) && input$growth_pval_time != "__all__") {
      df <- df |>
        dplyr::filter(.data$time_h == as.integer(input$growth_pval_time))
    }
    if (identical(view, "by_strain") && !is.null(input$growth_pval_strain) && input$growth_pval_strain != "__all__") {
      df <- df |>
        dplyr::filter(.data$strain == input$growth_pval_strain)
    }

    df |>
      dplyr::mutate(
        p_adj_overall = stats::p.adjust(.data$p_adj_within, method = padj_overall)
      ) |>
      dplyr::arrange(.data$p_adj_overall)
  })

      output$growth_pvals_table <- renderDT({
        req(isTRUE(input$growth_enable_pvals))
        df <- growth_pvals_filtered()
        coldefs <- dt_coldefs_numeric(df, digits = 4L)
        datatable(
          df,
          rownames = FALSE,
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            columnDefs = coldefs
          )
        )
      })

  output$dl_growth_pvals <- downloadHandler(
    filename = function() paste0("growth_pvalues_", input$growth_endpoint, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      readr::write_csv(growth_pvals_filtered(), file, na = "")
    }
  )

  output$dl_growth_table <- downloadHandler(
    filename = function() paste0("growth_summary_", input$growth_endpoint, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) writexl::write_xlsx(list(summary = growth_table_df()), path = file)
  )

  ## Objective 2: Tolerance
  output$tol_filter_ui <- renderUI({
    d <- data_all()
    if (input$tol_view == "agar") {
      choices <- sort(unique(d$tol_agar$analyte))
    } else {
      choices <- sort(unique(d$tol_broth$analyte))
    }
    selectInput("tol_analyte", "Analyte group", choices = choices, selected = choices[[1]] %||% NULL)
  })

  output$tol_rep_ui <- renderUI({
    if (!identical(input$tol_view, "broth")) return(NULL)
    tagList(
      hr(),
      h4("Technical duplicates"),
      selectInput(
        "tol_rep_handling",
        "Replicate handling (broth)",
        choices = c(
          "Mean(S-1, S-dup) (Recommended)" = "mean",
          "S-1 only" = "s1",
          "S-dup only" = "sdup",
          "Difference (S-1 − S-dup) (QC)" = "diff"
        ),
        selected = "mean"
      )
    )
  })

  tol_table_df <- reactive({
    d <- data_all()
    req(input$tol_analyte)
    if (input$tol_view == "agar") {
      df <- d$tol_agar |>
        filter(analyte == input$tol_analyte) |>
        mutate(value = agar_score)
      df |>
        select(strain, analyte, concentration, value)
    } else {
      df <- d$tol_broth |>
        filter(analyte == input$tol_analyte)

        # Wide summary (keeps zeros as real 0 for broth).
        wide <- df |>
          dplyr::select(strain, analyte, concentration, tech_rep, value = .data$value) |>
          dplyr::group_by(.data$strain, .data$analyte, .data$concentration, .data$tech_rep) |>
          dplyr::summarise(
            value = if (all(is.na(.data$value))) NA_real_ else mean(.data$value, na.rm = TRUE),
            .groups = "drop"
          ) |>
          tidyr::pivot_wider(names_from = tech_rep, values_from = value)

        wide |>
          dplyr::rowwise() |>
          dplyr::mutate(
            mean_value = mean(dplyr::c_across(dplyr::any_of(c("S-1", "S-dup"))), na.rm = TRUE),
            sd_value = {
              x <- dplyr::c_across(dplyr::any_of(c("S-1", "S-dup")))
              if (sum(!is.na(x)) < 2) NA_real_ else stats::sd(x, na.rm = TRUE)
            },
            diff_s1_sdup = .data$`S-1` - .data$`S-dup`,
            cv_pct = dplyr::if_else(!is.na(.data$mean_value) & .data$mean_value != 0, 100 * .data$sd_value / .data$mean_value, NA_real_),
            n_zero = sum(dplyr::c_across(dplyr::any_of(c("S-1", "S-dup"))) == 0, na.rm = TRUE)
          ) |>
          dplyr::ungroup() |>
          dplyr::arrange(.data$strain, .data$concentration)
      }
    })

  tol_heatmap_df <- reactive({
    d <- data_all()
    req(input$tol_analyte)
    if (input$tol_view == "agar") {
      return(d$tol_agar |>
        filter(analyte == input$tol_analyte) |>
        transmute(strain, analyte, concentration, value = agar_score))
    }

    df <- d$tol_broth |>
      filter(analyte == input$tol_analyte)
    rep_handling <- input$tol_rep_handling %||% "mean"
    df_use <- apply_techrep_handling(df, group_cols = c("strain", "analyte", "concentration"), handling = rep_handling, value_col = "value")
    df_use |>
      dplyr::transmute(strain, analyte, concentration, value = value_use)
  })

  tol_plot_obj <- reactive({
    df <- tol_heatmap_df()
    req(nrow(df) > 0)
    ggplot(df, aes(x = concentration, y = strain, fill = value)) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_y_discrete(limits = rev(sort(unique(df$strain)))) +
      scale_fill_viridis_c(option = "C", na.value = "#f1f1f1") +
      labs(
        title = paste0("Tolerance: ", input$tol_analyte, " (", input$tol_view, ")"),
        x = "Condition (concentration)",
        y = "Strain",
        fill = if (input$tol_view == "agar") "Agar score" else "Response"
      ) +
	      khadijat_theme() +
	      theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "bold"))
  })

  output$tol_plot <- renderPlot({
    tol_plot_obj()
  })

    output$tol_table <- renderDT({
      df <- tol_table_df()
      coldefs <- dt_coldefs_numeric(df, digits = 4L)
      datatable(
        df,
        rownames = FALSE,
        options = list(
          pageLength = 12,
          scrollX = TRUE,
          columnDefs = coldefs
        )
      )
    })

	  tol_dose_df <- reactive({
	    d <- data_all()
	    req(input$tol_analyte)
	    req(input$tol_view == "broth")
	    df <- d$tol_broth |>
	      filter(analyte == input$tol_analyte)

    rep_handling <- input$tol_rep_handling %||% "mean"
    df_use <- apply_techrep_handling(df, group_cols = c("strain", "analyte", "concentration"), handling = rep_handling, value_col = "value")

	    # Build a unique concentration -> numeric mapping to avoid many-to-many joins
	    # (each concentration is repeated across strains).
	    cond <- parse_condition_numeric(unique(df_use$concentration))
	    df_use |>
	      dplyr::left_join(cond, by = c("concentration" = "condition_raw")) |>
	      dplyr::mutate(
	        concentration_ord = dplyr::if_else(!is.na(.data$condition_value), sprintf("%g%s", .data$condition_value, .data$condition_unit), .data$concentration)
	      )
	  })

  tol_dose_plot_obj <- reactive({
    req(input$tol_view == "broth")
    df <- tol_dose_df()
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))
    pal <- strain_palette(df$strain)

    x_is_num <- !all(is.na(df$condition_value))
    if (identical(input$tol_dose_geom, "bar")) {
      pos <- if (identical(input$tol_bar_position, "stack")) "stack" else "dodge"
      p <- ggplot(df, aes(x = if (x_is_num) condition_value else factor(concentration, levels = unique(concentration)), y = value_use, fill = strain)) +
        geom_col(position = pos, width = 0.85, color = "white") +
        scale_fill_manual(values = pal) +
        labs(
          title = paste0("Dose-response: ", input$tol_analyte, " (broth)"),
          subtitle = paste0("Replicate handling: ", input$tol_rep_handling %||% "mean"),
          x = "Concentration / level",
          y = "Response",
          fill = "Strain"
        ) +
	        khadijat_theme() +
	        theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "bold"))
      if (x_is_num) p <- p + scale_x_continuous(breaks = sort(unique(df$condition_value)))
      p
    } else {
      p <- ggplot(df, aes(x = if (x_is_num) condition_value else factor(concentration, levels = unique(concentration)), y = value_use, color = strain, group = strain)) +
        geom_line(linewidth = 0.8, alpha = 0.9) +
        geom_point(size = 2) +
        scale_color_manual(values = pal) +
        labs(
          title = paste0("Dose-response: ", input$tol_analyte, " (broth)"),
          subtitle = paste0("Replicate handling: ", input$tol_rep_handling %||% "mean"),
          x = "Concentration / level",
          y = "Response",
          color = "Strain"
        ) +
	        khadijat_theme() +
	        theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "bold"))
      if (x_is_num) p <- p + scale_x_continuous(breaks = sort(unique(df$condition_value)))
      p
    }
  })

  output$tol_dose_plot <- renderPlot({ tol_dose_plot_obj() })

  tol_pvals_all <- reactive({
    req(isTRUE(input$tol_enable_pvals))
    req(identical(input$tol_view, "broth"))
    d <- data_all()
    req(input$tol_analyte)

    df <- d$tol_broth |>
      filter(analyte == input$tol_analyte) |>
      dplyr::filter(!is.na(.data$value))

    method <- input$tol_pval_method %||% "aov_tukey"
    out <- list()
    k <- 0L
    for (c0 in sort(unique(df$concentration))) {
      dat <- df |>
        filter(concentration == c0)
      if (nrow(dat) == 0) next
      res <- if (identical(method, "kw_wilcox")) safe_kw_wilcox(dat, response_col = "value", group_col = "strain") else safe_aov_tukey(dat, response_col = "value", group_col = "strain")
      if (nrow(res$pairwise) == 0) next
      k <- k + 1L
      out[[k]] <- res$pairwise |>
        dplyr::mutate(
          analyte = input$tol_analyte,
          concentration = c0,
          global_p = res$global_p,
          p_adj_within = .data$p_adj,
          method = method
        )
    }
    dplyr::bind_rows(out)
  })

  output$tol_pval_conc_ui <- renderUI({
    d <- data_all()
    if (!identical(input$tol_view, "broth")) return(NULL)
    req(input$tol_analyte)
    conc <- d$tol_broth |>
      filter(analyte == input$tol_analyte) |>
      distinct(concentration) |>
      arrange(concentration) |>
      pull(concentration)
    selectInput("tol_pval_concentration", "Filter concentration", choices = c("All" = "__all__", conc), selected = "__all__")
  })

  tol_pvals_filtered <- reactive({
    req(isTRUE(input$tol_enable_pvals))
    req(identical(input$tol_view, "broth"))
    padj_overall <- input$tol_padj_overall %||% "BH"
    df <- tol_pvals_all()
    if (nrow(df) == 0) return(df)
    if (!is.null(input$tol_pval_concentration) && input$tol_pval_concentration != "__all__") {
      df <- df |>
        dplyr::filter(.data$concentration == input$tol_pval_concentration)
    }
    df |>
      dplyr::mutate(p_adj_overall = stats::p.adjust(.data$p_adj_within, method = padj_overall)) |>
      dplyr::arrange(.data$p_adj_overall)
  })

      output$tol_pvals_table <- renderDT({
        req(isTRUE(input$tol_enable_pvals))
        req(identical(input$tol_view, "broth"))
        df <- tol_pvals_filtered()
        coldefs <- dt_coldefs_numeric(df, digits = 4L)
        datatable(
          df,
          rownames = FALSE,
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            columnDefs = coldefs
          )
        )
      })

  output$dl_tol_plot <- downloadHandler(
    filename = function() paste0("tolerance_", input$tol_view, "_", standardize_key(input$tol_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(tol_plot_obj(), file, width = 12, height = 6.5, dpi = 320)
  )

  output$dl_tol_dose_plot <- downloadHandler(
    filename = function() paste0("tolerance_dose_", standardize_key(input$tol_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(tol_dose_plot_obj(), file, width = 12, height = 6.5, dpi = 320)
  )

  output$dl_tol_pvals <- downloadHandler(
    filename = function() paste0("tolerance_broth_pvalues_", standardize_key(input$tol_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      readr::write_csv(tol_pvals_filtered(), file, na = "")
    }
  )

  output$dl_tol_table <- downloadHandler(
    filename = function() paste0("tolerance_", input$tol_view, "_", standardize_key(input$tol_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) writexl::write_xlsx(list(table = tol_table_df()), path = file)
  )

  ## Objective 2: Composite tolerance scoring (broth; exploratory)
  tol_index_presets <- list(
    balanced = c(mean = 0.40, min = 0.30, auc = 0.30, ic50 = 0.00),
    equal = c(mean = 0.25, min = 0.25, auc = 0.25, ic50 = 0.25),
    mean = c(mean = 0.70, min = 0.15, auc = 0.15, ic50 = 0.00),
    min = c(mean = 0.15, min = 0.70, auc = 0.15, ic50 = 0.00),
    auc = c(mean = 0.15, min = 0.15, auc = 0.70, ic50 = 0.00),
    ic50 = c(mean = 0.20, min = 0.10, auc = 0.20, ic50 = 0.50)
  )

  observeEvent(input$tol_index_preset, {
    p <- input$tol_index_preset %||% "custom"
    if (!p %in% names(tol_index_presets)) return()
    w <- tol_index_presets[[p]]
    updateSliderInput(session, "w_tol_mean", value = as.numeric(w[["mean"]]))
    updateSliderInput(session, "w_tol_min", value = as.numeric(w[["min"]]))
    updateSliderInput(session, "w_tol_auc", value = as.numeric(w[["auc"]]))
    updateSliderInput(session, "w_tol_ic50", value = as.numeric(w[["ic50"]]))
  }, ignoreInit = TRUE)

  tol_index_metrics <- reactive({
    req(identical(input$tol_view, "broth"))
    df <- tol_dose_df()
    if (nrow(df) == 0) return(tibble::tibble())

    unit <- unique(na.omit(df$condition_unit))
    numeric_ok <- !all(is.na(df$condition_value)) && length(unit) == 1
    unit0 <- if (length(unit) == 1) unit[[1]] else NA_character_

    df |>
      group_by(strain) |>
      summarise(
        analyte = input$tol_analyte,
        rep_handling = input$tol_rep_handling %||% "mean",
        n_levels = n_distinct(.data$concentration),
        mean_response = if (all(is.na(.data$value_use))) NA_real_ else mean(.data$value_use, na.rm = TRUE),
        min_response = if (all(is.na(.data$value_use))) NA_real_ else min(.data$value_use, na.rm = TRUE),
        auc = if (numeric_ok) trapz_auc(.data$condition_value, .data$value_use) else NA_real_,
        ic50 = if (numeric_ok) ic50_linear_interp(.data$condition_value, .data$value_use) else NA_real_,
        concentration_unit = unit0,
        numeric_concentration = numeric_ok,
        .groups = "drop"
      ) |>
      mutate(
        mean_response = ifelse(is.nan(.data$mean_response) | is.infinite(.data$mean_response), NA_real_, .data$mean_response),
        min_response = ifelse(is.nan(.data$min_response) | is.infinite(.data$min_response), NA_real_, .data$min_response),
        auc = ifelse(is.nan(.data$auc) | is.infinite(.data$auc), NA_real_, .data$auc),
        ic50 = ifelse(is.nan(.data$ic50) | is.infinite(.data$ic50), NA_real_, .data$ic50)
      )
  })

  tol_index_avail <- reactive({
    m <- tol_index_metrics()
    c(
      mean = nrow(m) > 0 && !all(is.na(m$mean_response)),
      min = nrow(m) > 0 && !all(is.na(m$min_response)),
      auc = nrow(m) > 0 && !all(is.na(m$auc)),
      ic50 = nrow(m) > 0 && !all(is.na(m$ic50))
    )
  })

  tol_index_weights_used <- reactive({
    avail <- tol_index_avail()
    w_in <- c(
      mean = input$w_tol_mean %||% 0.40,
      min = input$w_tol_min %||% 0.30,
      auc = input$w_tol_auc %||% 0.30,
      ic50 = input$w_tol_ic50 %||% 0.00
    )
    normalize_weights(w_in, avail)
  })

  tol_index_weights_tbl <- reactive({
    avail <- tol_index_avail()
    w_in <- c(
      mean = input$w_tol_mean %||% 0.40,
      min = input$w_tol_min %||% 0.30,
      auc = input$w_tol_auc %||% 0.30,
      ic50 = input$w_tol_ic50 %||% 0.00
    )
    w_used <- tol_index_weights_used()
    tibble::tibble(
      metric = names(w_in),
      available = as.logical(avail[names(w_in)]),
      weight_input = as.numeric(w_in),
      weight_used = as.numeric(w_used[names(w_in)])
    )
  })

  tol_index_scaled <- reactive({
    m <- tol_index_metrics()
    if (nrow(m) == 0) return(m)
    m |>
      mutate(
        s_mean = scale01_vec(.data$mean_response),
        s_min = scale01_vec(.data$min_response),
        s_auc = scale01_vec(.data$auc),
        s_ic50 = scale01_vec(.data$ic50),
        s_mean = ifelse(is.na(.data$s_mean), 0, .data$s_mean),
        s_min = ifelse(is.na(.data$s_min), 0, .data$s_min),
        s_auc = ifelse(is.na(.data$s_auc), 0, .data$s_auc),
        s_ic50 = ifelse(is.na(.data$s_ic50), 0, .data$s_ic50)
      )
  })

  tol_index_rank_df <- reactive({
    df <- tol_index_scaled()
    if (nrow(df) == 0) return(df)
    w <- tol_index_weights_used()
    df |>
      mutate(
        c_mean = w["mean"] * .data$s_mean,
        c_min = w["min"] * .data$s_min,
        c_auc = w["auc"] * .data$s_auc,
        c_ic50 = w["ic50"] * .data$s_ic50,
        score = .data$c_mean + .data$c_min + .data$c_auc + .data$c_ic50,
        w_mean = w["mean"],
        w_min = w["min"],
        w_auc = w["auc"],
        w_ic50 = w["ic50"],
        note = "Exploratory composite score (technical duplicates; no biological replication)."
      ) |>
      arrange(desc(.data$score)) |>
      mutate(rank = dplyr::dense_rank(desc(.data$score))) |>
      select(rank, everything())
  })

  tol_index_plot_obj <- reactive({
    req(identical(input$tol_view, "broth"))
    df <- tol_index_rank_df()
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))
    ggplot(df, aes(x = reorder(strain, score), y = score, fill = strain)) +
      geom_col(width = 0.75, color = "white") +
      coord_flip() +
      scale_fill_manual(values = strain_palette(df$strain), guide = "none") +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
      labs(
        title = paste0("Composite tolerance score (broth): ", input$tol_analyte),
        subtitle = paste0("Replicate handling: ", input$tol_rep_handling %||% "mean"),
        x = "Strain",
        y = "Score (0-1)"
      ) +
      khadijat_theme()
  })

  output$tol_index_plot <- renderPlot({ tol_index_plot_obj() })

  output$tol_index_table <- renderDT({
    req(identical(input$tol_view, "broth"))
    df <- tol_index_rank_df()
    df_disp <- df |>
      select(
        rank, strain, score,
        mean_response, min_response, auc, ic50,
        w_mean, w_min, w_auc, w_ic50,
        numeric_concentration, concentration_unit
      )
    datatable(
      df_disp,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df_disp, digits = 4L))
    )
  })

  output$tol_index_weights <- renderDT({
    req(identical(input$tol_view, "broth"))
    df <- tol_index_weights_tbl()
    datatable(
      df,
      rownames = FALSE,
      options = list(dom = "t", paging = FALSE, searching = FALSE, info = FALSE, ordering = FALSE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  tol_index_sensitivity_df <- reactive({
    base <- tol_index_scaled()
    if (nrow(base) == 0) return(base)
    avail <- tol_index_avail()
    presets <- c(names(tol_index_presets), "current")

    out <- list()
    k <- 0L
    for (p in presets) {
      w_in <- if (identical(p, "current")) {
        c(mean = input$w_tol_mean, min = input$w_tol_min, auc = input$w_tol_auc, ic50 = input$w_tol_ic50)
      } else {
        tol_index_presets[[p]]
      }
      w <- normalize_weights(w_in, avail)
      tmp <- base |>
        mutate(
          score = w["mean"] * .data$s_mean + w["min"] * .data$s_min + w["auc"] * .data$s_auc + w["ic50"] * .data$s_ic50,
          preset = p,
          w_mean = w["mean"],
          w_min = w["min"],
          w_auc = w["auc"],
          w_ic50 = w["ic50"]
        ) |>
        arrange(desc(.data$score)) |>
        mutate(rank = dplyr::dense_rank(desc(.data$score))) |>
        select(preset, rank, strain, score, w_mean, w_min, w_auc, w_ic50)
      k <- k + 1L
      out[[k]] <- tmp
    }

    dplyr::bind_rows(out) |>
      mutate(preset = factor(.data$preset, levels = presets))
  })

  tol_index_sensitivity_plot_obj <- reactive({
    req(identical(input$tol_view, "broth"))
    df <- tol_index_sensitivity_df()
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))
    ggplot(df, aes(x = preset, y = rank, group = strain, color = strain)) +
      geom_line(linewidth = 0.8, alpha = 0.85) +
      geom_point(size = 2.2, alpha = 0.9) +
      scale_color_manual(values = strain_palette(df$strain)) +
      scale_y_reverse(breaks = sort(unique(df$rank))) +
      labs(
        title = "Rank sensitivity across weight presets (exploratory)",
        x = "Preset",
        y = "Rank (1 = best)"
      ) +
      khadijat_theme() +
      theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"))
  })

  output$tol_index_sensitivity_plot <- renderPlot({ tol_index_sensitivity_plot_obj() })

  output$tol_index_sensitivity_table <- renderDT({
    req(identical(input$tol_view, "broth"))
    df <- tol_index_sensitivity_df()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 12, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  output$tol_index_metrics_table <- renderDT({
    req(identical(input$tol_view, "broth"))
    df <- tol_index_metrics()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  output$dl_tol_index_plot <- downloadHandler(
    filename = function() paste0("tolerance_index_", standardize_key(input$tol_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(tol_index_plot_obj(), file, width = 10, height = 6.5, dpi = 320)
  )

  output$dl_tol_index_table <- downloadHandler(
    filename = function() paste0("tolerance_index_", standardize_key(input$tol_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) {
      writexl::write_xlsx(
        list(
          metrics = tol_index_metrics(),
          weights = tol_index_weights_tbl(),
          ranking = tol_index_rank_df(),
          sensitivity = tol_index_sensitivity_df()
        ),
        path = file
      )
    }
  )

    ## Objective 3: Enzymes (restricted to production rate)
    enz_category_fixed <- "Enzyme production rate"

    enz_raw <- reactive({
      d <- data_all()
      d$enzymes |>
        dplyr::filter(.data$category == enz_category_fixed)
    })

    output$enz_select_ui <- renderUI({
      df <- enz_raw()
      measures <- sort(unique(df$measure))
      if (length(measures) == 0) {
        return(div(class = "app-note", tags$p("Not available in current dataset.")))
      }
      tagList(
        div(class = "app-note", tags$p(tags$b("Category:"), " Enzyme production rate (fixed).")),
        selectInput("enz_measure", "Measure", choices = measures, selected = measures[[1]] %||% NULL),
        hr(),
        h4("Technical duplicates"),
        selectInput(
          "enz_rep_handling",
          "Replicate handling",
          choices = c(
            "Mean(S-1, S-dup) (Recommended)" = "mean",
            "S-1 only" = "s1",
            "S-dup only" = "sdup",
            "Difference (S-1 − S-dup) (QC)" = "diff"
          ),
          selected = "mean"
        ),
        checkboxInput("enz_show_points", "Show technical duplicates (S-1/S-dup)", value = TRUE)
      )
    })

    enz_selected_raw <- reactive({
      df <- enz_raw()
      req(input$enz_measure)
      df |>
        dplyr::filter(.data$measure == input$enz_measure)
    })

    enz_df_use <- reactive({
      df <- enz_selected_raw()
      if (nrow(df) == 0) return(tibble::tibble())
      apply_techrep_handling(
        df,
        group_cols = c("strain", "measure"),
        handling = input$enz_rep_handling %||% "mean",
        value_col = "value"
      )
    })

    enz_plot_obj <- reactive({
      df_line <- enz_df_use()
      df_raw <- enz_selected_raw()
      if (nrow(df_raw) == 0) return(ggplot() + labs(title = "Not available in current dataset"))

      rep_handling <- input$enz_rep_handling %||% "mean"
      pal <- strain_palette(df_line$strain)

      p <- ggplot(df_line, aes(x = reorder(strain, value_use), y = value_use, fill = strain)) +
        geom_col(width = 0.75, color = "white") +
        coord_flip() +
        scale_fill_manual(values = pal, guide = "none") +
        labs(
          title = paste0("Enzymes: ", enz_category_fixed, " | ", input$enz_measure),
          subtitle = paste0("Replicate handling: ", rep_handling),
          x = "Strain",
          y = "Value"
        ) +
        khadijat_theme()

      if (!isTRUE(input$enz_show_points)) return(p)
      if (rep_handling != "mean") return(p)

      p +
        geom_point(
          data = df_raw,
          aes(x = strain, y = value, shape = tech_rep),
          inherit.aes = FALSE,
          size = 2.2,
          alpha = 0.65,
          color = "#1d2733"
        ) +
        scale_shape_manual(values = c("S-1" = 16, "S-dup" = 1)) +
        guides(shape = guide_legend(title = "Tech rep"))
    })

    output$enz_plot <- renderPlot({ enz_plot_obj() })

    enz_table_df <- reactive({
      df <- enz_raw()
      if (nrow(df) == 0) return(tibble::tibble())

      wide <- df |>
        dplyr::select(strain, measure, tech_rep, value = .data$value) |>
        dplyr::group_by(.data$strain, .data$measure, .data$tech_rep) |>
        dplyr::summarise(
          value = if (all(is.na(.data$value))) NA_real_ else mean(.data$value, na.rm = TRUE),
          .groups = "drop"
        ) |>
        tidyr::pivot_wider(names_from = tech_rep, values_from = value)

      wide |>
        dplyr::rowwise() |>
        dplyr::mutate(
          category = enz_category_fixed,
          mean_value = mean(dplyr::c_across(dplyr::any_of(c("S-1", "S-dup"))), na.rm = TRUE),
          sd_value = {
            x <- dplyr::c_across(dplyr::any_of(c("S-1", "S-dup")))
            if (sum(!is.na(x)) < 2) NA_real_ else stats::sd(x, na.rm = TRUE)
          },
          diff_s1_sdup = .data$`S-1` - .data$`S-dup`,
          cv_pct = dplyr::if_else(!is.na(.data$mean_value) & .data$mean_value != 0, 100 * .data$sd_value / .data$mean_value, NA_real_),
          n_detected = sum(!is.na(dplyr::c_across(dplyr::any_of(c("S-1", "S-dup")))))
        ) |>
        dplyr::ungroup() |>
        dplyr::arrange(.data$measure, dplyr::desc(.data$mean_value))
    })

    output$enz_table <- renderDT({
      df <- enz_table_df()
      datatable(
        df,
        rownames = FALSE,
        options = list(pageLength = 12, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
      )
    })

    enz_heatmap_df <- reactive({
      df <- enz_raw()
      if (nrow(df) == 0) return(tibble::tibble())
      rep_handling <- input$enz_rep_handling %||% "mean"
      df_use <- apply_techrep_handling(df, group_cols = c("strain", "measure"), handling = rep_handling, value_col = "value") |>
        dplyr::transmute(strain, measure, value = value_use)

      if (identical(input$enz_heatmap_scale, "z")) {
        df_use |>
          dplyr::group_by(.data$measure) |>
          dplyr::mutate(
            value = {
              m <- mean(.data$value, na.rm = TRUE)
              s <- stats::sd(.data$value, na.rm = TRUE)
              if (is.na(s) || s == 0) 0 else (.data$value - m) / s
            }
          ) |>
          dplyr::ungroup()
      } else {
        df_use
      }
    })

    enz_heatmap_plot_obj <- reactive({
      df <- enz_heatmap_df()
      if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))
      ggplot(df, aes(x = measure, y = strain, fill = value)) +
        geom_tile(color = "white", linewidth = 0.25) +
        scale_y_discrete(limits = rev(sort(unique(df$strain)))) +
        scale_fill_viridis_c(option = "C", na.value = "#f1f1f1") +
        labs(
          title = paste0("Enzyme production rate heatmap (", input$enz_heatmap_scale %||% "z", ")"),
          x = "Measure",
          y = "Strain",
          fill = if (identical(input$enz_heatmap_scale, "z")) "Z-score" else "Value"
        ) +
	        khadijat_theme() +
	        theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"))
    })

    output$enz_heatmap <- renderPlot({ enz_heatmap_plot_obj() })

    enz_pvals_all <- reactive({
      req(isTRUE(input$enz_enable_pvals))
      df <- enz_raw() |>
        dplyr::filter(!is.na(.data$value))

      method <- input$enz_pval_method %||% "aov_tukey"
      scope <- input$enz_pval_scope %||% "selected"
      measures <- if (identical(scope, "all")) sort(unique(df$measure)) else (input$enz_measure %||% character())

      out <- list()
      k <- 0L
      for (m in measures) {
        dat <- df |>
          dplyr::filter(.data$measure == m)
        if (nrow(dat) == 0) next
        res <- if (identical(method, "kw_wilcox")) safe_kw_wilcox(dat, response_col = "value", group_col = "strain") else safe_aov_tukey(dat, response_col = "value", group_col = "strain")
        if (nrow(res$pairwise) == 0) next
        k <- k + 1L
        out[[k]] <- res$pairwise |>
          dplyr::mutate(
            category = enz_category_fixed,
            measure = m,
            global_p = res$global_p,
            p_adj_within = .data$p_adj,
            method = method
          )
      }
      dplyr::bind_rows(out)
    })

    output$enz_pval_filter_ui <- renderUI({
      req(isTRUE(input$enz_enable_pvals))
      scope <- input$enz_pval_scope %||% "selected"
      if (!identical(scope, "all")) return(NULL)
      df <- enz_raw()
      measures <- sort(unique(df$measure))
      selectInput("enz_pval_measure", "Filter measure", choices = c("All" = "__all__", measures), selected = "__all__")
    })

    enz_pvals_filtered <- reactive({
      req(isTRUE(input$enz_enable_pvals))
      padj_overall <- input$enz_padj_overall %||% "BH"
      df <- enz_pvals_all()
      if (nrow(df) == 0) return(df)

      if (!is.null(input$enz_pval_measure) && input$enz_pval_measure != "__all__") {
        df <- df |>
          dplyr::filter(.data$measure == input$enz_pval_measure)
      }

      df |>
        dplyr::mutate(p_adj_overall = stats::p.adjust(.data$p_adj_within, method = padj_overall)) |>
        dplyr::arrange(.data$p_adj_overall)
    })

    output$enz_pvals_table <- renderDT({
      req(isTRUE(input$enz_enable_pvals))
      df <- enz_pvals_filtered()
      datatable(
        df,
        rownames = FALSE,
        options = list(pageLength = 15, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
      )
    })

    output$dl_enz_plot <- downloadHandler(
      filename = function() paste0("enzymes_production_rate_", standardize_key(input$enz_measure), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) save_plot(enz_plot_obj(), file, width = 10, height = 5.5, dpi = 320)
    )

    output$dl_enz_heatmap <- downloadHandler(
      filename = function() paste0("enzymes_production_rate_heatmap_", input$enz_heatmap_scale %||% "z", "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) save_plot(enz_heatmap_plot_obj(), file, width = 12, height = 6.5, dpi = 320)
    )

    output$dl_enz_pvals <- downloadHandler(
      filename = function() paste0("enzymes_production_rate_pvalues_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        readr::write_csv(enz_pvals_filtered(), file, na = "")
      }
    )

    output$dl_enz_table <- downloadHandler(
      filename = function() paste0("enzymes_production_rate_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
      content = function(file) writexl::write_xlsx(list(summary = enz_table_df()), path = file)
    )

    ## Objective 4: Volatiles (GC) - exploratory only (no replicates)
  gc_tidy <- reactive({
    d <- data_all()
    if (nrow(d$gc) == 0) return(tibble::tibble())
    gc_apply_duplicate_policy(d$gc, policy = input$gc_dup_policy %||% "sum")
  })

  output$gc_heatmap_strain_ui <- renderUI({
    g <- gc_tidy()
    if (nrow(g) == 0) return(NULL)
    strains <- sort(unique(g$strain))
    selectInput(
      "gc_heatmap_strain",
      "Heatmap strain focus",
      choices = c("All strains" = "__all__", strains),
      selected = "__all__"
    )
  })

  gc_heatmap_data <- reactive({
    g <- gc_tidy()
    req(nrow(g) > 0)

    focus <- input$gc_heatmap_strain %||% "__all__"
    g_focus <- if (identical(focus, "__all__")) g else dplyr::filter(g, .data$strain == focus)

    top <- g_focus |>
      group_by(compound_key, compound_std) |>
      summarise(m = if (all(is.na(.data$ra_ug_L))) NA_real_ else mean(.data$ra_ug_L, na.rm = TRUE), .groups = "drop") |>
      mutate(m = ifelse(is.nan(.data$m), NA_real_, .data$m)) |>
      arrange(desc(.data$m)) |>
      slice_head(n = as.integer(input$gc_top_n %||% 20L))

    out <- g |>
      semi_join(top, by = c("compound_key", "compound_std"))
    if (!identical(focus, "__all__")) {
      out <- out |>
        dplyr::filter(.data$strain == focus)
    }

    out <- out |>
      group_by(strain, time_h, compound_std) |>
      summarise(
        ra_ug_L = if (all(is.na(.data$ra_ug_L))) NA_real_ else mean(.data$ra_ug_L, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(ra_ug_L = ifelse(is.nan(.data$ra_ug_L) | is.infinite(.data$ra_ug_L), NA_real_, .data$ra_ug_L))

    # Stable compound ordering (high -> low) to improve readability.
    ord <- out |>
      group_by(compound_std) |>
      summarise(m = if (all(is.na(.data$ra_ug_L))) 0 else mean(.data$ra_ug_L, na.rm = TRUE), .groups = "drop") |>
      arrange(.data$m) |>
      pull(.data$compound_std)
    out |>
      mutate(compound_std = factor(.data$compound_std, levels = ord))
  })

  gc_heatmap_plot <- reactive({
    df <- gc_heatmap_data()
    req(nrow(df) > 0)

    focus <- input$gc_heatmap_strain %||% "__all__"
    ncol <- as.integer(input$gc_heatmap_ncol %||% 2)
    n_comp <- length(unique(df$compound_std))
    y_text <- dplyr::case_when(
      n_comp >= 80 ~ 4.5,
      n_comp >= 60 ~ 5,
      n_comp >= 40 ~ 5.5,
      n_comp >= 25 ~ 6,
      TRUE ~ 7
    )

    p <- ggplot(df, aes(x = factor(time_h), y = compound_std, fill = ra_ug_L)) +
      geom_tile(color = "white", linewidth = 0.12) +
      scale_fill_viridis_c(option = "C", na.value = "#f1f1f1", trans = "sqrt") +
      labs(
        title = paste0("Top volatiles heatmap (policy: ", input$gc_dup_policy, ")"),
        x = "Time (h)",
        y = "Compound",
        fill = "RA (ug/L)\n(sqrt scale)"
      ) +
      khadijat_theme(base_size = 11) +
      theme(
        axis.text.y = element_text(size = y_text, face = "bold"),
        strip.text = element_text(size = 11, face = "bold")
      )

    if (identical(focus, "__all__")) {
      p <- p + facet_wrap(~strain, ncol = ncol, scales = "free_y")
    }
    p
  })

  output$gc_heatmap <- renderPlot({ gc_heatmap_plot() })

  output$gc_compare_ui <- renderUI({
    g <- gc_tidy()
    if (nrow(g) == 0) {
      return(div(class = "app-note", tags$p("GC dataset not available in current dataset.")))
    }
    strains <- sort(unique(g$strain))
    times <- sort(unique(g$time_h))

    comp_rank <- g |>
      group_by(compound_std, compound_key) |>
      summarise(m = if (all(is.na(.data$ra_ug_L))) NA_real_ else mean(.data$ra_ug_L, na.rm = TRUE), .groups = "drop") |>
      mutate(m = ifelse(is.nan(.data$m), NA_real_, .data$m)) |>
      arrange(desc(.data$m))
    comp_choices <- comp_rank$compound_std
    comp_default <- head(comp_choices, 5)

    tagList(
      fluidRow(
        column(
          4,
          radioButtons(
            "gc_compare_facet",
            "Facet by",
            choices = c("Compound (Recommended)" = "compound", "Strain" = "strain"),
            selected = "compound",
            inline = TRUE
          )
        ),
        column(
          4,
          selectInput(
            "gc_compare_transform",
            "Transform",
            choices = c("sqrt (Recommended)" = "sqrt", "log1p" = "log1p", "None" = "none"),
            selected = "sqrt"
          )
        ),
        column(
          4,
          checkboxInput(
            "gc_compare_zero_fill",
            "Show non-detects as 0 (plot only)",
            value = TRUE
          )
        )
      ),
      fluidRow(
        column(
          6,
          selectInput(
            "gc_compare_strains",
            "Strains",
            choices = strains,
            selected = strains,
            multiple = TRUE
          )
        ),
        column(
          6,
          selectizeInput(
            "gc_compare_compounds",
            "Compounds",
            choices = comp_choices,
            selected = comp_default,
            multiple = TRUE,
            options = list(plugins = list("remove_button"), maxOptions = 500)
          )
        )
      ),
      fluidRow(
        column(
          6,
          selectInput(
            "gc_snapshot_time",
            "Snapshot time (h)",
            choices = times,
            selected = times[[length(times)]] %||% NULL
          )
        ),
        column(
          6,
          selectInput(
            "gc_snapshot_display",
            "Snapshot display",
            choices = c("Stacked bars" = "stack", "Dodged bars" = "dodge"),
            selected = "stack"
          )
        )
      ),
      div(
        class = "app-note",
        tags$p(tags$b("Note:"), "GC has no independent replicates in the current dataset; these comparisons are descriptive/exploratory only (no p-values/CIs).")
      )
    )
  })

  gc_compare_df <- reactive({
    g <- gc_tidy()
    req(nrow(g) > 0)

    strains_all <- sort(unique(g$strain))
    comp_rank <- g |>
      group_by(compound_std) |>
      summarise(m = if (all(is.na(.data$ra_ug_L))) NA_real_ else mean(.data$ra_ug_L, na.rm = TRUE), .groups = "drop") |>
      mutate(m = ifelse(is.nan(.data$m), NA_real_, .data$m)) |>
      arrange(desc(.data$m))
    comps_all <- comp_rank$compound_std

    strains_sel <- input$gc_compare_strains %||% strains_all
    comps_sel <- input$gc_compare_compounds %||% head(comps_all, 5)

    df <- g |>
      filter(.data$strain %in% strains_sel, .data$compound_std %in% comps_sel) |>
      group_by(strain, time_h, compound_std) |>
      summarise(
        ra_ug_L = if (all(is.na(.data$ra_ug_L))) NA_real_ else mean(.data$ra_ug_L, na.rm = TRUE),
        detected = any(.data$detected, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(ra_ug_L = ifelse(is.nan(.data$ra_ug_L) | is.infinite(.data$ra_ug_L), NA_real_, .data$ra_ug_L))

    v <- df$ra_ug_L
    if (isTRUE(input$gc_compare_zero_fill %||% TRUE)) {
      v <- dplyr::coalesce(v, 0)
    }

    tr <- input$gc_compare_transform %||% "sqrt"
    v_plot <- dplyr::case_when(
      identical(tr, "log1p") ~ log1p(pmax(v, 0)),
      identical(tr, "sqrt") ~ sqrt(pmax(v, 0)),
      TRUE ~ v
    )

    df |>
      mutate(value_plot = v_plot, transform = tr)
  })

  gc_compare_plot_obj <- reactive({
    df <- gc_compare_df()
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))

    facet_by <- input$gc_compare_facet %||% "compound"
    tr <- df$transform[[1]] %||% "sqrt"

    ylab <- dplyr::case_when(
      identical(tr, "log1p") ~ "log1p(RA (ug/L))",
      identical(tr, "sqrt") ~ "sqrt(RA (ug/L))",
      TRUE ~ "RA (ug/L)"
    )

    if (identical(facet_by, "strain")) {
      ggplot(df, aes(x = time_h, y = value_plot, color = compound_std, group = compound_std)) +
        geom_line(linewidth = 0.7, alpha = 0.85) +
        geom_point(size = 2, alpha = 0.9) +
        facet_wrap(~strain, scales = "free_y") +
        scale_color_viridis_d(option = "C") +
        scale_x_continuous(breaks = sort(unique(df$time_h))) +
        labs(
          title = "GC trajectories (selected compounds/strains; descriptive)",
          x = "Time (h)",
          y = ylab,
          color = "Compound"
        ) +
        khadijat_theme()
    } else {
      ggplot(df, aes(x = time_h, y = value_plot, color = strain, group = strain)) +
        geom_line(linewidth = 0.8, alpha = 0.85) +
        geom_point(size = 2.1, alpha = 0.9) +
        facet_wrap(~compound_std, scales = "free_y") +
        scale_color_manual(values = strain_palette(df$strain)) +
        scale_x_continuous(breaks = sort(unique(df$time_h))) +
        labs(
          title = "GC trajectories (selected compounds/strains; descriptive)",
          x = "Time (h)",
          y = ylab,
          color = "Strain"
        ) +
        khadijat_theme()
    }
  })

  output$gc_compare_plot <- renderPlot({ gc_compare_plot_obj() })

  gc_snapshot_plot_obj <- reactive({
    df <- gc_compare_df()
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))

    tp <- as.integer(input$gc_snapshot_time %||% max(df$time_h, na.rm = TRUE))
    df <- df |>
      filter(.data$time_h == tp)
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))

    tr <- df$transform[[1]] %||% "sqrt"
    ylab <- dplyr::case_when(
      identical(tr, "log1p") ~ "log1p(RA (ug/L))",
      identical(tr, "sqrt") ~ "sqrt(RA (ug/L))",
      TRUE ~ "RA (ug/L)"
    )

    pos <- if (identical(input$gc_snapshot_display %||% "stack", "dodge")) ggplot2::position_dodge2(width = 0.8, padding = 0.1) else "stack"

    ggplot(df, aes(x = strain, y = value_plot, fill = compound_std)) +
      geom_col(position = pos, width = 0.78, color = "white", linewidth = 0.2) +
      scale_fill_viridis_d(option = "C") +
      labs(
        title = paste0("Snapshot at ", tp, " h (descriptive)"),
        x = "Strain",
        y = ylab,
        fill = "Compound"
      ) +
	      khadijat_theme() +
	      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"))
  })

  output$gc_snapshot_plot <- renderPlot({ gc_snapshot_plot_obj() })

  gc_table_df <- reactive({
    g <- gc_tidy()
    if (nrow(g) == 0) return(tibble::tibble())
    g |>
      select(strain, time_h, compound_std, ra_ug_L, detected)
  })

  output$gc_table <- renderDT({
    df <- gc_table_df()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  output$gc_pca_ui <- renderUI({
    g <- gc_tidy()
    if (nrow(g) == 0) return(NULL)
    tagList(
      fluidRow(
        column(
          4,
          selectInput(
            "gc_pca_transform",
            "PCA transform",
            choices = c("log1p (Recommended)" = "log1p", "sqrt" = "sqrt", "None" = "none"),
            selected = "log1p"
          )
        ),
        column(
          4,
          selectInput(
            "gc_pca_features",
            "Features used",
            choices = c("All compounds" = "all", "Top by mean RA" = "top"),
            selected = "all"
          )
        ),
        column(4, checkboxInput("gc_pca_scale", "Scale features", value = TRUE))
      ),
      conditionalPanel(
        condition = "input.gc_pca_features === 'top'",
        sliderInput("gc_pca_top_n", "Top compounds (mean RA)", min = 10, max = 200, value = 50, step = 5)
      ),
      fluidRow(
        column(
          4,
          selectInput(
            "gc_pca_color_by",
            "Color by",
            choices = c("Strain" = "strain", "Time" = "time"),
            selected = "strain"
          )
        ),
        column(8, div(class = "app-note", tags$p("PCA uses 0-fill for non-detects for ordination only. No inferential claims.")))
      )
    )
  })

  gc_pca_fit <- reactive({
    g <- gc_tidy()
    req(nrow(g) > 0)

    df <- g |>
      mutate(feature = .data$compound_key) |>
      group_by(strain, time_h, feature) |>
      summarise(
        v = if (all(is.na(.data$ra_ug_L))) NA_real_ else mean(.data$ra_ug_L, na.rm = TRUE),
        compound_std = dplyr::first(.data$compound_std),
        .groups = "drop"
      ) |>
      mutate(v = ifelse(is.nan(.data$v) | is.infinite(.data$v), NA_real_, .data$v))

    feat_mode <- input$gc_pca_features %||% "all"
    if (identical(feat_mode, "top")) {
      top_n <- as.integer(input$gc_pca_top_n %||% 50L)
      top_tbl <- df |>
        group_by(feature, compound_std) |>
        summarise(m = if (all(is.na(.data$v))) NA_real_ else mean(.data$v, na.rm = TRUE), .groups = "drop") |>
        mutate(m = ifelse(is.nan(.data$m), NA_real_, .data$m)) |>
        arrange(desc(.data$m)) |>
        slice_head(n = top_n)
      df <- df |>
        semi_join(top_tbl, by = c("feature", "compound_std"))
    }

    wide <- df |>
      select(strain, time_h, feature, v) |>
      tidyr::pivot_wider(names_from = feature, values_from = v, values_fill = 0)

    mat <- as.matrix(wide[, setdiff(names(wide), c("strain", "time_h"))])
    mat[!is.finite(mat)] <- 0

    # Transform before PCA (common in compositional/abundance data).
    tr <- input$gc_pca_transform %||% "log1p"
    if (identical(tr, "sqrt")) mat <- sqrt(pmax(mat, 0))
    if (identical(tr, "log1p")) mat <- log1p(pmax(mat, 0))

    # Drop constant features to avoid scale() issues.
    sds <- apply(mat, 2, stats::sd)
    keep <- is.finite(sds) & sds > 0
    mat <- mat[, keep, drop = FALSE]
    if (ncol(mat) < 2 || nrow(mat) < 2) return(NULL)

    scale_flag <- isTRUE(input$gc_pca_scale %||% TRUE)
    pca <- stats::prcomp(mat, center = TRUE, scale. = scale_flag)
    ve <- (pca$sdev^2) / sum(pca$sdev^2)

    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    scores$strain <- wide$strain
    scores$time_h <- wide$time_h

    load <- as.data.frame(pca$rotation[, 1:2, drop = FALSE])
    load$feature <- rownames(load)
    rownames(load) <- NULL
    feat_map <- df |>
      distinct(feature, compound_std)
    load <- load |>
      left_join(feat_map, by = "feature") |>
      mutate(compound_std = .data$compound_std %||% .data$feature)

    list(
      pca = pca,
      scores = scores,
      var_explained = ve,
      loadings = load,
      transform = tr,
      scale = scale_flag
    )
  })

  gc_pca_loadings_df <- reactive({
    fit <- gc_pca_fit()
    if (is.null(fit)) return(tibble::tibble())
    fit$loadings |>
      dplyr::transmute(
        compound = .data$compound_std,
        PC1 = .data$PC1,
        PC2 = .data$PC2,
        abs_PC1 = abs(.data$PC1),
        abs_PC2 = abs(.data$PC2),
        max_abs = pmax(.data$abs_PC1, .data$abs_PC2)
      ) |>
      arrange(desc(.data$max_abs)) |>
      slice_head(n = 30) |>
      select(-.data$max_abs)
  })

  gc_pca_plot <- reactive({
    fit <- gc_pca_fit()
    if (is.null(fit)) {
      return(ggplot() + labs(title = "PCA not available", subtitle = "Not enough non-constant features after preprocessing."))
    }

    scores <- fit$scores
    ve <- fit$var_explained

    xlab <- paste0("PC1 (", scales::percent(ve[[1]] %||% 0, accuracy = 0.1), ")")
    ylab <- paste0("PC2 (", scales::percent(ve[[2]] %||% 0, accuracy = 0.1), ")")

    base_shapes <- c(16, 17, 15, 18, 3, 8, 4, 0, 1, 2)

    color_by <- input$gc_pca_color_by %||% "strain"
    if (identical(color_by, "time")) {
      strain_levels <- sort(unique(scores$strain))
      shape_map <- stats::setNames(base_shapes[seq_along(strain_levels)], strain_levels)
      ggplot(scores, aes(x = PC1, y = PC2, color = time_h, shape = factor(strain, levels = strain_levels))) +
        geom_point(size = 2.4, alpha = 0.85) +
        scale_color_viridis_c(option = "C") +
        scale_shape_manual(values = shape_map) +
        labs(
          title = "PCA of volatile profiles (exploratory)",
          subtitle = "0-fill for non-detects; no inferential claims.",
          x = xlab,
          y = ylab,
          shape = "Strain",
          color = "Time (h)"
        ) +
        khadijat_theme()
    } else {
      time_levels <- sort(unique(scores$time_h))
      shape_map <- stats::setNames(base_shapes[seq_along(time_levels)], as.character(time_levels))
      ggplot(scores, aes(x = PC1, y = PC2, color = strain, shape = factor(time_h, levels = time_levels))) +
        geom_point(size = 2.4, alpha = 0.85) +
        scale_color_manual(values = strain_palette(scores$strain)) +
        scale_shape_manual(values = shape_map) +
        labs(
          title = "PCA of volatile profiles (exploratory)",
          subtitle = "0-fill for non-detects; no inferential claims.",
          x = xlab,
          y = ylab,
          shape = "Time (h)",
          color = "Strain"
        ) +
        khadijat_theme()
    }
  })

  output$gc_pca <- renderPlot({ gc_pca_plot() })

  output$gc_pca_loadings <- renderDT({
    df <- gc_pca_loadings_df()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  gc_summary_df <- reactive({
    g <- gc_tidy()
    if (nrow(g) == 0) return(tibble::tibble())
    g |>
      group_by(strain, time_h) |>
      summarise(
        total_ra = sum(.data$ra_ug_L, na.rm = TRUE),
        richness = sum(.data$detected, na.rm = TRUE),
        .groups = "drop"
      )
  })

  gc_total_plot_obj <- reactive({
    df <- gc_summary_df()
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))
    pal <- strain_palette(df$strain)
    ggplot(df, aes(x = time_h, y = total_ra, color = strain, group = strain)) +
      geom_line(linewidth = 0.85, alpha = 0.9) +
      geom_point(size = 2) +
      scale_color_manual(values = pal) +
      scale_x_continuous(breaks = sort(unique(df$time_h))) +
      labs(
        title = "Total volatile load over time (descriptive)",
        x = "Time (h)",
        y = "Total RA (ug/L; detected values only)",
        color = "Strain"
      ) +
      khadijat_theme()
  })

  gc_richness_plot_obj <- reactive({
    df <- gc_summary_df()
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))
    pal <- strain_palette(df$strain)
    ggplot(df, aes(x = time_h, y = richness, color = strain, group = strain)) +
      geom_line(linewidth = 0.85, alpha = 0.9) +
      geom_point(size = 2) +
      scale_color_manual(values = pal) +
      scale_x_continuous(breaks = sort(unique(df$time_h))) +
      labs(
        title = "Compound richness over time (descriptive)",
        x = "Time (h)",
        y = "# compounds detected (>0)",
        color = "Strain"
      ) +
      khadijat_theme()
  })

  output$gc_total_plot <- renderPlot({ gc_total_plot_obj() })
  output$gc_richness_plot <- renderPlot({ gc_richness_plot_obj() })

  output$gc_map_coverage <- renderDT({
    d <- data_all()
    g <- gc_tidy()
    m <- d$compound_map

    if (nrow(g) == 0) {
      return(datatable(tibble::tibble(note = "GC dataset not available"), options = list(dom = "tip")))
    }
    if (nrow(m) == 0) {
      return(datatable(tibble::tibble(note = "Mapping file not available or empty"), options = list(dom = "tip")))
    }

    compounds <- g |>
      distinct(compound_key, compound_std) |>
      left_join(m |> distinct(compound_key, compound_class, common_name, sensory_note), by = "compound_key") |>
      mutate(mapped = !is.na(compound_class) | !is.na(common_name))

    cov_tbl <- tibble::tibble(
      total_compounds = nrow(compounds),
      mapped_compounds = sum(compounds$mapped, na.rm = TRUE),
      mapping_coverage = mapped_compounds / total_compounds
    )

    datatable(
      cov_tbl,
      rownames = FALSE,
      options = list(dom = "tip", columnDefs = dt_coldefs_numeric(cov_tbl, digits = 4L))
    )
  })

  output$dl_gc_heatmap <- downloadHandler(
    filename = function() paste0("gc_volatiles_heatmap_", input$gc_dup_policy, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(gc_heatmap_plot(), file, width = 10, height = 9, dpi = 320)
  )

  output$dl_gc_compare <- downloadHandler(
    filename = function() paste0("gc_volatiles_compare_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(gc_compare_plot_obj(), file, width = 11, height = 7.5, dpi = 320)
  )

  output$dl_gc_pca <- downloadHandler(
    filename = function() paste0("gc_volatiles_pca_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(gc_pca_plot(), file, width = 10, height = 6.5, dpi = 320)
  )

  output$dl_gc_loadings <- downloadHandler(
    filename = function() paste0("gc_volatiles_pca_loadings_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) readr::write_csv(gc_pca_loadings_df(), file, na = "")
  )

  output$dl_gc_table <- downloadHandler(
    filename = function() paste0("gc_volatiles_tidy_", input$gc_dup_policy, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) writexl::write_xlsx(list(table = gc_table_df()), path = file)
  )

  ## Objective 5: Metabolic dynamics (mono-culture; descriptive)
  output$met_select_ui <- renderUI({
    d <- data_all()
    src <- input$met_source %||% "met"
    if (identical(src, "aut")) {
      choices <- sort(unique(d$autolytic$analyte))
      selectInput("met_analyte", "Analyte (Autolytic-full)", choices = choices, selected = choices[[1]] %||% NULL)
    } else if (identical(src, "red24")) {
      choices <- sort(unique(d$reducing$analyte))
      selectInput("met_analyte", "Analyte (Reducing sugar-24h)", choices = choices, selected = choices[[1]] %||% NULL)
    } else {
      analytes <- d$metabolites |>
        distinct(type_unit, analyte) |>
        arrange(type_unit, analyte)
      choices <- paste0(analytes$type_unit, " | ", analytes$analyte)
      selectInput("met_analyte", "Analyte (Yeast metabolites)", choices = choices, selected = choices[[1]] %||% NULL)
    }
  })

  met_data <- reactive({
    req(input$met_analyte)
    d <- data_all()
    src <- input$met_source %||% "met"
    rep_handling <- input$met_rep_handling %||% "mean"

    if (identical(src, "aut")) {
      analyte_sel <- input$met_analyte
      df <- d$autolytic |>
        filter(.data$analyte == analyte_sel)
      y_axis <- analyte_sel
      title <- paste0("Autolytic dynamics: ", analyte_sel)
      group_cols <- c("strain", "time_h", "analyte")
      list(raw = df, src = src, type_unit = NA_character_, analyte = analyte_sel, title = title, y_axis = y_axis, group_cols = group_cols)
    } else if (identical(src, "red24")) {
      analyte_sel <- input$met_analyte
      df <- d$reducing |>
        filter(.data$analyte == analyte_sel)
      y_axis <- analyte_sel
      title <- paste0("Reducing sugar (24h): ", analyte_sel)
      group_cols <- c("strain", "time_h", "analyte")
      list(raw = df, src = src, type_unit = NA_character_, analyte = analyte_sel, title = title, y_axis = y_axis, group_cols = group_cols)
    } else {
      parts <- str_split_fixed(input$met_analyte, "\\s\\|\\s", 2)
      type_unit_sel <- parts[1, 1]
      analyte_sel <- parts[1, 2]

      df <- d$metabolites |>
        filter(.data$type_unit == type_unit_sel, .data$analyte == analyte_sel)

      y_axis <- paste0(analyte_sel, " (", type_unit_sel, ")")
      title <- paste0("Metabolite dynamics: ", analyte_sel, " (", type_unit_sel, ")")
      group_cols <- c("strain", "time_h", "type_unit", "analyte")
      list(raw = df, src = src, type_unit = type_unit_sel, analyte = analyte_sel, title = title, y_axis = y_axis, group_cols = group_cols)
    }
  })

  met_df_use <- reactive({
    md <- met_data()
    df <- md$raw
    if (nrow(df) == 0) return(tibble::tibble())
    handling <- input$met_rep_handling %||% "mean"
    apply_techrep_handling(
      df,
      group_cols = md$group_cols,
      handling = handling,
      value_col = "value"
    )
  })

  met_table_df <- reactive({
    md <- met_data()
    df <- md$raw
    if (nrow(df) == 0) return(tibble::tibble())
    summarise_tech_dups(df, value_col = "value", group_cols = md$group_cols) |>
      arrange(.data$strain, .data$time_h)
  })

  met_plot_obj <- reactive({
    md <- met_data()
    df_use <- met_df_use()
    if (nrow(df_use) == 0) return(ggplot() + labs(title = "Not available in current dataset"))

    rep_handling <- input$met_rep_handling %||% "mean"
    pal <- strain_palette(df_use$strain)
    ylab <- if (identical(rep_handling, "diff")) "Difference (S-1 - S-dup)" else (md$y_axis %||% "Value")

    p <- ggplot(df_use, aes(x = time_h, y = value_use, color = strain, group = strain)) +
      { if (length(unique(df_use$time_h)) >= 2) geom_line(linewidth = 0.85, alpha = 0.9) } +
      geom_point(size = 2) +
      scale_color_manual(values = pal) +
      scale_x_continuous(breaks = sort(unique(df_use$time_h))) +
      labs(
        title = md$title,
        subtitle = paste0("Replicate handling: ", rep_handling),
        x = "Time (h)",
        y = ylab,
        color = "Strain"
      ) +
      khadijat_theme()

    if (isTRUE(input$met_show_points %||% FALSE) && !identical(rep_handling, "diff")) {
      gd <- md$raw |>
        dplyr::filter(!is.na(.data$value)) |>
        dplyr::mutate(tech_rep = factor(.data$tech_rep, levels = c("S-1", "S-dup")))
      p <- p +
        geom_point(
          data = gd,
          aes(x = time_h, y = value, color = strain, shape = tech_rep),
          size = 1.6,
          alpha = 0.45,
          inherit.aes = FALSE
        ) +
        scale_shape_manual(values = c("S-1" = 16, "S-dup" = 1)) +
        guides(shape = guide_legend(title = "Tech rep"))
    }
    p
  })

  output$met_plot <- renderPlot({ met_plot_obj() })

  output$met_timepoint_ui <- renderUI({
    md <- met_data()
    df <- md$raw
    if (nrow(df) == 0) return(NULL)
    tp <- sort(unique(df$time_h))
    selectInput("met_snapshot_time", "Time (h)", choices = tp, selected = tp[[length(tp)]] %||% NULL)
  })

  output$met_delta_ui <- renderUI({
    md <- met_data()
    df <- md$raw
    if (nrow(df) == 0) return(NULL)
    tp <- sort(unique(df$time_h))
    if (length(tp) < 2) return(NULL)
    tagList(
      selectInput("met_delta_t1", "Delta t1 (h)", choices = tp, selected = tp[[1]] %||% NULL),
      selectInput("met_delta_t2", "Delta t2 (h)", choices = tp, selected = tp[[length(tp)]] %||% NULL)
    )
  })

  met_snapshot_plot_obj <- reactive({
    md <- met_data()
    df <- met_df_use()
    req(nrow(df) > 0)
    req(input$met_snapshot_time)
    t0 <- as.integer(input$met_snapshot_time)
    df_t <- df |>
      dplyr::filter(.data$time_h == t0)
    if (nrow(df_t) == 0) return(ggplot() + labs(title = "Selected timepoint not available"))

    pal <- strain_palette(df_t$strain)
    rep_handling <- input$met_rep_handling %||% "mean"

    if (identical(input$met_compare_geom %||% "dot", "bar")) {
      ggplot(df_t, aes(x = reorder(strain, value_use), y = value_use, fill = strain)) +
        geom_col(width = 0.8, color = "white") +
        coord_flip() +
        scale_fill_manual(values = pal, guide = "none") +
        labs(
          title = paste0(md$title, " (snapshot at ", t0, " h)"),
          subtitle = paste0("Replicate handling: ", rep_handling),
          x = "Strain",
          y = if (identical(rep_handling, "diff")) "Difference (S-1 - S-dup)" else (md$y_axis %||% "Value")
        ) +
        khadijat_theme()
    } else {
      ggplot(df_t, aes(x = reorder(strain, value_use), y = value_use, color = strain)) +
        geom_point(size = 3, alpha = 0.9) +
        coord_flip() +
        scale_color_manual(values = pal, guide = "none") +
        labs(
          title = paste0(md$title, " (snapshot at ", t0, " h)"),
          subtitle = paste0("Replicate handling: ", rep_handling),
          x = "Strain",
          y = if (identical(rep_handling, "diff")) "Difference (S-1 - S-dup)" else (md$y_axis %||% "Value")
        ) +
        khadijat_theme()
    }
  })

  output$met_snapshot_plot <- renderPlot({ met_snapshot_plot_obj() })

  met_delta_plot_obj <- reactive({
    md <- met_data()
    df <- met_df_use()
    req(nrow(df) > 0)
    req(input$met_delta_t1, input$met_delta_t2)
    t1 <- as.integer(input$met_delta_t1)
    t2 <- as.integer(input$met_delta_t2)
    if (is.na(t1) || is.na(t2) || t1 == t2) {
      return(ggplot() + labs(title = "Select two different timepoints"))
    }

    w <- df |>
      dplyr::select(strain, time_h, value_use) |>
      tidyr::pivot_wider(names_from = time_h, values_from = value_use)

    if (!(as.character(t1) %in% names(w)) || !(as.character(t2) %in% names(w))) {
      return(ggplot() + labs(title = "Selected timepoints not available"))
    }

    w <- w |>
      dplyr::mutate(delta = .data[[as.character(t2)]] - .data[[as.character(t1)]])

    pal <- strain_palette(w$strain)
    rep_handling <- input$met_rep_handling %||% "mean"
    ggplot(w, aes(x = reorder(strain, delta), y = delta, fill = strain)) +
      geom_col(width = 0.8, color = "white") +
      coord_flip() +
      scale_fill_manual(values = pal, guide = "none") +
      labs(
        title = paste0(md$title, " (delta: ", t2, " h ? ", t1, " h)"),
        subtitle = paste0("Replicate handling: ", rep_handling),
        x = "Strain",
        y = paste0("Delta (", if (identical(rep_handling, "diff")) "Difference (S-1 - S-dup)" else (md$y_axis %||% "Value"), ")")
      ) +
      khadijat_theme()
  })

  output$met_delta_plot <- renderPlot({ met_delta_plot_obj() })

  met_heatmap_plot_obj <- reactive({
    md <- met_data()
    df <- met_df_use()
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))
    rep_handling <- input$met_rep_handling %||% "mean"

    ggplot(df, aes(x = factor(time_h), y = strain, fill = value_use)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_viridis_c(option = "C", na.value = "#f1f1f1") +
      labs(
        title = paste0(md$title, " (heatmap)"),
        subtitle = paste0("Replicate handling: ", rep_handling),
        x = "Time (h)",
        y = "Strain",
        fill = if (identical(rep_handling, "diff")) "Difference (S-1 - S-dup)" else (md$y_axis %||% "Value")
      ) +
      khadijat_theme() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"))
  })

  output$met_heatmap_plot <- renderPlot({ met_heatmap_plot_obj() })

  ## Objective 5: Multivariate (PCA) (descriptive)
  output$met_pca_ui <- renderUI({
    src <- input$met_source %||% "met"
    rep_handling <- input$met_rep_handling %||% "mean"
    if (identical(rep_handling, "diff")) {
      return(div(class = "app-note", tags$p(tags$b("Not available:"), " PCA is disabled when replicate handling is Difference (QC).")))
    }
    if (identical(src, "red24")) {
      return(div(class = "app-note", tags$p(tags$b("Not available:"), " PCA requires multiple timepoints; Reducing sugar-24h has only one timepoint.")))
    }

    tagList(
      fluidRow(
        column(
          4,
          selectInput(
            "met_pca_transform",
            "PCA transform",
            choices = c("log1p (Recommended)" = "log1p", "sqrt" = "sqrt", "None" = "none"),
            selected = "log1p"
          )
        ),
        column(
          4,
          selectInput(
            "met_pca_features",
            "Features used",
            choices = c("All features" = "all", "Top by mean value" = "top"),
            selected = "top"
          )
        ),
        column(4, checkboxInput("met_pca_scale", "Scale features", value = TRUE))
      ),
      conditionalPanel(
        condition = "input.met_pca_features === 'top'",
        sliderInput("met_pca_top_n", "Top features (mean)", min = 5, max = 60, value = 20, step = 1)
      ),
      fluidRow(
        column(
          4,
          selectInput(
            "met_pca_color_by",
            "Color by",
            choices = c("Strain" = "strain", "Time" = "time"),
            selected = "strain"
          )
        ),
        column(4, checkboxInput("met_pca_zero_fill", "0-fill non-detects (ordination only)", value = TRUE)),
        column(4, div(class = "app-note", tags$p("Exploratory only; no inferential claims.")))
      )
    )
  })

  met_pca_fit <- reactive({
    src <- input$met_source %||% "met"
    rep_handling <- input$met_rep_handling %||% "mean"
    if (identical(rep_handling, "diff")) return(NULL)
    if (identical(src, "red24")) return(NULL)

    d <- data_all()
    raw <- if (identical(src, "aut")) d$autolytic else d$metabolites
    if (nrow(raw) == 0) return(NULL)

    df0 <- if (identical(src, "aut")) {
      raw |>
        transmute(strain, time_h, feature = analyte, tech_rep, value)
    } else {
      raw |>
        transmute(strain, time_h, feature = paste0(type_unit, " | ", analyte), tech_rep, value)
    }

    df_use <- apply_techrep_handling(df0, group_cols = c("strain", "time_h", "feature"), handling = rep_handling, value_col = "value") |>
      transmute(strain, time_h, feature, value_use)

    # Feature selection
    feat_mode <- input$met_pca_features %||% "top"
    if (identical(feat_mode, "top")) {
      top_n <- as.integer(input$met_pca_top_n %||% 20L)
      feat_tbl <- df_use |>
        group_by(feature) |>
        summarise(m = if (all(is.na(.data$value_use))) NA_real_ else mean(.data$value_use, na.rm = TRUE), .groups = "drop") |>
        mutate(m = ifelse(is.nan(.data$m) | is.infinite(.data$m), NA_real_, .data$m)) |>
        arrange(desc(.data$m)) |>
        slice_head(n = top_n)
      df_use <- df_use |>
        semi_join(feat_tbl, by = "feature")
    }

    zero_fill <- isTRUE(input$met_pca_zero_fill %||% TRUE)
    wide <- df_use |>
      tidyr::pivot_wider(names_from = feature, values_from = value_use, values_fill = if (zero_fill) 0 else NA_real_)

    mat <- as.matrix(wide[, setdiff(names(wide), c("strain", "time_h"))])
    if (!zero_fill && any(!is.finite(mat))) mat[!is.finite(mat)] <- NA_real_
    if (!zero_fill && any(is.na(mat))) {
      # With non-detects treated as NA, PCA is only possible on complete cases.
      keep_rows <- stats::complete.cases(mat)
      mat <- mat[keep_rows, , drop = FALSE]
      wide <- wide[keep_rows, , drop = FALSE]
    }

    if (nrow(mat) < 2 || ncol(mat) < 2) return(NULL)

    tr <- input$met_pca_transform %||% "log1p"
    if (identical(tr, "sqrt")) mat <- sqrt(pmax(mat, 0))
    if (identical(tr, "log1p")) mat <- log1p(pmax(mat, 0))

    mat[!is.finite(mat)] <- 0

    # Drop constant features
    sds <- apply(mat, 2, stats::sd)
    keep <- is.finite(sds) & sds > 0
    mat <- mat[, keep, drop = FALSE]
    if (nrow(mat) < 2 || ncol(mat) < 2) return(NULL)

    scale_flag <- isTRUE(input$met_pca_scale %||% TRUE)
    pca <- stats::prcomp(mat, center = TRUE, scale. = scale_flag)
    ve <- (pca$sdev^2) / sum(pca$sdev^2)

    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    scores$strain <- wide$strain
    scores$time_h <- wide$time_h

    load <- as.data.frame(pca$rotation[, 1:2, drop = FALSE])
    load$feature <- rownames(load)
    rownames(load) <- NULL

    list(
      pca = pca,
      scores = scores,
      var_explained = ve,
      loadings = load,
      transform = tr,
      scale = scale_flag,
      source = src,
      rep_handling = rep_handling,
      zero_fill = zero_fill
    )
  })

  met_pca_loadings_df <- reactive({
    fit <- met_pca_fit()
    if (is.null(fit)) return(tibble::tibble())
    fit$loadings |>
      transmute(
        feature = .data$feature,
        PC1 = .data$PC1,
        PC2 = .data$PC2,
        abs_PC1 = abs(.data$PC1),
        abs_PC2 = abs(.data$PC2),
        max_abs = pmax(.data$abs_PC1, .data$abs_PC2)
      ) |>
      arrange(desc(.data$max_abs)) |>
      slice_head(n = 30) |>
      select(-.data$max_abs)
  })

  met_pca_plot_obj <- reactive({
    fit <- met_pca_fit()
    if (is.null(fit)) {
      return(ggplot() + labs(title = "PCA not available", subtitle = "Not enough complete non-constant features after preprocessing."))
    }

    scores <- fit$scores
    ve <- fit$var_explained
    xlab <- paste0("PC1 (", scales::percent(ve[[1]] %||% 0, accuracy = 0.1), ")")
    ylab <- paste0("PC2 (", scales::percent(ve[[2]] %||% 0, accuracy = 0.1), ")")

    base_shapes <- c(16, 17, 15, 18, 3, 8, 4, 0, 1, 2)
    color_by <- input$met_pca_color_by %||% "strain"
    if (identical(color_by, "time")) {
      strain_levels <- sort(unique(scores$strain))
      shape_map <- stats::setNames(base_shapes[seq_along(strain_levels)], strain_levels)
      ggplot(scores, aes(x = PC1, y = PC2, color = time_h, shape = factor(strain, levels = strain_levels))) +
        geom_point(size = 2.4, alpha = 0.85) +
        scale_color_viridis_c(option = "C") +
        scale_shape_manual(values = shape_map) +
        labs(
          title = paste0("PCA (", fit$source, "; exploratory)"),
          subtitle = paste0("Replicate handling: ", fit$rep_handling, " | 0-fill: ", fit$zero_fill),
          x = xlab,
          y = ylab,
          shape = "Strain",
          color = "Time (h)"
        ) +
        khadijat_theme()
    } else {
      time_levels <- sort(unique(scores$time_h))
      shape_map <- stats::setNames(base_shapes[seq_along(time_levels)], as.character(time_levels))
      ggplot(scores, aes(x = PC1, y = PC2, color = strain, shape = factor(time_h, levels = time_levels))) +
        geom_point(size = 2.4, alpha = 0.85) +
        scale_color_manual(values = strain_palette(scores$strain)) +
        scale_shape_manual(values = shape_map) +
        labs(
          title = paste0("PCA (", fit$source, "; exploratory)"),
          subtitle = paste0("Replicate handling: ", fit$rep_handling, " | 0-fill: ", fit$zero_fill),
          x = xlab,
          y = ylab,
          shape = "Time (h)",
          color = "Strain"
        ) +
        khadijat_theme()
    }
  })

  output$met_pca_plot <- renderPlot({ met_pca_plot_obj() })

  output$met_pca_loadings <- renderDT({
    df <- met_pca_loadings_df()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  ## Objective 5: Yield/Efficiency summaries (Metabolites.xlsx / Yeast metabolites)
  output$met_yield_ui <- renderUI({
    d <- data_all()
    feats <- d$metabolites |>
      distinct(type_unit, analyte) |>
      arrange(type_unit, analyte)
    choices <- paste0(feats$type_unit, " | ", feats$analyte)
    tp <- sort(unique(d$metabolites$time_h))

    # Defaults: first Sugar as substrate; first Acid as product; fallback to first two.
    sub_default <- choices[[1]] %||% NULL
    prod_default <- choices[[min(2, length(choices))]] %||% NULL
    if (length(choices) > 0) {
      i_sugar <- which(str_detect(standardize_key(feats$type_unit), "sugar"))[1]
      if (!is.na(i_sugar)) sub_default <- choices[[i_sugar]]
      i_acid <- which(str_detect(standardize_key(feats$type_unit), "acid"))[1]
      if (!is.na(i_acid)) prod_default <- choices[[i_acid]]
    }

    tagList(
      fluidRow(
        column(6, selectInput("yield_substrate", "Substrate (type | analyte)", choices = choices, selected = sub_default)),
        column(6, selectInput("yield_product", "Product (type | analyte)", choices = choices, selected = prod_default))
      ),
      fluidRow(
        column(6, selectInput("yield_t1", "t1 (h)", choices = tp, selected = tp[[1]] %||% NULL)),
        column(6, selectInput("yield_t2", "t2 (h)", choices = tp, selected = tp[[length(tp)]] %||% NULL))
      ),
      div(class = "app-note", tags$p("Uses mean of detected values only (0 treated as non-detect/NA)."))
    )
  })

  met_yield_df <- reactive({
    d <- data_all()
    req(input$yield_substrate, input$yield_product, input$yield_t1, input$yield_t2)

    parts_s <- str_split_fixed(input$yield_substrate, "\\s\\|\\s", 2)
    parts_p <- str_split_fixed(input$yield_product, "\\s\\|\\s", 2)
    tu_s <- parts_s[1, 1]
    an_s <- parts_s[1, 2]
    tu_p <- parts_p[1, 1]
    an_p <- parts_p[1, 2]

    t1 <- as.integer(input$yield_t1)
    t2 <- as.integer(input$yield_t2)

    sub <- d$metabolites |>
      filter(.data$type_unit == tu_s, .data$analyte == an_s) |>
      summarise_tech_dups(value_col = "value", group_cols = c("strain", "time_h", "type_unit", "analyte")) |>
      select(strain, time_h, sub_value = mean_value)

    prod <- d$metabolites |>
      filter(.data$type_unit == tu_p, .data$analyte == an_p) |>
      summarise_tech_dups(value_col = "value", group_cols = c("strain", "time_h", "type_unit", "analyte")) |>
      select(strain, time_h, prod_value = mean_value)

    w_sub <- sub |>
      pivot_wider(names_from = time_h, values_from = sub_value)
    w_prod <- prod |>
      pivot_wider(names_from = time_h, values_from = prod_value)

    sub_t1 <- w_sub[[as.character(t1)]] %||% rep(NA_real_, nrow(w_sub))
    sub_t2 <- w_sub[[as.character(t2)]] %||% rep(NA_real_, nrow(w_sub))
    prod_t1 <- w_prod[[as.character(t1)]] %||% rep(NA_real_, nrow(w_prod))
    prod_t2 <- w_prod[[as.character(t2)]] %||% rep(NA_real_, nrow(w_prod))

    out <- tibble::tibble(strain = canonical_strain_id(union(w_sub$strain, w_prod$strain))) |>
      left_join(tibble::tibble(strain = w_sub$strain, sub_t1 = sub_t1, sub_t2 = sub_t2), by = "strain") |>
      left_join(tibble::tibble(strain = w_prod$strain, prod_t1 = prod_t1, prod_t2 = prod_t2), by = "strain") |>
      mutate(
        t1 = t1,
        t2 = t2,
        substrate = input$yield_substrate,
        product = input$yield_product,
        delta_substrate = .data$sub_t2 - .data$sub_t1,
        delta_product = .data$prod_t2 - .data$prod_t1,
        substrate_consumed = - .data$delta_substrate,
        yield_ratio = dplyr::if_else(is.finite(.data$substrate_consumed) & .data$substrate_consumed > 0 & is.finite(.data$delta_product), .data$delta_product / .data$substrate_consumed, NA_real_)
      ) |>
      arrange(.data$strain)

    out
  })

  met_yield_plot_obj <- reactive({
    df <- met_yield_df()
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))

    pal <- strain_palette(df$strain)
    d1 <- df |>
      filter(is.finite(.data$substrate_consumed) & is.finite(.data$delta_product))

    p_scatter <- ggplot(d1, aes(x = substrate_consumed, y = delta_product, color = strain, label = strain)) +
      geom_hline(yintercept = 0, linewidth = 0.4, color = "#333333") +
      geom_vline(xintercept = 0, linewidth = 0.4, color = "#333333") +
      geom_point(size = 3, alpha = 0.9) +
      geom_text(vjust = -0.8, fontface = "bold", size = 3.5, show.legend = FALSE) +
      scale_color_manual(values = pal, guide = "none") +
      labs(
        title = "Yield/Efficiency (descriptive)",
        subtitle = paste0("Delta computed as t2 - t1 (t1=", df$t1[[1]] %||% NA, " h; t2=", df$t2[[1]] %||% NA, " h)."),
        x = "Substrate consumed (-delta substrate)",
        y = "Product formed (delta product)"
      ) +
      khadijat_theme()

    d2 <- df |>
      filter(is.finite(.data$yield_ratio))
    p_ratio <- ggplot(d2, aes(x = reorder(strain, yield_ratio), y = yield_ratio, fill = strain)) +
      geom_col(width = 0.75, color = "white") +
      coord_flip() +
      scale_fill_manual(values = pal, guide = "none") +
      labs(
        title = "Yield ratio (delta product / consumed substrate)",
        x = "Strain",
        y = "Yield ratio"
      ) +
      khadijat_theme()

    p_scatter / p_ratio + plot_layout(heights = c(1.2, 1))
  })

  output$met_yield_plot <- renderPlot({ met_yield_plot_obj() })

  output$met_yield_table <- renderDT({
    df <- met_yield_df()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  met_pvals_by_time <- reactive({
    req(isTRUE(input$met_enable_pvals))
    md <- met_data()
    df <- md$raw
    if (nrow(df) == 0) return(tibble::tibble())

    method <- input$met_pval_method %||% "aov_tukey"
    out <- list()
    k <- 0L
    for (t in sort(unique(df$time_h))) {
      dat <- df |>
        dplyr::filter(.data$time_h == t) |>
        dplyr::filter(!is.na(.data$value))
      if (nrow(dat) == 0) next

      res <- if (identical(method, "kw_wilcox")) safe_kw_wilcox(dat, response_col = "value", group_col = "strain") else safe_aov_tukey(dat, response_col = "value", group_col = "strain")
      if (nrow(res$pairwise) == 0) next

      k <- k + 1L
      out[[k]] <- res$pairwise |>
        dplyr::mutate(
          context = "by_time",
          time_h = t,
          strain = NA_character_,
          global_p = res$global_p,
          p_adj_within = .data$p_adj,
          method = method
        )
    }
    dplyr::bind_rows(out)
  })

  met_pvals_by_strain <- reactive({
    req(isTRUE(input$met_enable_pvals))
    md <- met_data()
    df <- md$raw
    if (nrow(df) == 0) return(tibble::tibble())

    method <- input$met_pval_method %||% "aov_tukey"
    out <- list()
    k <- 0L
    for (s in sort(unique(df$strain))) {
      dat <- df |>
        dplyr::filter(.data$strain == s) |>
        dplyr::filter(!is.na(.data$value))
      if (nrow(dat) == 0) next

      res <- if (identical(method, "kw_wilcox")) safe_kw_wilcox(dat, response_col = "value", group_col = "time_h") else safe_aov_tukey(dat, response_col = "value", group_col = "time_h")
      if (nrow(res$pairwise) == 0) next

      k <- k + 1L
      out[[k]] <- res$pairwise |>
        dplyr::mutate(
          context = "by_strain",
          time_h = NA_integer_,
          strain = s,
          global_p = res$global_p,
          p_adj_within = .data$p_adj,
          method = method
        )
    }
    dplyr::bind_rows(out)
  })

  output$met_pval_filter_ui <- renderUI({
    md <- met_data()
    df <- md$raw
    if (nrow(df) == 0) return(NULL)
    view <- input$met_pval_view %||% "by_time"
    if (identical(view, "by_time")) {
      tp <- sort(unique(df$time_h))
      selectInput("met_pval_time", "Filter time (h)", choices = c("All" = "__all__", tp), selected = "__all__")
    } else {
      strains <- sort(unique(df$strain))
      selectInput("met_pval_strain", "Filter strain", choices = c("All" = "__all__", strains), selected = "__all__")
    }
  })

  met_pvals_filtered <- reactive({
    req(isTRUE(input$met_enable_pvals))
    view <- input$met_pval_view %||% "by_time"
    padj_overall <- input$met_padj_overall %||% "BH"

    df <- if (identical(view, "by_strain")) met_pvals_by_strain() else met_pvals_by_time()
    if (nrow(df) == 0) return(df)

    if (identical(view, "by_time") && !is.null(input$met_pval_time) && input$met_pval_time != "__all__") {
      df <- df |>
        dplyr::filter(.data$time_h == as.integer(input$met_pval_time))
    }
    if (identical(view, "by_strain") && !is.null(input$met_pval_strain) && input$met_pval_strain != "__all__") {
      df <- df |>
        dplyr::filter(.data$strain == input$met_pval_strain)
    }

    df |>
      dplyr::mutate(p_adj_overall = stats::p.adjust(.data$p_adj_within, method = padj_overall)) |>
      dplyr::arrange(.data$p_adj_overall)
  })

  output$met_pvals_table <- renderDT({
    req(isTRUE(input$met_enable_pvals))
    df <- met_pvals_filtered()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  output$met_table <- renderDT({
    df <- met_table_df()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  output$dl_met_plot <- downloadHandler(
    filename = function() paste0("metabolites_dynamics_", standardize_key(input$met_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(met_plot_obj(), file, width = 10, height = 5.5, dpi = 320)
  )

  output$dl_met_snapshot_plot <- downloadHandler(
    filename = function() paste0("metabolites_compare_", standardize_key(input$met_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(met_snapshot_plot_obj(), file, width = 10, height = 5.5, dpi = 320)
  )

  output$dl_met_heatmap_plot <- downloadHandler(
    filename = function() paste0("metabolites_heatmap_", standardize_key(input$met_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(met_heatmap_plot_obj(), file, width = 10, height = 5.8, dpi = 320)
  )

  output$dl_met_pca_plot <- downloadHandler(
    filename = function() paste0("met_pca_", input$met_source %||% "met", "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(met_pca_plot_obj(), file, width = 10, height = 6.5, dpi = 320)
  )

  output$dl_met_pca_loadings <- downloadHandler(
    filename = function() paste0("met_pca_loadings_", input$met_source %||% "met", "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) readr::write_csv(met_pca_loadings_df(), file, na = "")
  )

  output$dl_met_yield_plot <- downloadHandler(
    filename = function() paste0("met_yield_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(met_yield_plot_obj(), file, width = 10, height = 8, dpi = 320)
  )

  output$dl_met_yield_table <- downloadHandler(
    filename = function() paste0("met_yield_table_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) writexl::write_xlsx(list(table = met_yield_df()), path = file)
  )

  output$dl_met_pvals <- downloadHandler(
    filename = function() paste0("metabolites_pvalues_", standardize_key(input$met_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) readr::write_csv(met_pvals_filtered(), file, na = "")
  )

  output$dl_met_table <- downloadHandler(
    filename = function() paste0("metabolites_summary_", standardize_key(input$met_analyte), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) writexl::write_xlsx(list(summary = met_table_df()), path = file)
  )

  ## Objective 6: Best yeast (transparent, exploratory multi-criteria ranking)
  rank_presets <- list(
    balanced = c(growth = 0.35, tolerance = 0.25, enzymes = 0.20, volatiles = 0.20),
    equal = c(growth = 0.25, tolerance = 0.25, enzymes = 0.25, volatiles = 0.25),
    growth = c(growth = 0.55, tolerance = 0.15, enzymes = 0.15, volatiles = 0.15),
    tolerance = c(growth = 0.15, tolerance = 0.55, enzymes = 0.15, volatiles = 0.15),
    enzymes = c(growth = 0.15, tolerance = 0.15, enzymes = 0.55, volatiles = 0.15),
    volatiles = c(growth = 0.15, tolerance = 0.15, enzymes = 0.15, volatiles = 0.55)
  )

  observeEvent(input$rank_preset, {
    p <- input$rank_preset %||% "custom"
    if (!p %in% names(rank_presets)) return()
    w <- rank_presets[[p]]
    updateSliderInput(session, "w_growth", value = as.numeric(w[["growth"]]))
    updateSliderInput(session, "w_tolerance", value = as.numeric(w[["tolerance"]]))
    updateSliderInput(session, "w_enz", value = as.numeric(w[["enzymes"]]))
    updateSliderInput(session, "w_vol", value = as.numeric(w[["volatiles"]]))
  }, ignoreInit = TRUE)

  normalize_weights <- function(w, metric_avail) {
    w <- w[intersect(names(w), names(metric_avail))]
    w[is.na(w)] <- 0
    w[!metric_avail[names(w)]] <- 0
    if (sum(w) == 0) w <- stats::setNames(c(1, 0, 0, 0), names(w))
    w / sum(w)
  }

  scale01 <- function(x) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    rng <- range(x, na.rm = TRUE)
    if (rng[2] - rng[1] < .Machine$double.eps) return(rep(0.5, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }

  rank_metrics <- reactive({
    d <- data_all()
    gv <- gc_tidy()

    # Growth: microbial growth log CFU/mL (detected values only).
    g <- d$metabolites |>
      filter(standardize_key(type_unit) == standardize_key("microbial growth Log (CFU/mL)"))
    g_ts <- summarise_tech_dups(g, value_col = "value", group_cols = c("strain", "time_h", "analyte")) |>
      select(strain, time_h, mean_value)
    g_max <- g_ts |>
      group_by(strain) |>
      summarise(growth_max = max(mean_value, na.rm = TRUE), .groups = "drop")
    g_max$growth_max[is.infinite(g_max$growth_max)] <- NA_real_
    g_end <- g_ts |>
      group_by(strain) |>
      slice_max(.data$time_h, n = 1, with_ties = FALSE) |>
      transmute(strain, growth_end = .data$mean_value)

    # Survival: autolytic survival log CFU/mL (detected values only).
    a <- d$autolytic |>
      filter(standardize_key(analyte) == standardize_key("microbial-survival- Log (CFU/mL)"))
    a_ts <- summarise_tech_dups(a, value_col = "value", group_cols = c("strain", "time_h", "analyte")) |>
      select(strain, time_h, mean_value)
    surv_end <- a_ts |>
      group_by(strain) |>
      slice_max(.data$time_h, n = 1, with_ties = FALSE) |>
      transmute(strain, survival_end = .data$mean_value)

    # Tolerance (broth): mean + min across conditions (0 = no growth is real 0).
    tb <- d$tol_broth
    tb_sum <- summarise_tech_dups(tb, value_col = "value", group_cols = c("strain", "condition_key"))
    tb_agg <- tb_sum |>
      group_by(strain) |>
      summarise(
        tol_broth_mean = mean(.data$mean_value, na.rm = TRUE),
        tol_broth_min = min(.data$mean_value, na.rm = TRUE),
        .groups = "drop"
      )
    tb_agg$tol_broth_mean[is.infinite(tb_agg$tol_broth_mean)] <- NA_real_
    tb_agg$tol_broth_min[is.infinite(tb_agg$tol_broth_min)] <- NA_real_

    # Tolerance (agar): mean score across conditions (ordinal; descriptive).
    ta <- d$tol_agar
    ta_agg <- ta |>
      group_by(strain) |>
      summarise(
        tol_agar_mean = if (all(is.na(.data$agar_score))) NA_real_ else mean(.data$agar_score, na.rm = TRUE),
        tol_agar_min = if (all(is.na(.data$agar_score))) NA_real_ else min(.data$agar_score, na.rm = TRUE),
        .groups = "drop"
      )

    # Enzymes: mean across measures for production-rate category only (detected values only).
    ez <- d$enzymes |>
      filter(.data$category == "Enzyme production rate")
    ez_sum <- summarise_tech_dups(ez, value_col = "value", group_cols = c("strain", "measure_key")) |>
      group_by(strain) |>
      summarise(
        enz_mean = mean(.data$mean_value, na.rm = TRUE),
        enz_n_measures = sum(!is.na(.data$mean_value)),
        .groups = "drop"
      )
    ez_sum$enz_mean[is.infinite(ez_sum$enz_mean)] <- NA_real_

    # Volatiles: total RA + richness over time (exploratory; depends on duplicate policy).
    if (nrow(gv) == 0) {
      vol_agg <- tibble::tibble(
        strain = character(),
        vol_total_mean = numeric(),
        vol_total_end = numeric(),
        vol_richness_mean = numeric()
      )
    } else {
      vol_ts <- gv |>
        group_by(strain, time_h) |>
        summarise(
          vol_total = sum(.data$ra_ug_L, na.rm = TRUE),
          vol_richness = sum(.data$detected, na.rm = TRUE),
          .groups = "drop"
        )
      vol_agg <- vol_ts |>
        group_by(strain) |>
        summarise(
          vol_total_mean = mean(.data$vol_total, na.rm = TRUE),
          vol_richness_mean = mean(.data$vol_richness, na.rm = TRUE),
          .groups = "drop"
        ) |>
        left_join(
          vol_ts |>
            group_by(strain) |>
            slice_max(.data$time_h, n = 1, with_ties = FALSE) |>
            transmute(strain, vol_total_end = .data$vol_total),
          by = "strain"
        )
      vol_agg$vol_total_mean[is.infinite(vol_agg$vol_total_mean)] <- NA_real_
      vol_agg$vol_richness_mean[is.infinite(vol_agg$vol_richness_mean)] <- NA_real_
      vol_agg$vol_total_end[is.infinite(vol_agg$vol_total_end)] <- NA_real_
    }

    all_strains <- sort(unique(c(
      d$metabolites$strain,
      d$autolytic$strain,
      d$tol_broth$strain,
      d$tol_agar$strain,
      d$enzymes$strain,
      gv$strain
    )))

    out <- tibble::tibble(strain = canonical_strain_id(all_strains)) |>
      left_join(g_max, by = "strain") |>
      left_join(g_end, by = "strain") |>
      left_join(surv_end, by = "strain") |>
      left_join(tb_agg, by = "strain") |>
      left_join(ta_agg, by = "strain") |>
      left_join(ez_sum, by = "strain") |>
      left_join(vol_agg, by = "strain") |>
      mutate(
        gc_duplicate_policy = input$gc_dup_policy %||% NA_character_
      )

    out
  })

  rank_metric_avail <- reactive({
    m <- rank_metrics()
    c(
      growth = !all(is.na(m$growth_max)),
      tolerance = !all(is.na(m$tol_broth_mean)),
      enzymes = !all(is.na(m$enz_mean)),
      volatiles = !all(is.na(m$vol_total_mean))
    )
  })

  rank_weights_used <- reactive({
    avail <- rank_metric_avail()
    w_in <- c(growth = input$w_growth, tolerance = input$w_tolerance, enzymes = input$w_enz, volatiles = input$w_vol)
    normalize_weights(w_in, avail)
  })

  rank_scaled <- reactive({
    m <- rank_metrics()
    m |>
      mutate(
        s_growth = scale01(.data$growth_max),
        s_tol = scale01(.data$tol_broth_mean),
        s_enz = scale01(.data$enz_mean),
        s_vol = scale01(.data$vol_total_mean),
        s_growth = ifelse(is.na(.data$s_growth), 0, .data$s_growth),
        s_tol = ifelse(is.na(.data$s_tol), 0, .data$s_tol),
        s_enz = ifelse(is.na(.data$s_enz), 0, .data$s_enz),
        s_vol = ifelse(is.na(.data$s_vol), 0, .data$s_vol)
      )
  })

  rank_df <- reactive({
    df <- rank_scaled()
    avail <- rank_metric_avail()
    w <- rank_weights_used()

    df <- df |>
      mutate(
        c_growth = w["growth"] * .data$s_growth,
        c_tol = w["tolerance"] * .data$s_tol,
        c_enz = w["enzymes"] * .data$s_enz,
        c_vol = w["volatiles"] * .data$s_vol,
        score = .data$c_growth + .data$c_tol + .data$c_enz + .data$c_vol,
        w_growth = w["growth"],
        w_tolerance = w["tolerance"],
        w_enzymes = w["enzymes"],
        w_volatiles = w["volatiles"],
        note = paste0(
          "Exploratory ranking (technical duplicates; no biological replication). Metrics used: ",
          paste(names(avail)[avail], collapse = ", "),
          ". No co-culture/synergy in v1."
        )
      ) |>
      arrange(desc(.data$score)) |>
      mutate(rank = dplyr::dense_rank(desc(.data$score))) |>
      select(rank, everything())

    df
  })

  rank_table_df <- reactive({
    df <- rank_df()
    df |>
      select(
        rank, strain, score,
        growth_max, tol_broth_mean, enz_mean, vol_total_mean,
        w_growth, w_tolerance, w_enzymes, w_volatiles,
        gc_duplicate_policy
      )
  })

  rank_components_df <- reactive({
    df <- rank_df()
    df |>
      select(
        rank, strain, score,
        growth_max, growth_end, survival_end,
        tol_broth_mean, tol_broth_min, tol_agar_mean, tol_agar_min,
        enz_mean, enz_n_measures,
        vol_total_mean, vol_total_end, vol_richness_mean,
        s_growth, s_tol, s_enz, s_vol,
        c_growth, c_tol, c_enz, c_vol,
        w_growth, w_tolerance, w_enzymes, w_volatiles,
        gc_duplicate_policy
      )
  })

  rank_weights_tbl <- reactive({
    avail <- rank_metric_avail()
    w_in <- c(growth = input$w_growth, tolerance = input$w_tolerance, enzymes = input$w_enz, volatiles = input$w_vol)
    w_used <- rank_weights_used()
    tibble::tibble(
      metric = names(w_in),
      available = as.logical(avail[names(w_in)]),
      weight_input = as.numeric(w_in),
      weight_used = as.numeric(w_used[names(w_in)])
    )
  })

  output$rank_weights <- renderDT({
    df <- rank_weights_tbl()
    datatable(
      df,
      rownames = FALSE,
      options = list(dom = "t", paging = FALSE, searching = FALSE, info = FALSE, ordering = FALSE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  rank_plot_obj <- reactive({
    df <- rank_df()
    req(nrow(df) > 0)
    ggplot(df, aes(x = reorder(strain, score), y = score, fill = strain)) +
      geom_col(width = 0.75, color = "white") +
      coord_flip() +
      scale_fill_manual(values = strain_palette(df$strain), guide = "none") +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
      labs(title = "Composite score by strain (exploratory)", x = "Strain", y = "Score (0-1)") +
      khadijat_theme()
  })

  output$rank_plot <- renderPlot({ rank_plot_obj() })

  output$rank_table <- renderDT({
    df <- rank_table_df()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  rank_contrib_plot_obj <- reactive({
    df <- rank_df()
    req(nrow(df) > 0)
    long <- df |>
      select(strain, c_growth, c_tol, c_enz, c_vol) |>
      tidyr::pivot_longer(cols = -strain, names_to = "component", values_to = "contrib") |>
      mutate(
        component = recode(
          component,
          c_growth = "Growth",
          c_tol = "Tolerance",
          c_enz = "Enzymes",
          c_vol = "Volatiles"
        )
      )
    ggplot(long, aes(x = reorder(strain, contrib, sum), y = contrib, fill = component)) +
      geom_col(width = 0.75, color = "white", linewidth = 0.2) +
      coord_flip() +
      scale_fill_viridis_d(option = "C") +
      labs(title = "Score contributions by metric (exploratory)", x = "Strain", y = "Weighted contribution") +
      khadijat_theme()
  })

  output$rank_contrib_plot <- renderPlot({ rank_contrib_plot_obj() })

  output$rank_components_table <- renderDT({
    df <- rank_components_df()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  rank_sensitivity_df <- reactive({
    base <- rank_scaled()
    avail <- rank_metric_avail()

    presets <- c(names(rank_presets), "current")
    out <- list()
    k <- 0L
    for (p in presets) {
      w_in <- if (identical(p, "current")) {
        c(growth = input$w_growth, tolerance = input$w_tolerance, enzymes = input$w_enz, volatiles = input$w_vol)
      } else {
        rank_presets[[p]]
      }
      w <- normalize_weights(w_in, avail)
      tmp <- base |>
        mutate(
          score = w["growth"] * .data$s_growth + w["tolerance"] * .data$s_tol + w["enzymes"] * .data$s_enz + w["volatiles"] * .data$s_vol,
          preset = p,
          w_growth = w["growth"],
          w_tolerance = w["tolerance"],
          w_enzymes = w["enzymes"],
          w_volatiles = w["volatiles"]
        ) |>
        arrange(desc(.data$score)) |>
        mutate(rank = dplyr::dense_rank(desc(.data$score))) |>
        select(preset, rank, strain, score, w_growth, w_tolerance, w_enzymes, w_volatiles)
      k <- k + 1L
      out[[k]] <- tmp
    }
    dplyr::bind_rows(out) |>
      mutate(preset = factor(.data$preset, levels = presets))
  })

  rank_sensitivity_plot_obj <- reactive({
    df <- rank_sensitivity_df()
    if (nrow(df) == 0) return(ggplot() + labs(title = "Not available in current dataset"))
    ggplot(df, aes(x = preset, y = rank, group = strain, color = strain)) +
      geom_line(linewidth = 0.8, alpha = 0.85) +
      geom_point(size = 2.2, alpha = 0.9) +
      scale_color_manual(values = strain_palette(df$strain)) +
      scale_y_reverse(breaks = sort(unique(df$rank))) +
      labs(
        title = "Rank sensitivity across weight presets (exploratory)",
        x = "Preset",
        y = "Rank (1 = best)"
      ) +
      khadijat_theme() +
      theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"))
  })

  output$rank_sensitivity_plot <- renderPlot({ rank_sensitivity_plot_obj() })

  output$rank_sensitivity_table <- renderDT({
    df <- rank_sensitivity_df()
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 12, scrollX = TRUE, columnDefs = dt_coldefs_numeric(df, digits = 4L))
    )
  })

  output$rank_trade_ui <- renderUI({
    df <- rank_metrics()
    if (nrow(df) == 0) return(NULL)
    choices <- c(
      "Growth max (CFU)" = "growth_max",
      "Growth endpoint (CFU)" = "growth_end",
      "Survival endpoint (CFU)" = "survival_end",
      "Tolerance broth mean" = "tol_broth_mean",
      "Tolerance broth min" = "tol_broth_min",
      "Tolerance agar mean" = "tol_agar_mean",
      "Enzymes mean" = "enz_mean",
      "Volatiles total mean" = "vol_total_mean",
      "Volatiles total endpoint" = "vol_total_end",
      "Volatiles richness mean" = "vol_richness_mean"
    )
    tagList(
      fluidRow(
        column(5, selectInput("rank_trade_x", "X metric", choices = choices, selected = "growth_max")),
        column(5, selectInput("rank_trade_y", "Y metric", choices = choices, selected = "vol_total_mean")),
        column(2, checkboxInput("rank_trade_scaled", "Scale 0-1", value = TRUE))
      ),
      div(class = "app-note", tags$p("Trade-off plots are descriptive only and use the current dataset (technical duplicates)."))
    )
  })

  rank_trade_plot_obj <- reactive({
    df <- rank_metrics()
    req(nrow(df) > 0)
    req(input$rank_trade_x, input$rank_trade_y)
    xcol <- input$rank_trade_x
    ycol <- input$rank_trade_y

    plot_df <- df |>
      transmute(
        strain = .data$strain,
        x = suppressWarnings(as.numeric(.data[[xcol]])),
        y = suppressWarnings(as.numeric(.data[[ycol]]))
      )

    if (isTRUE(input$rank_trade_scaled %||% TRUE)) {
      plot_df <- plot_df |>
        mutate(
          x = scale01(.data$x),
          y = scale01(.data$y)
        )
    }

    ggplot(plot_df, aes(x = x, y = y, color = strain)) +
      geom_point(size = 3, alpha = 0.9) +
      geom_text(aes(label = strain), nudge_y = 0.02, size = 3.5, show.legend = FALSE) +
      scale_color_manual(values = strain_palette(plot_df$strain), guide = "none") +
      labs(
        title = "Metric trade-offs (descriptive)",
        x = xcol,
        y = ycol
      ) +
      khadijat_theme()
  })

  output$rank_trade_plot <- renderPlot({ rank_trade_plot_obj() })

  output$dl_rank_plot <- downloadHandler(
    filename = function() paste0("best_yeast_ranking_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(rank_plot_obj(), file, width = 10, height = 6, dpi = 320)
  )

  output$dl_rank_contrib_plot <- downloadHandler(
    filename = function() paste0("best_yeast_contributions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(rank_contrib_plot_obj(), file, width = 10, height = 6, dpi = 320)
  )

  output$dl_rank_sensitivity_plot <- downloadHandler(
    filename = function() paste0("best_yeast_sensitivity_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(rank_sensitivity_plot_obj(), file, width = 11, height = 6, dpi = 320)
  )

  output$dl_rank_table <- downloadHandler(
    filename = function() paste0("best_yeast_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) {
      writexl::write_xlsx(
        list(
          ranking = rank_table_df(),
          weights = rank_weights_tbl(),
          components = rank_components_df(),
          sensitivity = rank_sensitivity_df()
        ),
        path = file
      )
    }
  )

  ## QC tab: S-1 vs S-dup scatter (technical duplicates)
  output$qc_analyte_ui <- renderUI({
    d <- data_all()
    if (input$qc_dataset == "met") {
      choices <- sort(unique(d$metabolites$analyte))
    } else if (input$qc_dataset == "aut") {
      choices <- sort(unique(d$autolytic$analyte))
    } else if (input$qc_dataset == "tol") {
      choices <- sort(unique(d$tol_broth$concentration))
    } else {
      choices <- sort(unique(d$enzymes$measure))
    }
    selectInput("qc_analyte", "Analyte/measure", choices = choices, selected = choices[[1]] %||% NULL)
  })

  qc_plot_obj <- reactive({
    d <- data_all()
    req(input$qc_analyte)
    if (input$qc_dataset == "met") {
      df <- d$metabolites |>
        filter(analyte == input$qc_analyte) |>
        select(strain, time_h, tech_rep, value_raw) |>
        mutate(key = paste(strain, time_h, sep = "_"))
    } else if (input$qc_dataset == "aut") {
      df <- d$autolytic |>
        filter(analyte == input$qc_analyte) |>
        select(strain, time_h, tech_rep, value_raw) |>
        mutate(key = paste(strain, time_h, sep = "_"))
    } else if (input$qc_dataset == "tol") {
      df <- d$tol_broth |>
        filter(concentration == input$qc_analyte) |>
        select(strain, analyte, concentration, tech_rep, value_raw) |>
        mutate(key = paste(strain, analyte, concentration, sep = "_"))
    } else {
      df <- d$enzymes |>
        filter(measure == input$qc_analyte) |>
        select(strain, category, measure, tech_rep, value_raw) |>
        mutate(key = paste(strain, category, measure, sep = "_"))
    }

    wide <- df |>
      select(key, tech_rep, value_raw) |>
      group_by(.data$key, .data$tech_rep) |>
      summarise(value_raw = if (all(is.na(.data$value_raw))) NA_real_ else mean(.data$value_raw, na.rm = TRUE), .groups = "drop") |>
      pivot_wider(names_from = tech_rep, values_from = value_raw)

    req(all(c("S-1", "S-dup") %in% names(wide)))

    dat <- wide |>
      dplyr::transmute(x = .data[["S-1"]], y = .data[["S-dup"]]) |>
      dplyr::filter(is.finite(.data$x) & is.finite(.data$y))

    n <- nrow(dat)
    r <- if (n >= 2) suppressWarnings(stats::cor(dat$x, dat$y, method = "pearson")) else NA_real_
    fit <- if (n >= 2 && stats::sd(dat$x) > 0) stats::lm(y ~ x, data = dat) else NULL
    a <- if (!is.null(fit)) unname(stats::coef(fit)[[1]]) else NA_real_
    b <- if (!is.null(fit)) unname(stats::coef(fit)[[2]]) else NA_real_

    fmt4 <- function(x) {
      if (is.na(x) || !is.finite(x)) return("NA")
      if (abs(x) < 0.5 * 10^(-4)) x <- 0
      sprintf("%.4f", x)
    }
    lbl <- paste0(
      "n = ", n,
      "\nPearson r = ", fmt4(r),
      "\ny = ", fmt4(a), " + ", fmt4(b), " x"
    )

    if (n >= 1) {
      x_rng <- range(dat$x, na.rm = TRUE)
      y_rng <- range(dat$y, na.rm = TRUE)
      x_pos <- x_rng[[1]] + 0.05 * diff(x_rng)
      y_pos <- y_rng[[2]] - 0.05 * diff(y_rng)
    } else {
      x_pos <- 0
      y_pos <- 0
    }

    ggplot(dat, aes(x = x, y = y)) +
      geom_point(alpha = 0.7, size = 2, color = "#0b5cab") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.85, color = "#d95f02") +
      annotate(
        "label",
        x = x_pos,
        y = y_pos,
        label = lbl,
        hjust = 0,
        vjust = 1,
        size = 3.6,
        label.size = 0.25,
        alpha = 0.92,
        fill = "white"
      ) +
      labs(
        title = "Technical duplicate agreement (S-1 vs S-dup)",
        subtitle = "Points on the dashed line indicate perfect agreement. Large deviations flag measurement/QC issues.",
        x = "S-1 (raw)",
        y = "S-dup (raw)"
      ) +
      khadijat_theme()
  })

  output$qc_plot <- renderPlot({ qc_plot_obj() })

  output$dl_qc_plot <- downloadHandler(
    filename = function() paste0("qc_s1_vs_sdup_", input$qc_dataset, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_plot(qc_plot_obj(), file, width = 7.5, height = 6, dpi = 320)
  )
}

shinyApp(ui, server)
