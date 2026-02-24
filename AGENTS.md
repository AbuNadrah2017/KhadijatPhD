# Analysis Plan and Peer-Review Checklist (KhadijatPhD)

This repo's goal is to build an **R Shiny** app that **loads the provided datasets internally** (no user upload required) and exposes **peer-review-grade analyses** as interactive tabs. Each Objective below should become a Results tab containing multiple analyses, figures, and exportable tables.

## Objectives (from `Objectives/Objectives- data analysis.docx`)
Aim: Investigate yeast growth dynamics, enzymatic activities, and metabolic interactions in cocoa pulp simulation media, focusing on volatile compound production and microbial synergy during fermentation.

1. Determine the growth characteristics of each yeast species/strain through mono- and co-culture fermentation experiments.
2. Evaluate the physiological properties of each isolated yeast and their combinations regarding tolerance to pH, temperature, acetic acid, citric acid, maltose, glucose, fructose, total reducing sugar by DNS, and ethanol.
3. Assess the enzymatic activities of each isolated yeast in cocoa pulp medium (invertase, pectinase, polygalacturonase, cellulase).
4. Identify and characterize the volatile compounds produced by each yeast species during mono- and co-culture fermentation experiments.
5. Analyse the metabolic interactions in co-cultures to identify synergistic or antagonistic effects on metabolite production (e.g., ethanol, organic acids, flavour compounds).
6. Determine the best performing yeast species.

Scope note (current datasets): there are no co-culture samples/identifiers in the reviewed files, so analyses will be **mono-culture only**. Co-culture comparison and synergy/antagonism claims should be disabled/deferred until co-culture data is added.
Hard constraint (v1): do not include any co-culture UI, co-culture comparisons, or synergy/antagonism scoring in the app until co-culture samples exist and are explicitly encoded in the data/metadata.
Hard constraint (data integrity): all figures/tables/statistics must be computed from the loaded datasets. No hardcoded results, no simulated/fake data, and no placeholder outputs. If an assay/condition is missing, show "Not available in current dataset" and do not generate that analysis.

## Current Shiny App (v1): Implemented Analyses
The Shiny app is implemented in `R/app.R` and provides an Overview tab, Objective 1-6 tabs, and a QC tab.

- Data sources are loaded internally from:
- `R/Metabolites.xlsx`
- `R/Yeast tolerance test.xlsx`
- `R/GC_VolatilesYeast.xlsx`
- `R/exampleMapData.xlsx` (mapping template; incomplete by design)
- Each Objective tab includes a visible "Data sources (file / sheet)" note so users can trace which spreadsheet + sheet(s) are used.
- Each Objective tab includes concise on-page explanations (what the analysis shows, what data are used, and key limitations). Any p-values are explicitly labeled as technical-duplicate / exploratory (no biological replication).
- Figures are styled for publication using `khadijat_theme()` (bold axis labels/ticks/legend; centered bold titles).
- All numeric tables are formatted to 4 decimal places (p-values show `<0.0001` when smaller than the display threshold).
- Downloads are available per tab (plots as PNG; tables as XLSX; p-values as CSV where enabled), plus session info and analysis parameters.
- Co-culture comparisons and synergy/antagonism scoring are explicitly deferred (no co-culture IDs in the current dataset).
- Any p-values in v1 are opt-in and explicitly labeled as "technical duplicates / exploratory" (no biological replication).

## Audit Log (Implementation + Tests)
Date: 2026-02-08

Test command:
- `Rscript scripts/run_tests.R`

Result:
- Failures: 0
- Errors: 0
- Warnings: present (non-fatal; 58 in the latest test run)

Findings (addressed in code):
- Objective 2 dose-response: fixed a many-to-many join that could duplicate rows when parsing numeric concentrations (mapping is now built from unique concentration strings).
- Objective 2 tolerance index: added server-side default weights when slider inputs are missing (robustness in non-interactive/test contexts).
- Objective 4 data-source note: clarified that `QC-All` is present in the GC workbook but is not used in v1 (no background subtraction).
- Objective 5 downloads: added download handlers for PCA (plot + loadings) and yield/efficiency (plot + table).
- Tests: expanded coverage to include Objective 1 growth-model fitting, Objective 2 tolerance index, and Objective 5 PCA/yield plus source switching.

Warnings observed (do not change computed results, but should be kept in mind for maintenance):
- Packages built under a different R version (environment-level warning only).
- tidyselect deprecation warnings (future-proofing task).
- ggplot warnings about removed rows (expected when non-detects are treated as `NA` and therefore dropped from geoms).
- nls/growth-model warnings (step factor/singular gradient) can appear for some model fits; failed fits are flagged as non-converged and excluded from "best model" display.

Detailed audit checks (peer-review + integrity):
- Data provenance: all analytic datasets are loaded from the provided Excel workbooks via `readxl::read_excel()` wrappers (`R/load_data.R`, `R/gc_parse.R`). No built-in/demo datasets (e.g., `iris`, `mtcars`) are used.
- No simulated/fake results: codebase scan found no random data generation (`rnorm`, `runif`, `sample`, `set.seed`) and no placeholder result tables. Any fixed constants are limited to UI styling (themes/palettes) and transparent scoring presets (weights), not computed outcomes.
- No writing back to inputs: the app does not modify the source spreadsheets. All file output is via explicit user-triggered downloads (PNG/CSV/XLSX/TXT/JSON).
- Co-culture exclusion enforced: no co-culture analyses or synergy/antagonism scoring are implemented in v1; the UI and ranking logic explicitly defer these until co-culture IDs exist.
- Replicate structure enforced: `S-1`/`S-dup` are treated as technical duplicates (repeat measurements of the same sample/run), not independent fermentations. Any p-values are opt-in and explicitly labeled as technical-duplicate/exploratory only.
- Non-detect policy enforced: for concentration-like sheets, `0` is treated as non-detect (converted to `NA`) and is not imputed; detection rates are reported. Exception: tolerance broth/agar retain `0` as real "no growth".
- Download completeness: all `downloadButton()` entries in `R/app.R` have matching `downloadHandler()` implementations (37/37).

Q1 peer-review readiness statement (honest constraint):
- The app is engineered to be peer-review-safe as an exploratory/descriptive analysis environment with transparent data handling and reproducible exports. However, Q1 papers that make biological between-strain claims typically require independent fermentation replicates encoded as `run_id`/`batch_id` (not present in the current dataset). Until such replication exists, inferential outputs (p-values, CIs, "best strain" claims) must be reported as exploratory/technical-duplicate only and interpreted accordingly.

## Objective Coverage (current dataset)
The app includes an "Objective coverage" panel (Overview tab) that lists what is supported vs not available in the current dataset.

- Objective 1 (growth): supported for mono-cultures via CFU, OD, and survival time series. Co-culture comparisons are not available.
- Objective 2 (tolerance): supported for the conditions present in the tolerance sheets. Any objective items not present as conditions (e.g., maltose/citric acid, if absent) must be shown as "Not available in current dataset".
- Objective 3 (enzymes): supported for the assays present in `Enzyme production` (e.g., PG/PME/PNL/beta-glucosidase/protease). If invertase/cellulase are not measured directly, show "Not available in current dataset" (do not infer).
- Objective 4 (volatiles): supported for mono-culture volatile profiles over time. Co-culture volatile comparisons are not available. Mapping coverage to the example compound mapping is reported; full sensory/class summaries require a complete mapping table (deferred).
- Objective 5 (metabolic dynamics): supported for mono-culture metabolic dynamics using `R/Metabolites.xlsx / Yeast metabolites`, plus additional physiology panels from `R/Yeast tolerance test.xlsx / Autolytic-full` and `R/Metabolites.xlsx / Reducing sugar-24h`. Co-culture synergy/antagonism is not available.
- Objective 6 (best yeast): supported as a transparent, multi-criteria ranking using only available endpoints, with sensitivity analysis.

## Data Inventory (current)

### `R/Metabolites.xlsx`
Sheet: `Yeast metabolites` (time series; replicates `S-1`, `S-dup`)

- Columns: `Sample-Name`, `Type-unit`, `Analyte-Acid`, `sample-h`, `S-1`, `S-dup`
- Samples seen (normalize any legacy `6B` to `Y6B`): `SC`, `Y1`, `Y1A`, `Y3A`, `Y4`, `Y5`, `Y6B`
- Timepoints in `sample-h`: `0, 12, 24, 36, 48, 72, 96` hours
- Key analytes (as `Analyte-Acid` + `Type-unit`):
- Sugars (g/L): `Glucose`, `Fructose`, `Sucrose`
- Organic acids (g/L): `citric`, `lactic`, `malic`, `succinic`
- `pH`
- Microbial growth: `microbial growth Log (CFU/mL)`

Sheet: `Reducing sugar-24h` (replicates `S-1`, `S-dup`)

- Columns: `Analyte`, `sample`, `S-1`, `S-dup`
- Analytes: `Extracellular_g_per_L`, `Intracellular_g_per_L_`, `Total_g_per_L`

### `R/Yeast tolerance test.xlsx`
Sheet: `tolerance test-agar` (ordinal plate score)

- Columns: `sample`, `analyte`, `concentration`, `agar-score`
- `agar-score` values observed: `0, 2, 3, 5`
- Stress/analyte groups observed include: `osmotolerant`, `Ethanol`, `Thermotolerance`, `Acetic acid`, `Acid`, `Combine stress`, plus enzyme-screen style entries.

Sheet: `tolerance test-broth` (quantitative; replicates `S-1`, `S-dup`)

- Columns: `sample`, `analyte`, `concentration`, `S-1`, `S-dup`
- Stress/analyte groups observed include: `Sugar`, `Ethanol`, `Thermotolerance`, `Acid/pH`, `Acetic Acid` (note inconsistent spelling/casing)

Sheet: `Autolytic-full` (time series; replicates `S-1`, `S-dup`)

- Columns: `Sample-Name`, `Analyte-Acid`, `sample-h`, `S-1`, `S-dup`
- Timepoints in `sample-h`: `0, 24, 48, 72, 96` hours
- Analytes observed include: `OD`, `microbial-survival- Log (CFU/mL)`, `FAN-N(mg/L)`, `Protein concentration (mg/L, BSA equivalents)`, `Reducing sugar-DNS, mg/L only`

Sheet: `Enzyme production` (replicates `S-1`, `S-dup`)

- Columns: `sample`, `analyte-Log(CFU/g )`, `concentration`, `S-1`, `S-dup`
- Categories in `analyte-Log(CFU/g )`: `Enzyme production`, `Enzyme degradation`, `Enzyme production rate`
- Measures in `concentration` include: `PG-*`, `PME-*`, `PNL-*`, `beta glucosidase-*`, `protease-*`, and `protein (mg/mL CE)`

### `R/GC_VolatilesYeast.xlsx`
Sheets: `QC-All`, `SC-Volatiles`, `Y1-Volatiles`, `Y1A-Volatiles`, `Y3A-Volatiles`, `Y4-Volatiles`, `Y5-volatiles`, `Y6B-Volatiles`

- Structure: each sheet is multiple timepoint blocks laid out side-by-side.
- For each timepoint, there is a column for `RA (ug/L)` and a column for compound `Names` (or `Common (sensory) name` in some places).
- Timepoints present across strains: `0, 12, 24, 36, 48, 72, 96` hours.
- `QC-All` appears to be a blank/matrix/QC profile (background compounds).
- Known header issues to handle during parsing:
- In `Y5-volatiles`, the 0h RA header is `Y5-RA_ug_per_L` (time is implied by adjacent `Y5-0-Names`).
- In `Y6B-Volatiles`, the 72h/96h RA headers contain extra text like `A1:B31` / `A1:B25` that can break naive time extraction.

### `R/exampleMapData.xlsx`
Provided as a sample mapping of compound class and common (sensory) names.

- Sheets:
- `compound_map` with columns: `Compound class`, `Common (sensory) name`, `Sensory note / comment`
- `examplegrouping` showing an example grouped output with columns: `Time`, `Volatiles`, `Compound class`, `Common (sensory) name`, `relative concentration`, `Unit`
- Important: `R/exampleMapData.xlsx` is an example only and does not contain a complete mapping for all compounds. A detailed, project-specific mapping table must be generated and maintained (covering all compounds and any sensor groupings you want), and the Shiny app should use that complete mapping for class/sensory summaries and reporting.
- Minimum recommended columns for the full mapping table:
- `raw_compound_name` (as it appears in `R/GC_VolatilesYeast.xlsx`)
- `canonical_compound_name` (after standardization; used as the join key)
- `compound_class` (e.g., Alcohol, Ester, Acid)
- `common_sensory_name` (optional, if different from canonical)
- `sensory_note_comment` (optional)
- `sensor_group` (optional; only if you define sensors/groupings)
- The Shiny app should treat `compound_map` as the template for a mapping dictionary from standardized compound names to compound class and sensory notes. (If you later want sensor hardware mappings, add columns explicitly; they are not present in the current file.)

## Metadata Needed (for peer review)
These items are not fully encoded in the current spreadsheets but are typically required to defend the analysis in peer review.

- Strain metadata: a table mapping `SC`, `Y1`, `Y1A`, `Y3A`, `Y4`, `Y5`, `Y6B` to species/strain names and any relevant isolation notes.
- Experimental design (current dataset): `S-1` and `S-dup` are duplicate measurements of the same fermentation sample/run (technical duplicates, not independent fermentations). If independent fermentation reruns exist in the future, they must be encoded explicitly (e.g., `run_id`, `batch_id`) so the app can enable replicate-based inference.
- Replicate interpretation (v1): treat `S-1` and `S-dup` as technical duplicates used for QC and for computing a single summary value (mean by default) per analyte/time/strain; do not treat them as independent biological replicates.
- Fermentation conditions: media composition, inoculum, temperature, agitation/aeration, pH control, and any deviations.
- Measurement metadata: method descriptions for each assay (e.g., GC method, DNS method). If detection/quantification limits (LOD/LOQ) are not available, the app must state that explicitly in outputs.
- Objective coverage metadata: a simple "what was measured vs not measured" list for each objective item (e.g., maltose tolerance, invertase activity). The app should surface this as "Not available in current dataset" rather than implying results.

## Cross-Cutting Data Issues to Resolve (before any statistics)
These are mandatory for peer review.

1. Canonical sample IDs: treat `Y6B` consistently across all datasets (normalize any legacy variant to `Y6B`). Keep a lookup table with columns `raw_id`, `canonical_id`, `notes`.
2. Canonical timepoints: parse times from `sample-h` and GC headers robustly. Standardize to integer `time_h` with units documented.
3. Technical duplicates: `S-1` and `S-dup` must be retained in long form. Use them for QC (agreement, CV%) and to compute a single summary value per analyte/time/strain (mean by default). Do not treat `S-1` and `S-dup` as independent biological replicates.
3a. Replicate agreement QC: provide per-measure `S-1` vs `S-dup` comparisons (scatter + correlation, and/or Bland-Altman) and a list of rows failing predefined agreement thresholds.
4. Units and variable naming: normalize `Type-unit` variants (e.g., `Sugar (g/L)` vs `sugar (g/L)`). Maintain a dictionary: `analyte`, `unit`, `method`, `sheet_source`.
5. Volatile compound name standardization: fix systematic inconsistencies (case, spacing, typos, non-ASCII hyphens). Implement a controlled mapping table and apply it before aggregation/testing.
6. Duplicates within a timepoint: some GC timepoints contain duplicate compound names. This must be a user-selectable option (sum, mean, or keep distinct) and must be recorded in outputs for reproducibility.
7. No co-culture samples in current dataset: hide/disable co-culture comparison analyses for now. If co-cultures are added later, encode them in sample IDs (e.g., `Y1+Y3A`) or provide separate metadata.
8. Non-detects / zeros (v1 policy; no LOD/LOQ used):

- For concentration-like endpoints (metabolites, acids/sugars, enzymes, volatiles RA), `0` means "not detected" (non-detect), not a true measured zero. This must not be imputed. Minimum requirements:
  - Treat non-detects as missing (NA) for numeric concentration summaries/models and clearly report `n_detected` and `n_total`.
  - Always include a detection view (percent detected) alongside concentration summaries.
  - For log-scale plots/models, apply log transforms to detected (positive) values only. If non-detects must be shown on the same figure, show them as a separate symbol/category rather than forcing them onto a log scale.
  - State explicitly in outputs that LOD/LOQ values are not available for this dataset and are not used.
Exception (tolerance broth/agar): in `tolerance test-broth`, `0` means "no growth" (real 0) and must be retained as 0 (not treated as NA/non-detect). In `tolerance test-agar`, `agar-score = 0` is also a real "no growth" score.

LOD/LOQ definitions (for reporting context only):

- LOD (limit of detection): smallest amount distinguishable from measurement noise.
- LOQ (limit of quantification): smallest amount that can be quantified with acceptable accuracy/precision.

9. Statistical inference gating (current dataset): independent fermentation replicates (biological repeats) are not encoded; `S-1`/`S-dup` are technical duplicates. Therefore, the app must not present p-values/CIs as evidence of biological between-strain differences. If users request p-values, the app may compute them only as an explicitly labeled, opt-in "technical duplicates / exploratory" view (measurement-level repeatability), with strong wording that biological inference requires independent replicate IDs. If future datasets include `run_id`/`batch_id`/biological replicate IDs, inferential statistics can be enabled conditionally.

Name standardization rules (must be explicit for peer review):

- Unicode normalize to avoid look-alike characters (e.g., replace non-ASCII hyphens with `-` so `1-Heptanol` and `1\u2011Heptanol` are treated as the same compound).
- Trim leading/trailing whitespace and collapse repeated spaces.
- Case normalization for matching (e.g., treat `Ethanol` and `ethanol` as the same), while preserving a canonical display name from the mapping table.
- Apply a curated correction map for known typos/synonyms (examples observed in the GC sheets: `2-Ocatanone` -> `2-Octanone`, `Isoamylalcohol` -> `Isoamyl alcohol`, `Ethyl Acetate` -> `Ethyl acetate`).

User-selectable analysis options (must be recorded in outputs for reproducibility):

- GC duplicate handling within a timepoint: allow the user to select one of:
- `sum`: sum duplicate values for the same canonical compound within a timepoint
- `mean`: mean of duplicate values within a timepoint
- `keep_distinct`: keep duplicates as separate features (requires a deterministic naming rule such as appending `#1`, `#2` based on row order)
- The chosen option must be shown in the UI and written into downloads/reports.

## Analysis Tabs (Shiny app v1: implemented analyses)

### Tab 1: Growth Dynamics (Objective 1)
Scope: mono-culture only (no co-culture samples available).

Primary endpoints:

- Microbial growth: `microbial growth Log (CFU/mL)` (from `R/Metabolites.xlsx / Yeast metabolites`)
- Survival: `microbial-survival- Log (CFU/mL)` (from `R/Yeast tolerance test.xlsx / Autolytic-full`)
- Biomass proxy: `OD` (from `Autolytic-full`)

Implemented analyses:

1. Technical duplicate handling (required): `S-1` and `S-dup` are treated as repeat measurements of the same sample/run (technical duplicates), not independent fermentations.
- User-selectable replicate handling for plots/tables: `Mean(S-1,S-dup)` (default), `S-1 only`, `S-dup only`, `Difference (S-1 - S-dup)` (QC).
- Optional overlay of raw technical-duplicate points (detected only) when plotting the mean.
2. Core visualizations:
- Time-series line plot: endpoint vs time with strain colors.
- Bar chart over time with user-selectable stacking (`stack` vs `dodge`) for strain comparison.
- Snapshot comparison at a selected timepoint (dot plot or bar plot).
- Delta plot between two selected timepoints (`t2 - t1`) per strain.
- Heatmap (strain x time).
3. Summary table (downloadable): wide table including `S-1`, `S-dup`, `mean_value`, `sd_value`, `diff_s1_sdup`, and `cv_pct` by `strain` and `time_h`.
4. Optional p-values (opt-in; technical duplicates only; exploratory, not biological inference):
- Compare strains at each timepoint: global test + all-pairs comparisons, with within-method adjustment and an overall adjustment across the displayed pairwise tests.
- Compare timepoints within each strain: global test + all-pairs comparisons, with the same adjustment policy.
- Methods: `ANOVA + Tukey` or `Kruskal + Wilcoxon` (user-selectable).
5. Downloads: time-series plot (PNG), snapshot plot (PNG), p-values (CSV), summary table (XLSX), plus growth-model fit plot (PNG) and model outputs (XLSX).

6. Growth model fitting (exploratory; no biological replication in current dataset):
- Fit standard growth models per strain (user-selectable): Logistic, Gompertz, Baranyi.
- Report fit metrics and derived curve summaries per strain/model: AIC (where applicable), RMSE, convergence status, `max_pred`, `mu_max`, `lag_tangent`, `t_inflect`, `t_95%`.
- Provide prediction plots (observed vs fitted) and residual/diagnostic plots (residuals vs time and vs fitted) and flag failed fits.
- Peer-review constraint: with the current dataset (technical duplicates only), treat model fits and any between-strain comparisons as descriptive/exploratory only. Confidence intervals / parameter testing are deferred until independent fermentation replicate IDs (e.g., `run_id`) are added.
- Implementation note (v1): growth-model fitting is disabled when replicate handling is `Difference (S-1 - S-dup)`.

### Tab 2: Stress Tolerance (Objective 2)
Data sources:

- `tolerance test-agar` (ordinal score)
- `tolerance test-broth` (quantitative)

Analyses:

1. Agar tolerance (ordinal): heatmap of `agar-score` by strain and concentration for the selected analyte group (descriptive only; no p-values in v1).
2. Broth tolerance (quantitative; 0 = real "no growth"):
- Heatmap by strain and concentration for the selected analyte group.
- Dose-response plot with user controls:
  - Geometry: line chart or bar chart.
  - If bars: stacked vs dodged for strain comparison.
  - If concentrations can be parsed numerically, the x-axis uses numeric concentration; otherwise it is treated as an ordered factor.
3. Technical duplicates (broth only): user-selectable replicate handling: mean / S-1 / S-dup / diff (QC).
4. Optional broth p-values (opt-in; technical duplicates only; exploratory):
- For each concentration (each heatmap column): global test + all-pairs comparisons across strains.
- Methods: `ANOVA + Tukey` or `Kruskal + Wilcoxon` (user-selectable).
- Multiple-testing: within-method adjustment plus an overall adjustment across displayed pairwise tests; optional concentration filter for readability.
5. Downloads: heatmap (PNG), dose-response (PNG; broth), tolerance index plot (PNG; broth), tolerance index outputs (XLSX), p-values (CSV; broth), and table (XLSX).
6. Composite tolerance scoring (transparent "tolerance index"; broth-only; exploratory):
- For broth (quantitative) per selected analyte group, compute strain-level summary metrics:
  - Mean response across all concentrations/levels.
  - Worst-case performance (minimum response across concentrations).
  - AUC across concentrations/levels (only if concentrations can be parsed numerically, with a single unit).
  - IC50 (linear interpolation; only if numeric concentrations are available and response crosses 50% of the strain's max).
- Provide strain rankings with:
  - Explicit weights (preset options + custom sliders).
  - Automatic gating: if AUC/IC50 are not computable for the selected analyte group, those metrics are excluded and weights are renormalized.
  - Sensitivity analysis: rank changes across presets (plot + table).
- Peer-review constraints:
  - Composite scores are descriptive/exploratory only in v1 (technical duplicates; no biological replication).
  - All rules + weights used are visible in the UI and exported in downloads.

Peer-review notes:

- Objective mentions citric acid, maltose, fructose specifically; confirm whether those conditions exist. If not, document as a data gap rather than inferring.

### Tab 3: Enzymatic Activity (Objective 3)
Data sources:

- `Enzyme production` sheet (quantitative; includes activity per mL and per mg protein)
- `tolerance test-agar` contains some enzyme-screen style entries (qualitative/ordinal)

Analyses:

v1 scope decision: only the `Enzyme production rate` category is shown in the app. `Enzyme production` and `Enzyme degradation` are excluded from Objective 3 UI/outputs.

1. Production-rate comparisons across strains (per selected measure):
- Bar chart (ordered) for the selected measure.
- Optional overlay of technical duplicate points (S-1/S-dup) when plotting the mean.
2. Multi-measure overview:
- Strain x measure heatmap for the production-rate category.
- User-selectable scaling: `Raw` or `Z-score per measure` (recommended for different magnitudes/units).
3. Optional exploratory p-values (technical duplicates only; opt-in):
- Global test + all-pairs comparisons across strains for:
  - the selected measure, or
  - all measures (with multiple-testing adjustment applied and recorded).
4. Time comparisons are not available in v1 because no time variable is encoded in the enzyme sheet.
5. Downloads: plot (PNG), heatmap (PNG), p-values (CSV), and summary table (XLSX).

Peer-review notes:

- Objective lists invertase, pectinase, polygalacturonase, cellulase. Map the available assays (`PG`, `PNL`, `PME`, beta-glucosidase, protease) to these categories explicitly. If invertase/cellulase are not directly measured, state this and avoid over-claiming.

### Tab 4: Volatile Compounds (GC) (Objective 4)
Scope: mono-culture only (no co-culture samples available).

Data sources:

- `R/GC_VolatilesYeast.xlsx` (time series per strain; strain sheets)
- `R/exampleMapData.xlsx` / `compound_map` (mapping coverage only; this file is an example and incomplete)

Implementation note (v1): `QC-All` exists as a sheet in `R/GC_VolatilesYeast.xlsx`, but the v1 loader excludes sheets starting with `QC` and no background subtraction is performed.

Analyses:

1. Parsing + standardization (required for peer review):
- Robust parsing of wide GC sheets into a tidy table (`strain`, `time_h`, `compound_std`, `ra_ug_L`, `detected`).
- Compound-name standardization: casing/whitespace/unicode-hyphen normalization and a curated typo/synonym map (e.g., `Ethanol` = `ethanol`; `1\u2011Heptanol` = `1-Heptanol`).
2. Duplicate compound names within a timepoint (user-selectable; recorded in outputs):
- `sum`: sum duplicates
- `mean`: mean duplicates
- `keep_distinct`: keep duplicates as separate features with deterministic suffix `#1/#2`
3. Descriptive profiling (no replicates; no p-values/CIs):
- Heatmap of Top N compounds per timepoint (Top N slider: 1-100), with strain focus:
  - Single strain (cleaner), or
  - All strains (faceted; user-selectable facet column count).
- Compare view:
  - Trajectories over time for selected strains/compounds (facet by compound or by strain).
  - Transform options for plotting: `sqrt` (default), `log1p`, or `none`.
  - Optional "show non-detects as 0" for plotting only (explicitly labeled; not imputation for inference).
  - Snapshot at a selected timepoint with stacked or dodged bars.
- PCA (exploratory only):
  - 0-fill for non-detects for ordination only; drop constant features.
  - Transform options: `log1p`, `sqrt`, `none`; optional feature scaling.
  - Features: all compounds or Top N by mean RA.
  - Outputs: PCA scatter with variance explained + a "Top loadings" table.
- Summary:
  - Total volatile load over time (sum RA; detected values only).
  - Compound richness over time (# detected > 0).
- Mapping coverage:
  - Mapping coverage table using `R/exampleMapData.xlsx` (coverage only; mapping is incomplete by design and must be expanded for full sensory/class reporting).
4. Downloads: heatmap (PNG), compare plot (PNG), PCA (PNG), PCA loadings (CSV), and table (XLSX).

Important constraint for peer review (current dataset):

- GC volatiles are single measurements per strain/timepoint (no replicates). The app must not present compound-level p-values or confidence intervals. Restrict this tab to descriptive and exploratory analyses (profiles, trajectories, PCA/clustering) and clearly label conclusions as exploratory.

### Tab 5: Metabolic Dynamics (Objective 5; mono-culture only)
Scope: the current dataset does not contain co-culture samples, so co-culture synergy/antagonism is deferred. v1 focuses on within-strain, mono-culture metabolic dynamics across the sheets available in the provided workbooks.

Implemented analyses:

1. Data source selection (mono-culture only):
- `R/Metabolites.xlsx / Yeast metabolites` (metabolites time series; technical duplicates).
- `R/Yeast tolerance test.xlsx / Autolytic-full` (physiology panel over time; technical duplicates).
- `R/Metabolites.xlsx / Reducing sugar-24h` (24h reducing sugar summaries; technical duplicates).
2. Selection: user selects an analyte appropriate to the chosen data source.
3. Non-detect handling (v1 policy): for concentration-like values, `0` is treated as non-detect and excluded from means (no imputation). Summary tables report detection counts/rates.
4. Technical duplicates: user-selectable replicate handling for plots/tables: mean / S-1 / S-dup / diff (QC). Optional raw-point overlay (detected only) when plotting the mean.
5. Core visualizations (per selected analyte):
- Dynamics time series (mean detected values by strain over time; points-only when there is only one timepoint).
- Timepoint compare: snapshot across strains at a selected timepoint (dot or bar).
- Delta plot: change between two selected timepoints (`t2 - t1`) by strain.
- Heatmap (strain x time) for the selected analyte.
6. Multivariate (PCA; exploratory):
- Available for `Yeast metabolites` and `Autolytic-full` only; disabled for `Reducing sugar-24h` (single timepoint) and when replicate handling is `Difference`.
- Options: transform (`log1p`, `sqrt`, none), feature selection (all vs Top N by mean), scaling, and optional 0-fill for non-detects (ordination only).
- Outputs: PCA scatter (variance explained) + Top loadings table.
7. Yield/Efficiency summaries (descriptive; `Yeast metabolites` only):
- User selects substrate and product analytes (from `Type-unit | Analyte`) and two timepoints (`t1`, `t2`).
- Outputs:
  - Scatter: product formed (`delta product`) vs substrate consumed (`-delta substrate`).
  - Bar chart: yield ratio = (`delta product`) / (`consumed substrate`) where defined.
- Peer-review constraint: uses detected values only (no imputation). If required analytes/timepoints are missing or non-detect, results show as “Not available in current dataset” (no placeholders).
8. Optional p-values (opt-in; technical duplicates only; exploratory):
- Compare strains at each timepoint, or compare timepoints within each strain.
- Global test + all-pairs comparisons, with within-method adjustment and overall adjustment across displayed tests.
- Methods: `ANOVA + Tukey` or `Kruskal + Wilcoxon` (user-selectable).
9. Downloads: dynamics (PNG), compare (PNG), heatmap (PNG), PCA (PNG), PCA loadings (CSV), yield plot (PNG), yield table (XLSX), p-values (CSV), and summary table (XLSX).

### Tab 6: Best Performing Yeast Decision Framework (Objective 6)
This tab must make the selection criteria explicit and reproducible.

Analyses (expanded; peer-review-safe):

1. Performance rubric (explicit): define the exact metrics computed from available datasets and show which are available vs missing.
   - Primary composite metrics (v1): growth (max CFU), tolerance (broth mean across conditions), enzymes (mean across production-rate measures), volatiles (mean total RA across time; depends on GC duplicate policy).
   - Supporting context metrics (v1, descriptive): growth endpoint, survival endpoint, broth worst-case (min), agar mean/min score, enzyme measure count, volatile richness and endpoint total.
2. Multi-criteria ranking (exploratory):
   - Scale each metric to a 0-1 score across strains.
   - Compute a weighted-sum composite score and rank strains.
   - Provide a full component breakdown: raw metrics, scaled metrics, weights used (after normalization), weighted contributions, and total score (all downloadable).
3. Weight presets + sensitivity analysis (required):
   - Provide preset weight sets (Balanced, Equal, Growth-focused, Tolerance-focused, Enzymes-focused, Volatiles-focused) plus custom sliders.
   - Sensitivity outputs: rank changes across presets (plot + table), recorded in downloads.
4. Trade-off visualization:
   - Scatter plot of two user-selected metrics (raw or scaled 0-1) to visualize competing objectives without collapsing everything into a single score.
5. Peer-review constraints (v1):
   - Current datasets contain technical duplicates only (no independent fermentation replicates encoded). Ranking is exploratory and must not be presented with statistical confidence claims.
   - No co-culture comparisons or synergy/antagonism scoring in v1.

### QC & Downloads Tab (technical duplicates)
The app includes a dedicated QC tab for assessing agreement between technical duplicates (`S-1` vs `S-dup`) across datasets.

Implemented analyses:

1. Dataset selector: Metabolites, Autolytic, Tolerance broth, or Enzymes.
2. Technical-duplicate agreement plot:
- Scatter plot of `S-1` vs `S-dup` for the selected endpoint/measure.
- Annotation on the plot: `n`, Pearson correlation (`r`), and linear model equation `y = a + b x` (intercept `a`, slope `b`).
- Reference lines: fitted regression line and `y = x` agreement line.
3. Download: QC plot (PNG).

## Peer-Review Readiness Checklist (must satisfy before releasing results)
1. Methods reproducibility: every figure/table generated from scripted analysis with versioned dependencies; exportable session info and analysis parameters.
2. Transparent data handling: document all cleaning rules (ID/time/name mapping, duplicate handling, non-detect policy). If QC subtraction is used in the future, document it explicitly (v1 does not subtract `QC-All` for GC).
3. Appropriate statistical design: account for repeated measures over time and replicate structure; check assumptions; use robust/non-parametric alternatives where required.
4. Multiplicity control (when inferential testing is enabled): for compound-level testing, apply FDR correction and report adjusted p-values.
5. Effect sizes + uncertainty (when inferential testing is enabled): report effect sizes with confidence intervals, not p-values alone.
6. Diagnostics (when inferential models are used): provide model diagnostics (residuals, influence/outlier checks) for inferential models.
7. No over-claiming: avoid causal language unless the experimental design supports it; treat supervised classification carefully (no claims of \"separation\" without validation).

## Deferred / Future Enhancements (not implemented in v1)
These are intentionally deferred either because the current dataset does not support them (e.g., no co-cultures; no biological replicate IDs; no GC replicates), or because they require additional design decisions.

1. Co-culture analyses (Objectives 1/4/5): any mono vs co-culture comparisons, synergy/antagonism claims, and interaction scoring are deferred until co-culture samples/IDs are added.
2. Biological replicate-aware inference: enable mixed-effects / repeated-measures models only if independent fermentation replicate IDs (e.g., `run_id`/`batch_id`) are added. v1 p-values are technical-duplicate/exploratory only.
3. Objective 1 growth model inference upgrades: model fitting is implemented in v1 as descriptive/exploratory. Parameter uncertainty (CIs), between-strain parameter comparisons, and replicate-aware models are deferred until independent fermentation replicate IDs are added.
4. Objective 2 tolerance index upgrades: the composite tolerance index is implemented in v1 as descriptive/exploratory. Model-based dose-response fitting (nonlinear models, IC50 CIs) and replicate-aware ranking inference are deferred until independent fermentation replicate IDs are added.
5. Objective 4 GC QC/normalization upgrades: `QC-All` background subtraction, internal-standard normalization, and mapping-driven sensory/class aggregation require explicit rules plus a complete mapping table (the current `R/exampleMapData.xlsx` is an example only).
6. Objective 5 cross-omics and advanced multivariate: correlations between metabolites and volatiles, joint trajectory/clustering across multiple panels, and replicate-aware time-series models are deferred until (a) objective-specific metadata is encoded and (b) independent fermentation replicate IDs exist.
7. Report bundling: a one-click report (Quarto/Rmd) that packages all figures/tables plus parameter logs (and vector formats like PDF/SVG) is deferred (v1 supports per-tab PNG/CSV/XLSX downloads and exports session info + analysis parameters).
