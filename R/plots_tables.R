## Shared plot/table helpers for publication-ready outputs.

khadijat_theme <- function(base_size = 12, base_family = "Segoe UI") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(face = "bold", color = "#1d2733"),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = grid::unit(3, "pt"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold")
    )
}

strain_palette <- function(strains) {
  # Stable, colorblind-friendly-ish palette.
  pal <- c(
    SC = "#1b9e77",
    Y1 = "#d95f02",
    Y1A = "#7570b3",
    Y3A = "#e7298a",
    Y4 = "#66a61e",
    Y5 = "#e6ab02",
    Y6B = "#a6761d"
  )
  s <- unique(canonical_strain_id(strains))
  out <- pal[s]
  out[is.na(out)] <- "#666666"
  stats::setNames(out, s)
}

summarise_tech_dups <- function(df, value_col = "value", group_cols) {
  # Summarise technical duplicates (S-1/S-dup) as mean of detected values only.
  # 0 has already been converted to NA in loaders.
  stopifnot(value_col %in% names(df))
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_detected = sum(!is.na(.data[[value_col]])),
      mean_value = if (all(is.na(.data[[value_col]]))) NA_real_ else mean(.data[[value_col]], na.rm = TRUE),
      sd_value = if (sum(!is.na(.data[[value_col]])) < 2) NA_real_ else stats::sd(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(p_detected = dplyr::if_else(.data$n_total > 0, .data$n_detected / .data$n_total, NA_real_))
}

download_table_csv <- function(df, filename_stub) {
  shiny::downloadHandler(
    filename = function() paste0(filename_stub, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      readr::write_csv(df, file, na = "")
    }
  )
}

download_table_xlsx <- function(df, filename_stub, sheet = "results") {
  shiny::downloadHandler(
    filename = function() paste0(filename_stub, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content = function(file) {
      writexl::write_xlsx(setNames(list(df), sheet), path = file)
    }
  )
}

save_plot <- function(p, file, width = 9, height = 5, dpi = 300) {
  ext <- tolower(tools::file_ext(file))
  if (ext == "png") {
    ggplot2::ggsave(file, plot = p, width = width, height = height, dpi = dpi, device = ragg::agg_png, bg = "white")
  } else if (ext == "svg") {
    ggplot2::ggsave(file, plot = p, width = width, height = height, device = svglite::svglite, bg = "white")
  } else if (ext == "pdf") {
    ggplot2::ggsave(file, plot = p, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  } else {
    ggplot2::ggsave(file, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  }
}

download_plot_handler <- function(plot_reactive, filename_stub, width = 9, height = 5) {
  shiny::downloadHandler(
    filename = function() paste0(filename_stub, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      p <- plot_reactive()
      save_plot(p, file = file, width = width, height = height, dpi = 320)
    }
  )
}
