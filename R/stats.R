## Statistical helpers.
##
## IMPORTANT (v1):
## - The current dataset contains only technical duplicates (S-1/S-dup) for most endpoints.
## - Any p-values computed from these are "technical duplicate" statistics and must not be
##   presented as biological inference unless independent fermentation replicate IDs exist.

apply_techrep_handling <- function(df, group_cols, handling = c("mean", "s1", "sdup", "diff"), value_col = "value") {
  handling <- match.arg(handling)
  stopifnot(is.data.frame(df), all(group_cols %in% names(df)), "tech_rep" %in% names(df), value_col %in% names(df))

  if (handling == "s1") {
    # If a dataset accidentally contains duplicates within the same group+tech_rep,
    # collapse them deterministically (mean) so downstream pivots/plots remain stable.
    return(df |>
      dplyr::filter(.data$tech_rep == "S-1") |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), tech_rep = .data$tech_rep) |>
      dplyr::summarise(
        value_use = if (all(is.na(.data[[value_col]]))) NA_real_ else mean(.data[[value_col]], na.rm = TRUE),
        .groups = "drop"
      ))
  }
  if (handling == "sdup") {
    return(df |>
      dplyr::filter(.data$tech_rep == "S-dup") |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), tech_rep = .data$tech_rep) |>
      dplyr::summarise(
        value_use = if (all(is.na(.data[[value_col]]))) NA_real_ else mean(.data[[value_col]], na.rm = TRUE),
        .groups = "drop"
      ))
  }
  if (handling == "mean") {
    return(df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarise(
        value_use = if (all(is.na(.data[[value_col]]))) NA_real_ else mean(.data[[value_col]], na.rm = TRUE),
        sd_use = if (sum(!is.na(.data[[value_col]])) < 2) NA_real_ else stats::sd(.data[[value_col]], na.rm = TRUE),
        n_total = dplyr::n(),
        n_nonmissing = sum(!is.na(.data[[value_col]])),
        .groups = "drop"
      ))
  }
  # diff: S-1 - S-dup (QC/diagnostic)
  wide <- df |>
    dplyr::select(dplyr::all_of(group_cols), tech_rep = .data$tech_rep, v = .data[[value_col]]) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), tech_rep = .data$tech_rep) |>
    dplyr::summarise(v = if (all(is.na(.data$v))) NA_real_ else mean(.data$v, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = tech_rep, values_from = v)
  if (!all(c("S-1", "S-dup") %in% names(wide))) {
    wide$`S-1` <- NA_real_
    wide$`S-dup` <- NA_real_
  }
  wide |>
    dplyr::mutate(value_use = .data$`S-1` - .data$`S-dup`)
}

as_factor_safe <- function(x) {
  if (is.factor(x)) x else factor(x, levels = unique(x))
}

safe_aov_tukey <- function(df, response_col, group_col) {
  stopifnot(is.data.frame(df), response_col %in% names(df), group_col %in% names(df))
  dat <- df |>
    dplyr::transmute(
      y = suppressWarnings(as.numeric(.data[[response_col]])),
      g = as_factor_safe(.data[[group_col]])
    ) |>
    dplyr::filter(!is.na(.data$y) & !is.na(.data$g))

  if (nrow(dat) < 3 || length(unique(dat$g)) < 2) {
    return(list(global_p = NA_real_, pairwise = tibble::tibble()))
  }

  fit <- stats::aov(y ~ g, data = dat)
  an <- summary(fit)[[1]]
  global_p <- as.numeric(an["g", "Pr(>F)"] %||% NA_real_)

  tk <- tryCatch(stats::TukeyHSD(fit, which = "g"), error = function(e) NULL)
  if (is.null(tk) || is.null(tk$g)) {
    return(list(global_p = global_p, pairwise = tibble::tibble()))
  }

  pw <- as.data.frame(tk$g, stringsAsFactors = FALSE)
  pw$comparison <- rownames(pw)
  rownames(pw) <- NULL

  pw <- tibble::as_tibble(pw) |>
    dplyr::rename(
      diff = .data$diff,
      lwr = .data$lwr,
      upr = .data$upr,
      p_adj = .data$`p adj`
    ) |>
    dplyr::mutate(
      group1 = stringr::str_split_fixed(.data$comparison, "-", 2)[, 1],
      group2 = stringr::str_split_fixed(.data$comparison, "-", 2)[, 2]
    )

  list(global_p = global_p, pairwise = pw)
}

safe_kw_wilcox <- function(df, response_col, group_col, p_adjust_method = "BH") {
  stopifnot(is.data.frame(df), response_col %in% names(df), group_col %in% names(df))
  dat <- df |>
    dplyr::transmute(
      y = suppressWarnings(as.numeric(.data[[response_col]])),
      g = as_factor_safe(.data[[group_col]])
    ) |>
    dplyr::filter(!is.na(.data$y) & !is.na(.data$g))

  if (nrow(dat) < 3 || length(unique(dat$g)) < 2) {
    return(list(global_p = NA_real_, pairwise = tibble::tibble()))
  }

  global_p <- tryCatch(stats::kruskal.test(y ~ g, data = dat)$p.value, error = function(e) NA_real_)

  pw <- tryCatch(
    stats::pairwise.wilcox.test(dat$y, dat$g, p.adjust.method = p_adjust_method, exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(pw) || is.null(pw$p.value)) {
    return(list(global_p = global_p, pairwise = tibble::tibble()))
  }

  m <- pw$p.value
  rn <- rownames(m)
  cn <- colnames(m)
  out <- list()
  k <- 0L
  for (i in seq_along(rn)) {
    for (j in seq_along(cn)) {
      pv <- m[i, j]
      if (!is.na(pv)) {
        k <- k + 1L
        out[[k]] <- tibble::tibble(
          group1 = rn[[i]],
          group2 = cn[[j]],
          p_adj = as.numeric(pv)
        )
      }
    }
  }

  list(global_p = global_p, pairwise = dplyr::bind_rows(out))
}

parse_condition_numeric <- function(x) {
  # Parse numeric concentration from strings like "Glucose-15%", "Temperature45oC", "Acidity-pH2", "Acidity-0.6%".
  x <- as.character(x)
  x <- normalize_hyphens(x)
  x_trim <- stringr::str_trim(x)

  # pH
  m_ph <- stringr::str_match(x_trim, "(?i)ph\\s*([0-9]+(?:\\.[0-9]+)?)")
  ph_val <- suppressWarnings(as.numeric(m_ph[, 2]))

  # percent
  m_pct <- stringr::str_match(x_trim, "([0-9]+(?:\\.[0-9]+)?)\\s*%")
  pct_val <- suppressWarnings(as.numeric(m_pct[, 2]))

  # temperature
  m_t <- stringr::str_match(x_trim, "(?i)temp(?:erature)?\\s*([0-9]+(?:\\.[0-9]+)?)")
  t_val <- suppressWarnings(as.numeric(m_t[, 2]))

  # fallback: first numeric token
  m_any <- stringr::str_match(x_trim, "([0-9]+(?:\\.[0-9]+)?)")
  any_val <- suppressWarnings(as.numeric(m_any[, 2]))

  unit <- dplyr::case_when(
    !is.na(ph_val) ~ "pH",
    !is.na(t_val) ~ "C",
    !is.na(pct_val) ~ "%",
    !is.na(any_val) ~ "num",
    TRUE ~ NA_character_
  )

  val <- dplyr::case_when(
    !is.na(ph_val) ~ ph_val,
    !is.na(t_val) ~ t_val,
    !is.na(pct_val) ~ pct_val,
    !is.na(any_val) ~ any_val,
    TRUE ~ NA_real_
  )

  tibble::tibble(condition_raw = x_trim, condition_value = val, condition_unit = unit)
}

trapz_auc <- function(x, y) {
  # Trapezoidal AUC for numeric x/y, ignoring NA pairs.
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  o <- order(x)
  x <- x[o]
  y <- y[o]
  sum(diff(x) * (y[-length(y)] + y[-1]) / 2)
}

ic50_linear_interp <- function(x, y, target = 0.5) {
  # Compute a simple IC50/EC50-like threshold by linear interpolation on normalized response.
  # Returns NA if x is non-numeric, if there is no crossing, or if data are insufficient.
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) return(NA_real_)

  o <- order(x)
  x <- x[o]
  y <- y[o]

  y_max <- max(y, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) return(NA_real_)
  y_norm <- y / y_max

  if (any(abs(y_norm - target) < 1e-12)) {
    return(x[which.min(abs(y_norm - target))])
  }

  i <- which(y_norm <= target)[1]
  if (is.na(i) || i <= 1) return(NA_real_)
  j <- i - 1

  x1 <- x[j]
  x2 <- x[i]
  y1 <- y_norm[j]
  y2 <- y_norm[i]
  if (!is.finite(x1) || !is.finite(x2) || !is.finite(y1) || !is.finite(y2) || y2 == y1) return(NA_real_)

  x1 + (target - y1) * (x2 - x1) / (y2 - y1)
}

scale01_vec <- function(x) {
  # Map numeric vector to [0,1] range with stable handling of constants/NA.
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[[1]]) || !is.finite(rng[[2]])) return(rep(NA_real_, length(x)))
  if (rng[[2]] - rng[[1]] < .Machine$double.eps) return(rep(0.5, length(x)))
  (x - rng[[1]]) / (rng[[2]] - rng[[1]])
}

normalize_weights <- function(w, metric_avail) {
  # Normalize named numeric weights to sum to 1, zeroing unavailable metrics.
  w <- w[intersect(names(w), names(metric_avail))]
  w[is.na(w)] <- 0
  w[!metric_avail[names(w)]] <- 0
  if (sum(w) == 0) {
    # Deterministic fallback if everything is unavailable or zeroed.
    w <- stats::setNames(rep(0, length(w)), names(w))
    if (length(w) > 0) w[[1]] <- 1
  }
  w / sum(w)
}

derive_curve_metrics <- function(t_grid, y_grid, y_baseline = NA_real_, frac = 0.95) {
  t_grid <- suppressWarnings(as.numeric(t_grid))
  y_grid <- suppressWarnings(as.numeric(y_grid))
  ok <- is.finite(t_grid) & is.finite(y_grid)
  t_grid <- t_grid[ok]
  y_grid <- y_grid[ok]
  if (length(t_grid) < 3) {
    return(tibble::tibble(
      max_pred = NA_real_,
      t_inflect = NA_real_,
      mu_max = NA_real_,
      lag_tangent = NA_real_,
      t_frac = NA_real_
    ))
  }

  o <- order(t_grid)
  t_grid <- t_grid[o]
  y_grid <- y_grid[o]

  max_pred <- max(y_grid, na.rm = TRUE)
  dy <- diff(y_grid) / diff(t_grid)
  mu_max <- if (length(dy) == 0) NA_real_ else max(dy, na.rm = TRUE)
  i_inf <- if (length(dy) == 0) NA_integer_ else which.max(dy) + 1L
  t_inflect <- if (is.na(i_inf) || i_inf < 1 || i_inf > length(t_grid)) NA_real_ else t_grid[[i_inf]]
  y_inflect <- if (is.na(i_inf) || i_inf < 1 || i_inf > length(y_grid)) NA_real_ else y_grid[[i_inf]]

  if (!is.finite(y_baseline)) y_baseline <- min(y_grid, na.rm = TRUE)
  lag_tangent <- if (is.finite(mu_max) && mu_max > 0 && is.finite(t_inflect) && is.finite(y_inflect)) {
    t_inflect - (y_inflect - y_baseline) / mu_max
  } else {
    NA_real_
  }

  t_frac <- NA_real_
  if (is.finite(max_pred) && max_pred > 0 && is.finite(frac) && frac > 0 && frac < 1) {
    thr <- frac * max_pred
    idx <- which(y_grid >= thr)[1]
    if (!is.na(idx)) t_frac <- t_grid[[idx]]
  }

  tibble::tibble(
    max_pred = max_pred,
    t_inflect = t_inflect,
    mu_max = mu_max,
    lag_tangent = lag_tangent,
    t_frac = t_frac
  )
}

baranyi_fun <- function(t, y0, ymax, mu, lag) {
  # Baranyi & Roberts growth model (log-scale form).
  # This implementation is for descriptive fits only in v1 (technical duplicates only).
  t <- suppressWarnings(as.numeric(t))
  mu <- pmax(suppressWarnings(as.numeric(mu)), 1e-8)
  lag <- pmax(suppressWarnings(as.numeric(lag)), 0)

  A <- t + (1 / mu) * log(exp(-mu * t) + exp(-mu * lag) - exp(-mu * (t + lag)))
  y0 + mu * A - log(1 + (exp(mu * A) - 1) / exp(ymax - y0))
}

fit_growth_models <- function(df, models = c("logistic", "gompertz", "baranyi"), grid_n = 200L) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c("strain", "time_h", "y") %in% names(df)))
  models <- intersect(models, c("logistic", "gompertz", "baranyi"))
  if (length(models) == 0) {
    return(list(fits = tibble::tibble(), pred = tibble::tibble(), resid = tibble::tibble()))
  }

  fits_out <- list()
  pred_out <- list()
  resid_out <- list()
  k_fit <- 0L
  k_pred <- 0L
  k_resid <- 0L

  for (s in sort(unique(df$strain))) {
    dat <- df |>
      dplyr::filter(.data$strain == s) |>
      dplyr::transmute(t = suppressWarnings(as.numeric(.data$time_h)), y = suppressWarnings(as.numeric(.data$y))) |>
      dplyr::filter(is.finite(.data$t) & is.finite(.data$y)) |>
      dplyr::arrange(.data$t)
    if (nrow(dat) < 4) next

    y_base <- min(dat$y, na.rm = TRUE)
    t_min <- min(dat$t, na.rm = TRUE)
    t_max <- max(dat$t, na.rm = TRUE)
    t_grid <- seq(t_min, t_max, length.out = as.integer(grid_n))

    for (m in models) {
      fit <- NULL
      params <- list(
        Asym = NA_real_,
        xmid = NA_real_,
        scal = NA_real_,
        b2 = NA_real_,
        b3 = NA_real_,
        y0 = NA_real_,
        ymax = NA_real_,
        mu = NA_real_,
        lag = NA_real_
      )

      if (identical(m, "logistic")) {
        fit <- tryCatch(
          stats::nls(y ~ stats::SSlogis(t, Asym, xmid, scal), data = dat, control = stats::nls.control(warnOnly = TRUE, maxiter = 200)),
          error = function(e) NULL
        )
      } else if (identical(m, "gompertz")) {
        fit <- tryCatch(
          stats::nls(y ~ stats::SSgompertz(t, Asym, b2, b3), data = dat, control = stats::nls.control(warnOnly = TRUE, maxiter = 200)),
          error = function(e) NULL
        )
      } else if (identical(m, "baranyi")) {
        y0_start <- min(dat$y, na.rm = TRUE)
        ymax_start <- max(dat$y, na.rm = TRUE)
        mu_start <- max((ymax_start - y0_start) / (t_max - t_min + 1e-6), 1e-3)
        lag_start <- 0.1 * (t_max - t_min)
        fit <- tryCatch(
          stats::nls(
            y ~ baranyi_fun(t, y0, ymax, mu, lag),
            data = dat,
            start = list(y0 = y0_start, ymax = ymax_start, mu = mu_start, lag = lag_start),
            algorithm = "port",
            lower = c(y0 = y0_start - 2, ymax = ymax_start, mu = 1e-6, lag = 0),
            upper = c(y0 = ymax_start, ymax = ymax_start + 5, mu = 10, lag = t_max),
            control = stats::nls.control(warnOnly = TRUE, maxiter = 300)
          ),
          error = function(e) NULL
        )
      }

      if (is.null(fit)) next

      co <- tryCatch(stats::coef(fit), error = function(e) NULL)
      if (!is.null(co) && length(co) > 0) {
        for (nm in names(co)) params[[nm]] <- as.numeric(co[[nm]])
      }

      y_hat <- tryCatch(stats::predict(fit, newdata = data.frame(t = dat$t)), error = function(e) rep(NA_real_, nrow(dat)))
      rmse <- sqrt(mean((dat$y - y_hat)^2, na.rm = TRUE))
      aic <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)

      y_grid <- tryCatch(stats::predict(fit, newdata = data.frame(t = t_grid)), error = function(e) rep(NA_real_, length(t_grid)))
      met <- derive_curve_metrics(t_grid, y_grid, y_baseline = y_base, frac = 0.95)

      k_fit <- k_fit + 1L
      fits_out[[k_fit]] <- tibble::tibble(
        strain = s,
        model = m,
        n_points = nrow(dat),
        aic = as.numeric(aic),
        rmse = as.numeric(rmse),
        Asym = params$Asym,
        xmid = params$xmid,
        scal = params$scal,
        b2 = params$b2,
        b3 = params$b3,
        y0 = params$y0,
        ymax = params$ymax,
        mu = params$mu,
        lag = params$lag
      ) |>
        dplyr::bind_cols(met)

      k_pred <- k_pred + 1L
      pred_out[[k_pred]] <- tibble::tibble(
        strain = s,
        model = m,
        time_h = t_grid,
        y_pred = y_grid
      )

      k_resid <- k_resid + 1L
      resid_out[[k_resid]] <- tibble::tibble(
        strain = s,
        model = m,
        time_h = dat$t,
        y_obs = dat$y,
        y_fit = y_hat,
        resid = dat$y - y_hat
      )
    }
  }

  list(
    fits = dplyr::bind_rows(fits_out),
    pred = dplyr::bind_rows(pred_out),
    resid = dplyr::bind_rows(resid_out)
  )
}
