# Refit the mortality model and reproduce the three figures.
# Run: Rscript code/Figure_5.R
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript code/Figure_5.R")
project <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/")))
required <- c("data.table", "fixest", "ggplot2", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
if (!capabilities("cairo")) stop("This workflow requires Cairo graphics support.")
suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE, scipen = 999)
run_figure <- function() {
  work <- tempfile("figure5_")
  dir.create(work)
  old_wd <- getwd()
  on.exit({setwd(old_wd); unlink(work, recursive = TRUE)}, add = TRUE)
  if (!file.copy(file.path(project, "input"), work, recursive = TRUE)) stop("Could not copy inputs.")
  setwd(work)
  for (dir in c("analysis", "figures", "source_data")) dir.create(file.path("output", dir), recursive = TRUE)
  local({
    dat <- data.table::fread("input/historical_panel.csv", data.table = FALSE)
    dat <- dat[order(dat$iso3, dat$t), , drop = FALSE]
    stopifnot(nrow(dat) == 450L, length(unique(dat$iso3)) == 30L,
              !anyDuplicated(dat[c("iso3", "t")]))
    dat$iso3 <- factor(dat$iso3)
    dat$region <- factor(dat$region)
    dat$log_death <- ifelse(dat$death > 0, log(dat$death), NA_real_)
    dat$t_center <- dat$t - min(dat$t, na.rm = TRUE)

    primary_max_lag <- 6L
    forecast_cutoff <- 2018L
    test_years <- 2019:2021
    controls <- c("ln_gdp", "ln_pop", "ln_health_spend", "temp", "pm25")

    for (k in 0:primary_max_lag) {
      nm <- paste0("hamr_l", k)
      dat[[nm]] <- unlist(lapply(split(dat$hamr, dat$iso3), function(x) {
        if (k == 0L) return(x)
        c(rep(NA_real_, k), x[seq_len(max(0L, length(x) - k))])
      }), use.names = FALSE)
    }
    lag_cols <- function(k) paste0("hamr_l", 0:k)
    model_data <- function(data, k) {
      req <- c("death", "log_death", lag_cols(k), controls, "region")
      data[complete.cases(data[, req, drop = FALSE]), , drop = FALSE]
    }
    formula_string <- function(outcome, k, include_covariates = TRUE) {
      rhs <- lag_cols(k)
      if (include_covariates) rhs <- c(rhs, controls)
      paste0(outcome, " ~ ", paste(rhs, collapse = " + "),
             " + fixest::i(region, t_center) | iso3")
    }
    forecast_fit <- function(k, include_covariates = TRUE) {
      req <- c("death", "log_death", lag_cols(k), controls, "region")
      train <- dat[dat$t <= forecast_cutoff, , drop = FALSE]
      test <- dat[dat$t %in% test_years, , drop = FALSE]
      train <- train[complete.cases(train[, req, drop = FALSE]), , drop = FALSE]
      test <- test[complete.cases(test[, req, drop = FALSE]), , drop = FALSE]
      m <- fixest::feols(as.formula(formula_string("log_death", k, include_covariates)),
                         data = train, cluster = ~iso3, notes = FALSE)
      test$pred_log <- as.numeric(stats::predict(m, newdata = test))
      smear <- mean(exp(stats::resid(m)), na.rm = TRUE)
      test$pred_death <- exp(test$pred_log) * smear
      metrics <- data.frame(model = if (include_covariates) "adjusted" else "common-sample unadjusted",
                            lag_order = k, train_n = nrow(train), test_n = nrow(test),
                            rmse_death = sqrt(mean((test$death - test$pred_death)^2)),
                            mae_death = mean(abs(test$death - test$pred_death)),
                            rmse_log = sqrt(mean((test$log_death - test$pred_log)^2)),
                            mae_log = mean(abs(test$log_death - test$pred_log)),
                            smear_factor = smear, stringsAsFactors = FALSE)
      list(model = m, test = test, metrics = metrics)
    }

    forecast_adjusted <- lapply(0:primary_max_lag, forecast_fit, include_covariates = TRUE)
    forecast_unadjusted <- lapply(0:primary_max_lag, forecast_fit, include_covariates = FALSE)
    forecast_metrics <- rbind(do.call(rbind, lapply(forecast_adjusted, function(x) x$metrics)),
                              do.call(rbind, lapply(forecast_unadjusted, function(x) x$metrics)))
    write.csv(forecast_metrics, "output/analysis/covariate_forecast_backtest_metrics.csv", row.names = FALSE, fileEncoding = "UTF-8")
    forecast_predictions <- do.call(rbind, c(
      lapply(forecast_adjusted, function(x) transform(x$test, model = "adjusted", lag_order = x$metrics$lag_order)),
      lapply(forecast_unadjusted, function(x) transform(x$test, model = "common-sample unadjusted", lag_order = x$metrics$lag_order))
    ))
    write.csv(forecast_predictions, "output/analysis/covariate_forecast_backtest_predictions.csv", row.names = FALSE, fileEncoding = "UTF-8")
    best_row <- forecast_metrics[forecast_metrics$model == "adjusted", ][which.min(forecast_metrics[forecast_metrics$model == "adjusted", "rmse_log"]), ]
    best_k <- as.integer(best_row$lag_order)
    write.csv(best_row, "output/analysis/selected_adjusted_forecast_model.csv", row.names = FALSE, fileEncoding = "UTF-8")

    final_data <- model_data(dat, best_k)
    final_model <- fixest::feols(as.formula(formula_string("log_death", best_k)),
                                 data = final_data, cluster = ~iso3, notes = FALSE)
    final_smear <- mean(exp(stats::resid(final_model)), na.rm = TRUE)
    stopifnot(stats::nobs(final_model) == 330L)
    saveRDS(list(model = final_model, smear_factor = final_smear, max_lag = best_k,
                 controls = controls, base_year = min(dat$t, na.rm = TRUE)),
            "output/analysis/final_adjusted_forecast_bundle.rds")

    plot_theme <- ggplot2::theme_classic(base_size = 10, base_family = "Arial") +
      ggplot2::theme(legend.position = "bottom", axis.title = ggplot2::element_text(size = 10),
                     plot.title = ggplot2::element_text(face = "bold", size = 11))
    p_forecast <- ggplot2::ggplot(forecast_metrics, ggplot2::aes(x = lag_order, y = rmse_log, color = model)) +
      ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = c("adjusted" = "#B2182B", "common-sample unadjusted" = "#2166AC")) +
      ggplot2::labs(x = "Maximum HAMR lag (years)", y = "Test RMSE of log(death)",
                    color = NULL) + plot_theme

    grDevices::cairo_pdf("output/figures/covariate_forecast_backtest.pdf", width = 7.2, height = 4.8)
    print(p_forecast)
    invisible(grDevices::dev.off())
  })

  local({
    bundle_file <- "output/analysis/final_adjusted_forecast_bundle.rds"
    history_file <- "input/historical_panel.csv"
    hamr_file <- "input/hamr_scenario_trajectories.csv"
    scenarios <- c("SSP1", "SSP2", "SSP3", "SSP5")
    hamr_map <- c(
      "HAMR-anti" = "HAMR_anti",
      "HAMR-current" = "HAMR_current",
      "HAMR-strength-5y" = "HAMR_strength_5y",
      "HAMR-strength-10y" = "HAMR_strength_10y"
    )
    years_all <- 2007:2050

    bundle <- readRDS(bundle_file)
    model <- bundle$model
    if (as.integer(bundle$max_lag) != 4L) {
      stop("The supplied final model does not use HAMR lag 0-4.")
    }
    country_levels <- names(fixest::fixef(model)$iso3)

    history <- data.table::fread(history_file, data.table = FALSE)
    history$iso3 <- toupper(trimws(as.character(history$iso3)))
    history$year <- as.integer(history$t)
    history <- history[history$iso3 %in% country_levels & history$year %in% 2007:2021, , drop = FALSE]
    history <- history[, c("iso3", "year", "death", "hamr", "ln_gdp", "ln_pop",
                           "ln_health_spend", "temp", "pm25", "region"), drop = FALSE]
    history$region <- as.character(history$region)
    if (anyNA(history$region)) stop("Historical covariate data contain missing regions.")
    region_levels <- sort(unique(history$region))

    # Equivalent dummy-variable version of the final model, used only to obtain
    # confidence intervals. The point predictions still come from the saved final
    # fixed-effects model above.
    ci_history <- history
    ci_history$iso3 <- factor(ci_history$iso3, levels = country_levels)
    ci_history$region <- factor(ci_history$region, levels = region_levels)
    ci_history$log_death <- ifelse(ci_history$death > 0, log(ci_history$death), NA_real_)
    ci_history$t_center <- ci_history$year - as.numeric(bundle$base_year)
    for (k in 0:as.integer(bundle$max_lag)) {
      nm <- paste0("hamr_l", k)
      ci_history[[nm]] <- ave(ci_history$hamr, ci_history$iso3, FUN = function(x) {
        if (k == 0L) x else c(rep(NA_real_, k), x[seq_len(max(0L, length(x) - k))])
      })
    }
    ci_formula <- stats::as.formula(
      "log_death ~ hamr_l0 + hamr_l1 + hamr_l2 + hamr_l3 + hamr_l4 + ln_gdp + ln_pop + ln_health_spend + temp + pm25 + fixest::i(region, t_center) + fixest::i(iso3)"
    )
    ci_required <- c("log_death", bundle$controls, "region", "iso3", paste0("hamr_l", 0:as.integer(bundle$max_lag)))
    ci_history <- ci_history[complete.cases(ci_history[, ci_required, drop = FALSE]), , drop = FALSE]
    ci_model <- fixest::feols(ci_formula, data = ci_history, cluster = ~iso3, notes = FALSE)

    # HAMR trajectories: use the conditional country-specific trajectory supplied
    # by the source file, then retain the countries identifiable by the final model.
    hamr_raw <- data.table::fread(hamr_file, data.table = FALSE)
    hamr_raw$ISO3 <- toupper(trimws(as.character(hamr_raw$ISO3)))
    hamr_raw$year <- as.integer(hamr_raw$year)
    hamr_raw$scenario <- toupper(trimws(as.character(hamr_raw$scenario)))
    hamr_raw <- hamr_raw[
      hamr_raw$scenario %in% scenarios & hamr_raw$year %in% years_all &
        hamr_raw$ISO3 %in% country_levels, , drop = FALSE
    ]
    if (!nrow(hamr_raw)) stop("No usable HAMR rows remain after filtering.")

    hamr_long <- do.call(rbind, lapply(names(hamr_map), function(label) {
      z <- hamr_raw[, c("ISO3", "year", "scenario", hamr_map[[label]])]
      names(z)[ncol(z)] <- "hamr"
      z$hamr_type <- label
      z
    }))
    hamr_long$iso3 <- toupper(trimws(as.character(hamr_long$ISO3)))
    hamr_long$year <- as.integer(hamr_long$year)
    hamr_long$hamr <- as.numeric(hamr_long$hamr)
    hamr_long <- hamr_long[, c("iso3", "year", "scenario", "hamr_type", "hamr"), drop = FALSE]
    hamr_long <- hamr_long[!is.na(hamr_long$hamr), , drop = FALSE]

    expected_hamr <- length(country_levels) * length(years_all) * length(scenarios) * length(hamr_map)
    if (nrow(hamr_long) != expected_hamr) {
      stop("Unexpected HAMR rows: found ", nrow(hamr_long), "; expected ", expected_hamr,
           ". Check scenario, prediction_type, or country coverage.")
    }


    future_cov <- data.table::fread("input/future_covariates.csv", data.table = FALSE)
    stopifnot(nrow(future_cov) == 30L * 29L * 4L,
              !anyDuplicated(future_cov[c("iso3", "year", "scenario")]),
              all(complete.cases(future_cov[c(bundle$controls, "region")])))
    historical_cov <- merge(
      history[, c("iso3", "year", "death", "ln_gdp", "ln_pop", "ln_health_spend", "temp", "pm25", "region")],
      data.frame(scenario = scenarios, stringsAsFactors = FALSE),
      by = NULL, sort = FALSE
    )
    historical_cov <- historical_cov[, c("iso3", "year", "scenario", "death", "ln_gdp", "ln_pop",
                                         "ln_health_spend", "temp", "pm25", "region"), drop = FALSE]
    cov_panel <- rbind(historical_cov, future_cov[, names(historical_cov), drop = FALSE])
    cov_panel <- cov_panel[order(cov_panel$scenario, cov_panel$iso3, cov_panel$year), , drop = FALSE]

    panel <- merge(hamr_long, cov_panel, by = c("iso3", "year", "scenario"), all = FALSE, sort = FALSE)
    panel$period <- ifelse(panel$year <= 2021, "historical_fit", "future_projection")
    panel <- panel[order(panel$scenario, panel$hamr_type, panel$iso3, panel$year), , drop = FALSE]
    group_id <- interaction(panel$scenario, panel$hamr_type, panel$iso3, drop = TRUE)
    for (k in 0:as.integer(bundle$max_lag)) {
      nm <- paste0("hamr_l", k)
      panel[[nm]] <- ave(panel$hamr, group_id, FUN = function(x) {
        if (k == 0L) x else c(rep(NA_real_, k), x[seq_len(max(0L, length(x) - k))])
      })
    }

    panel$iso3 <- factor(panel$iso3, levels = country_levels)
    panel$region <- factor(panel$region, levels = region_levels)
    panel$t_center <- panel$year - as.numeric(bundle$base_year)
    panel$log_death <- ifelse(!is.na(panel$death) & panel$death > 0,
                              log(panel$death), NA_real_)
    model_vars <- c(bundle$controls, "region", paste0("hamr_l", 0:as.integer(bundle$max_lag)))
    available <- complete.cases(panel[, model_vars, drop = FALSE])
    panel$pred_log_death <- NA_real_
    panel$pred_death <- NA_real_
    panel$pred_death_low95 <- NA_real_
    panel$pred_death_high95 <- NA_real_
    if (any(available)) {
      newdata <- panel[available, , drop = FALSE]
      panel$pred_log_death[available] <- as.numeric(stats::predict(model, newdata = newdata))
      ci_pred <- stats::predict(ci_model, newdata = newdata, se.fit = TRUE,
                                interval = "confidence", level = 0.95)
      smear <- as.numeric(bundle$smear_factor)
      panel$pred_death[available] <- exp(panel$pred_log_death[available]) * smear
      panel$pred_death_low95[available] <- exp(ci_pred$ci_low) * smear
      panel$pred_death_high95[available] <- exp(ci_pred$ci_high) * smear
    }
    message("Prediction unavailable for the first ", bundle$max_lag,
            " years because pre-2007 HAMR lags are not supplied.")
    panel$prediction_available <- available
    panel$iso3 <- as.character(panel$iso3)
    panel$region <- as.character(panel$region)

    country_out <- panel[, c("iso3", "year", "period", "scenario", "hamr_type", "hamr",
                             "hamr_l0", "hamr_l1", "hamr_l2", "hamr_l3", "hamr_l4",
                             "ln_gdp", "ln_pop", "ln_health_spend", "temp", "pm25", "region",
                             "death", "prediction_available", "pred_log_death", "pred_death",
                             "pred_death_low95", "pred_death_high95"), drop = FALSE]
    country_out <- country_out[order(country_out$scenario, country_out$hamr_type,
                                     country_out$iso3, country_out$year), , drop = FALSE]
    write.csv(country_out, "output/analysis/death_predictions_2007_2050_by_country.csv", row.names = FALSE,
              fileEncoding = "UTF-8")

    annual <- data.table::as.data.table(country_out)[,
      .(pred_death_total = if (all(is.na(pred_death))) NA_real_ else sum(pred_death, na.rm = TRUE),
        pred_death_low95 = if (all(is.na(pred_death_low95))) NA_real_ else sum(pred_death_low95, na.rm = TRUE),
        pred_death_high95 = if (all(is.na(pred_death_high95))) NA_real_ else sum(pred_death_high95, na.rm = TRUE),
        pred_death_mean = if (all(is.na(pred_death))) NA_real_ else mean(pred_death, na.rm = TRUE),
        observed_death_total = if (all(is.na(death))) NA_real_ else sum(death, na.rm = TRUE),
        n_country_predicted = sum(!is.na(pred_death)),
        n_country_total = data.table::uniqueN(iso3)),
      by = .(scenario, hamr_type, year, period)]
    annual <- as.data.frame(annual[order(annual$scenario, annual$hamr_type, annual$year), ])
    write.csv(annual, "output/analysis/death_predictions_2007_2050_annual_summary.csv", row.names = FALSE,
              fileEncoding = "UTF-8")
  })

  local({
    annual_file <- "output/analysis/death_predictions_2007_2050_annual_summary.csv"
    annual <- data.table::fread(annual_file, encoding = "UTF-8")
    annual[, year := as.integer(year)]
    annual[, pred_death_total := as.numeric(pred_death_total)]

    scenarios <- c("SSP1", "SSP2", "SSP3", "SSP5")
    hamr_levels <- c("HAMR-anti", "HAMR-current", "HAMR-strength-5y", "HAMR-strength-10y")

    baseline <- annual[scenario %in% scenarios & hamr_type == "HAMR-current" & year == 2019,
                       .(baseline_death_2019 = pred_death_total), by = scenario]
    target <- annual[scenario %in% scenarios & hamr_type %in% hamr_levels & year == 2050,
                     .(scenario, hamr_type, death_2050 = pred_death_total)]
    delta <- merge(target, baseline, by = "scenario", all.x = TRUE)
    delta[, `:=`(
      delta_death = death_2050 - baseline_death_2019,
      delta_percent = 100 * (death_2050 - baseline_death_2019) / baseline_death_2019
    )]
    delta[, hamr_type := factor(hamr_type, levels = hamr_levels)]
    data.table::setorder(delta, scenario, hamr_type)
    data.table::fwrite(delta, "output/source_data/death_2050_change_vs_2019_summary.csv",
                       bom = TRUE)

    method_labels <- c("HAMR-anti" = "anti", "HAMR-current" = "current",
                       "HAMR-strength-5y" = "strength 5y",
                       "HAMR-strength-10y" = "strength 10y")

    theme_nature <- ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
      ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
        axis.text = ggplot2::element_text(size = 6.5, colour = "black"),
        axis.title = ggplot2::element_text(size = 7),
        axis.text.x = ggplot2::element_text(size = 6.2, colour = "black"),
        panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, colour = "grey85",
                                                   linetype = "dotted"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 6.5),
        legend.text = ggplot2::element_text(size = 6.2),
        plot.title = ggplot2::element_text(size = 8, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 6.5),
        plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "pt")
      )

    min_delta <- min(delta$delta_percent, na.rm = TRUE)
    max_delta <- max(delta$delta_percent, na.rm = TRUE)
    y_min <- floor(min(min_delta, 0) / 5) * 5 - 5
    y_max <- ceiling(max(max_delta, 0) / 5) * 5 + 5

    # Single-panel version: four SSP groups, each with the four HAMR trajectories.
    # This follows the grouped-bar logic of delta_SSP_onepanel.pdf.
    ssp_labels <- c("SSP1" = "SSP1", "SSP2" = "SSP2",
                    "SSP3" = "SSP3", "SSP5" = "SSP5")
    one_panel_palette <- c(
      "HAMR-anti" = "#5B6472",
      "HAMR-current" = "#277DA1",
      "HAMR-strength-10y" = "#D9852B",
      "HAMR-strength-5y" = "#149C7E"
    )
    one_panel_levels <- c(
      "HAMR-current", "HAMR-anti", "HAMR-strength-10y", "HAMR-strength-5y"
    )
    delta_one_panel <- data.table::copy(delta)
    delta_one_panel[, scenario := factor(scenario, levels = scenarios,
                                        labels = unname(ssp_labels[scenarios]))]
    delta_one_panel[, hamr_type := factor(as.character(hamr_type),
                                          levels = one_panel_levels)]
    bar_position <- ggplot2::position_dodge(width = 0.76)

    fig_one_panel <- ggplot2::ggplot(delta_one_panel,
                                     ggplot2::aes(x = scenario, y = delta_percent,
                                                  fill = hamr_type)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.45, colour = "grey25") +
      ggplot2::geom_col(position = bar_position, width = 0.76,
                        colour = "black", linewidth = 0.35) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%+.1f%%", delta_percent),
                     vjust = ifelse(delta_percent >= 0, -0.35, 1.25)),
        position = bar_position, size = 2.7, colour = "grey25"
      ) +
      ggplot2::scale_fill_manual(values = one_panel_palette, breaks = one_panel_levels,
                                 labels = unname(method_labels[one_panel_levels]), drop = FALSE) +
      ggplot2::scale_y_continuous(limits = c(y_min, y_max),
                                  breaks = seq(y_min, y_max, by = 10),
                                  labels = function(x) format(x, trim = TRUE, scientific = FALSE),
                                  expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
      ggplot2::labs(
        x = "SSP pathway", y = "Change in predicted deaths (%)",
        fill = "Policy strength"
      ) +
      theme_nature +
      ggplot2::theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        plot.title = ggplot2::element_text(size = 10, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 7)
      )

    one_panel_stem <- "output/figures/death_change_2050_vs_2019"
    grDevices::cairo_pdf(paste0(one_panel_stem, ".pdf"), width = 183 / 25.4,
                         height = 120 / 25.4, family = "Arial")
    print(fig_one_panel)
    grDevices::dev.off()
  })

  local({
    annual_file <- "output/analysis/death_predictions_2007_2050_annual_summary.csv"
    d <- data.table::fread(annual_file, encoding = "UTF-8")
    d[, year := as.integer(year)]
    d[, c("pred_death_total", "pred_death_low95", "pred_death_high95") :=
        lapply(.SD, as.numeric), .SDcols = c("pred_death_total", "pred_death_low95", "pred_death_high95")]
    d <- d[year >= 2011 & is.finite(pred_death_total)]

    scenarios <- c("SSP1", "SSP2", "SSP3", "SSP5")
    policies <- c("HAMR-anti", "HAMR-current", "HAMR-strength-5y", "HAMR-strength-10y")
    if (data.table::uniqueN(d$scenario) != length(scenarios) ||
        data.table::uniqueN(d$hamr_type) != length(policies)) {
      stop("The annual prediction file must contain all four SSP pathways and HAMR trajectories.")
    }

    policy_colours <- c("HAMR-anti" = "#5B6472", "HAMR-current" = "#277DA1",
                        "HAMR-strength-10y" = "#D9852B", "HAMR-strength-5y" = "#149C7E")
    ssp_shapes <- c("SSP1" = 16, "SSP2" = 15, "SSP3" = 17, "SSP5" = 18)
    # Show markers only in the projection period, at biennial years plus 2050.
    markers <- d[year >= 2022 & (year %% 2 == 0 | year == 2050)]
    d_bottom <- d[hamr_type %in% c("HAMR-anti", "HAMR-strength-5y")]
    d_strength10 <- d[hamr_type == "HAMR-strength-10y"]
    d_current <- d[hamr_type == "HAMR-current"]
    markers_bottom <- markers[hamr_type %in% c("HAMR-anti", "HAMR-strength-5y")]
    markers_strength10 <- markers[hamr_type == "HAMR-strength-10y"]
    markers_current <- markers[hamr_type == "HAMR-current"]

    # The horizontal target is the same reference used in the original four-panel
    # no-CI figure: 90% of each SSP's HAMR-current predicted total in 2019.
    target_lines <- d[year == 2019 & hamr_type == "HAMR-current",
                      .(target_death = 0.9 * pred_death_total), by = scenario]

    y_min <- floor(min(d$pred_death_total, na.rm = TRUE) / 50000) * 50000
    y_max <- ceiling(max(d$pred_death_total, na.rm = TRUE) / 50000) * 50000
    ci_y_min <- floor(min(c(d$pred_death_total, d$pred_death_low95), na.rm = TRUE) / 50000) * 50000
    ci_y_max <- ceiling(max(c(d$pred_death_total, d$pred_death_high95), na.rm = TRUE) / 50000) * 50000

    theme_figure <- ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
      ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
        axis.title = ggplot2::element_text(size = 8),
        axis.text = ggplot2::element_text(size = 7, colour = "black"),
        legend.title = ggplot2::element_text(size = 7.2, face = "bold"),
        legend.text = ggplot2::element_text(size = 6.8),
        legend.key.width = grid::unit(15, "pt"),
        legend.key.height = grid::unit(10, "pt"),
        plot.title = ggplot2::element_text(size = 9.5, face = "bold", hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 7.2, colour = "#4D4D4D", hjust = 0),
        plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 6.5, "pt")
      )

    base_plot <- ggplot2::ggplot(d, ggplot2::aes(
      x = year, y = pred_death_total, group = interaction(scenario, hamr_type)
    )) +
      ggplot2::geom_hline(data = target_lines,
                          ggplot2::aes(yintercept = target_death),
                          inherit.aes = FALSE, linetype = "dashed", linewidth = 0.4,
                          colour = "grey35") +
      ggplot2::geom_vline(xintercept = 2022, linewidth = 0.4, linetype = "dashed", colour = "#667085") +
      ggplot2::annotate("text", x = 2022.3,
                        y = y_max - 0.03 * (y_max - y_min), label = "Projection",
                        hjust = 0, vjust = 1, size = 2.5, colour = "#4D596A") +
      ggplot2::scale_x_continuous(limits = c(2011, 2050),
                                  breaks = c(2011, 2020, 2022, 2030, 2040, 2050),
                                  expand = ggplot2::expansion(mult = c(0.01, 0.015))) +
      ggplot2::labs(x = "Year", y = "Predicted deaths (30-country total)") +
      theme_figure

    # Policy-colour / SSP-shape version; the same explicit layering is used.
    fig_policy_colour_ssp_shape <- base_plot +
      ggplot2::geom_line(data = d_bottom, ggplot2::aes(colour = hamr_type), linewidth = 0.62, alpha = 0.82) +
      ggplot2::geom_point(data = markers_bottom, ggplot2::aes(colour = hamr_type, shape = scenario), size = 1.4, stroke = 0.25, alpha = 0.82) +
      ggplot2::geom_line(data = d_strength10, ggplot2::aes(colour = hamr_type), linewidth = 0.62, alpha = 0.82) +
      ggplot2::geom_point(data = markers_strength10, ggplot2::aes(colour = hamr_type, shape = scenario), size = 1.4, stroke = 0.25, alpha = 0.82) +
      ggplot2::geom_line(data = d_current, ggplot2::aes(colour = hamr_type), linewidth = 0.62, alpha = 0.82) +
      ggplot2::geom_point(data = markers_current, ggplot2::aes(colour = hamr_type, shape = scenario), size = 1.4, stroke = 0.25, alpha = 0.82) +
      ggplot2::scale_colour_manual(values = policy_colours, name = "Policy strength") +
      ggplot2::scale_shape_manual(values = ssp_shapes, name = "SSP pathway") +
      ggplot2::guides(
        colour = ggplot2::guide_legend(order = 1, override.aes = list(shape = NA, linewidth = 1.0)),
        shape = ggplot2::guide_legend(order = 2, override.aes = list(colour = "#343434", size = 2.3))
      )

    # Matched 95% CI versions. Ribbons are added beneath every line, marker and
    # reference line; the explicit line layering above is therefore retained.
    ci_data <- d[is.finite(pred_death_low95) & is.finite(pred_death_high95)]
    ci_y_scale <- ggplot2::scale_y_continuous(
      limits = c(ci_y_min, ci_y_max), breaks = seq(ci_y_min, ci_y_max, by = 50000),
      labels = scales::label_number(big.mark = ",", accuracy = 1),
      expand = ggplot2::expansion(mult = c(0.03, 0.075))
    )
    policy_ci_ribbon <- ggplot2::geom_ribbon(
      data = ci_data,
      ggplot2::aes(x = year, ymin = pred_death_low95, ymax = pred_death_high95,
                   fill = hamr_type, group = interaction(scenario, hamr_type)),
      inherit.aes = FALSE, alpha = 0.075, colour = NA
    )

    fig_policy_colour_ssp_shape_ci <- fig_policy_colour_ssp_shape
    fig_policy_colour_ssp_shape_ci$layers <- c(list(policy_ci_ribbon), fig_policy_colour_ssp_shape_ci$layers)
    fig_policy_colour_ssp_shape_ci <- fig_policy_colour_ssp_shape_ci +
      ci_y_scale +
      ggplot2::scale_fill_manual(values = policy_colours, guide = "none")


    data.table::fwrite(d, "output/source_data/death_prediction_trajectories.csv")
    ggplot2::ggsave("output/figures/death_predictions_2007_2050i.pdf",
                    fig_policy_colour_ssp_shape_ci, width = 183, height = 120,
                units = "mm", bg = "white", device = grDevices::cairo_pdf, family = "Arial")
  })

  output <- file.path(project, "output")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  figures <- c("covariate_forecast_backtest.pdf", "death_change_2050_vs_2019.pdf", "death_predictions_2007_2050i.pdf")
  if (!all(file.copy(file.path("output", "figures", figures), output, overwrite = TRUE))) stop("Could not copy all figures.")
  cat("Created three figures in output/.\n")
}
run_figure()
