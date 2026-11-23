# Refit the model and create Figure 4a and 4b from the supplied inputs.
# Run: Rscript code/Figure_4.R
args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(args) != 1L) stop("Run with Rscript code/Figure_4.R")
project <- dirname(dirname(normalizePath(sub("^--file=", "", args), winslash = "/")))
run_figure <- function() {
  root <- tempfile("figure4_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  if (!file.copy(file.path(project, "input"), root, recursive = TRUE)) stop("Could not copy input files.")
  options(stringsAsFactors = FALSE, digits = 12)

  # Refit the full linear-year within-between random-intercept model and create
  # the scenario-year estimates used by both article figures.
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("Package 'nlme' is required. Install it before running this script.")
  }

  input_dir <- file.path(root, "input")
  output_dir <- file.path(root, "output")
  analysis_dir <- file.path(output_dir, "analysis")
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

  historical <- read.csv(file.path(input_dir, "historical_analysis_data_2007_2021.csv"), check.names = FALSE)
  future_raw <- read.csv(file.path(input_dir, "future_scenario_covariates_2022_2050.csv"), check.names = FALSE)

  raw_vars <- c("lnGDP", "pop", "lnhealth", "PM25", "temp")
  required_historical <- c("country", "ISO3", "region", "year", "HAMR_anti", "HAMR_current",
                           "policy_year", "policy_time", "ATT", raw_vars)
  required_future <- c("country", "ISO3", "region", "year", "scenario", raw_vars)
  if (!all(required_historical %in% names(historical))) stop("Historical input is missing required columns.")
  if (!all(required_future %in% names(future_raw))) stop("Future input is missing required columns.")
  if (anyDuplicated(historical[c("country", "year")])) stop("Historical input has duplicate country-year records.")
  if (anyDuplicated(future_raw[c("country", "year", "scenario")])) stop("Future input has duplicate country-year-scenario records.")

  country_levels <- sort(unique(historical$country))
  region_levels <- c("Western Europe", "Central-Eastern Europe", "Northern Europe", "Southern Europe")
  model_formula <- HAMR_anti ~ lnGDP_between + lnGDP_within + pop_between + pop_within +
    lnhealth_between + lnhealth_within + PM25_between + PM25_within +
    temp_between + temp_within + region + year_centered

  make_design <- function(data, mean_years) {
    means <- aggregate(data[data$year %in% mean_years, raw_vars],
                       list(country = data$country[data$year %in% mean_years]), mean)
    names(means)[-1] <- paste0(raw_vars, "_between")
    data <- merge(data, means, by = "country", all.x = TRUE, sort = FALSE)
    for (variable in raw_vars) {
      data[[paste0(variable, "_within")]] <- data[[variable]] - data[[paste0(variable, "_between")]]
    }
    data$country <- factor(data$country, levels = country_levels)
    data$region <- factor(data$region, levels = region_levels)
    data$year_centered <- as.integer(data$year) - 2007L
    list(data = data[order(data$country, data$year), ], means = means)
  }

  fit_model <- function(data) {
    nlme::lme(fixed = model_formula, random = ~ 1 | country, data = data,
              method = "REML", na.action = na.fail,
              control = nlme::lmeControl(opt = "optim", returnObject = TRUE, msMaxIter = 300))
  }

  predict_parts <- function(model, newdata) {
    X <- model.matrix(delete.response(terms(model)), newdata)
    coefficients <- nlme::fixed.effects(model)
    X <- X[, names(coefficients), drop = FALSE]
    fixed <- as.numeric(X %*% coefficients)
    random_effects <- nlme::ranef(model)
    random_lookup <- setNames(as.numeric(random_effects[, 1]), row.names(random_effects))
    list(X = X, fixed = fixed,
         conditional = fixed + unname(random_lookup[as.character(newdata$country)]))
  }

  metric <- function(observed, predicted, type) {
    data.frame(prediction_type = type, N = length(observed),
               RMSE = sqrt(mean((observed - predicted)^2)),
               MAE = mean(abs(observed - predicted)), R2 = cor(observed, predicted)^2,
               bias = mean(predicted - observed), observed_SD = sd(observed))
  }

  # Temporal validation is retained because it is the gate used before the final fit.
  validation_design <- make_design(historical, 2007:2018)
  training <- validation_design$data[validation_design$data$year <= 2018, ]
  validation <- validation_design$data[validation_design$data$year >= 2019, ]
  validation_model <- fit_model(training)
  validation_prediction <- predict_parts(validation_model, validation)
  validation_metrics <- rbind(
    metric(validation$HAMR_anti, validation_prediction$conditional, "conditional_country_RE"),
    metric(validation$HAMR_anti, validation_prediction$fixed, "fixed_only")
  )
  validation_metrics$passes_rule <- validation_metrics$RMSE < validation_metrics$observed_SD & validation_metrics$R2 >= 0.50
  write.csv(validation_metrics, file.path(analysis_dir, "temporal_validation_metrics_2019_2021.csv"), row.names = FALSE)

  # Final model: all 2007-2021 observations and full-period country means.
  final_design <- make_design(historical, 2007:2021)
  model <- fit_model(final_design$data)
  saveRDS(model, file.path(analysis_dir, "HAMR_anti_full_linear_year_REML_2007_2021.rds"))
  writeLines(capture.output(summary(model)), file.path(analysis_dir, "model_summary.txt"))

  # The future file holds analysis-ready scenario covariates. Recalculate the
  # within-between decomposition from the final historical country means.
  future <- future_raw[, c("country", "ISO3", "region", "year", "scenario", raw_vars)]
  future <- merge(future, final_design$means, by = "country", all.x = TRUE, sort = FALSE)
  for (variable in raw_vars) {
    future[[paste0(variable, "_within")]] <- future[[variable]] - future[[paste0(variable, "_between")]]
  }
  future$country <- factor(future$country, levels = country_levels)
  future$region <- factor(future$region, levels = region_levels)
  future$year_centered <- as.integer(future$year) - 2007L
  future <- future[order(future$country, match(future$scenario, c("SSP1", "SSP2", "SSP3", "SSP5")), future$year), ]

  prediction <- predict_parts(model, future)
  X <- prediction$X
  fixed_variance <- rowSums((X %*% vcov(model)) * X)
  variance_components <- nlme::VarCorr(model)
  random_intercept_variance <- as.numeric(variance_components[1, "Variance"])
  residual_variance <- as.numeric(variance_components[nrow(variance_components), "Variance"])
  observations_per_country <- table(final_design$data$country)
  posterior_random_variance <- 1 / (1 / random_intercept_variance +
    as.numeric(observations_per_country[as.character(future$country)]) / residual_variance)
  future$pred_fixed_only <- prediction$fixed
  future$pred_conditional_country_RE <- prediction$conditional

  # Existing policy effect and strengthening schedules. ATT is an externally
  # specified input, rather than a parameter estimated in this script.
  att_values <- c(`1` = 1.516142, `2` = 2.427547, `3` = 2.306026, `4` = 3.115974,
                  `5` = 4.078662, `6` = 5.692036, `7` = 7.141248, `8` = 8.077958,
                  `9` = 5.846664, `10` = 6.988011, `11` = 8.27435, `12` = 10.37837,
                  `13` = 6.766425)
  att_function <- function(policy_time) {
    result <- rep(0, length(policy_time))
    active <- !is.na(policy_time) & policy_time >= 1
    exact <- active & policy_time <= 13
    result[exact] <- unname(att_values[as.character(policy_time[exact])])
    result[active & policy_time >= 14] <- 8.101789
    result
  }
  reinforcement_total <- function(year, starts) {
    Reduce(`+`, lapply(starts, function(start) att_function(ifelse(year < start, 0, year - start + 1))),
           init = rep(0, length(year)))
  }
  policy_lookup <- unique(historical[c("country", "policy_year")])
  future <- merge(future, policy_lookup, by = "country", all.x = TRUE, sort = FALSE)
  future$policy_time <- ifelse(is.na(future$policy_year) | future$year < future$policy_year,
                               0L, future$year - future$policy_year + 1L)
  future$ATT_base <- att_function(future$policy_time)
  future$ATT_strength_5y <- future$ATT_base + reinforcement_total(future$year, c(2026, 2031, 2036, 2041, 2046))
  future$ATT_strength_10y <- future$ATT_base + reinforcement_total(future$year, c(2026, 2036, 2046))

  make_future_output <- function(prediction_type) {
    predicted <- if (prediction_type == "conditional_country_RE") future$pred_conditional_country_RE else future$pred_fixed_only
    variance <- fixed_variance + if (prediction_type == "conditional_country_RE") posterior_random_variance else 0
    out <- future[c("country", "ISO3", "region", "year", "scenario", "policy_year", "policy_time",
                    "ATT_base", "ATT_strength_5y", "ATT_strength_10y")]
    out$prediction_type <- prediction_type
    out$HAMR_anti <- predicted
    out$HAMR_current <- predicted + out$ATT_base
    out$HAMR_strength_5y <- predicted + out$ATT_strength_5y
    out$HAMR_strength_10y <- predicted + out$ATT_strength_10y
    standard_error <- sqrt(variance)
    for (name in c("HAMR_anti", "HAMR_current", "HAMR_strength_5y", "HAMR_strength_10y")) {
      out[[paste0(name, "_CI_low")]] <- out[[name]] - 1.96 * standard_error
      out[[paste0(name, "_CI_high")]] <- out[[name]] + 1.96 * standard_error
    }
    out
  }
  future_output <- rbind(make_future_output("conditional_country_RE"), make_future_output("fixed_only"))

  # Historical values are observed/reconstructed values; confidence intervals are intentionally absent.
  historical_base <- historical[c("country", "ISO3", "region", "year", "policy_year", "policy_time", "ATT",
                                  "HAMR_anti", "HAMR_current")]
  names(historical_base)[names(historical_base) == "ATT"] <- "ATT_base"
  historical_base$ATT_strength_5y <- historical_base$ATT_base
  historical_base$ATT_strength_10y <- historical_base$ATT_base
  historical_base$HAMR_strength_5y <- historical_base$HAMR_current
  historical_base$HAMR_strength_10y <- historical_base$HAMR_current
  historical_output <- do.call(rbind, lapply(c("conditional_country_RE", "fixed_only"), function(type) {
    do.call(rbind, lapply(c("SSP1", "SSP2", "SSP3", "SSP5"), function(scenario) {
      out <- historical_base
      out$scenario <- scenario
      out$prediction_type <- type
      out
    }))
  }))
  for (name in c("HAMR_anti", "HAMR_current", "HAMR_strength_5y", "HAMR_strength_10y")) {
    historical_output[[paste0(name, "_CI_low")]] <- NA_real_
    historical_output[[paste0(name, "_CI_high")]] <- NA_real_
  }

  summarise_future <- function(data) {
    groups <- split(seq_len(nrow(data)), interaction(data$prediction_type, data$scenario, data$year, drop = TRUE))
    do.call(rbind, lapply(groups, function(index) {
      type <- data$prediction_type[index[1]]
      scenario <- data$scenario[index[1]]
      year <- data$year[index[1]]
      rows <- which(future$scenario == scenario & future$year == year)
      average_X <- colMeans(X[rows, , drop = FALSE])
      mean_variance <- as.numeric(t(average_X) %*% vcov(model) %*% average_X) +
        if (type == "conditional_country_RE") sum(posterior_random_variance[rows]) / length(rows)^2 else 0
      standard_error <- sqrt(mean_variance)
      out <- data.frame(prediction_type = type, scenario = scenario, year = year)
      for (name in c("HAMR_anti", "HAMR_current", "HAMR_strength_5y", "HAMR_strength_10y")) {
        mean_value <- mean(data[[name]][index])
        out[[name]] <- mean_value
        out[[paste0(name, "_CI_low")]] <- mean_value - 1.96 * standard_error
        out[[paste0(name, "_CI_high")]] <- mean_value + 1.96 * standard_error
      }
      out
    }))
  }
  future_summary <- summarise_future(future_output)
  historical_summary <- aggregate(cbind(HAMR_anti, HAMR_current, HAMR_strength_5y, HAMR_strength_10y) ~
    prediction_type + scenario + year, historical_output, mean)
  for (name in c("HAMR_anti", "HAMR_current", "HAMR_strength_5y", "HAMR_strength_10y")) {
    historical_summary[[paste0(name, "_CI_low")]] <- NA_real_
    historical_summary[[paste0(name, "_CI_high")]] <- NA_real_
  }
  scenario_year_summary <- rbind(historical_summary[, names(future_summary)], future_summary)
  scenario_year_summary <- scenario_year_summary[order(scenario_year_summary$prediction_type,
    match(scenario_year_summary$scenario, c("SSP1", "SSP2", "SSP3", "SSP5")), scenario_year_summary$year), ]
  write.csv(scenario_year_summary, file.path(analysis_dir, "scenario_year_mean_HAMR_2007_2050.csv"), row.names = FALSE)

  fixed_effects <- summary(model)$tTable
  fixed_effects <- data.frame(term = row.names(fixed_effects), estimate = fixed_effects[, "Value"],
    std_error = fixed_effects[, "Std.Error"], df = fixed_effects[, "DF"],
    CI_low = fixed_effects[, "Value"] - qt(0.975, fixed_effects[, "DF"]) * fixed_effects[, "Std.Error"],
    CI_high = fixed_effects[, "Value"] + qt(0.975, fixed_effects[, "DF"]) * fixed_effects[, "Std.Error"],
    p_value = fixed_effects[, "p-value"], row.names = NULL)
  write.csv(fixed_effects, file.path(analysis_dir, "model_fixed_effect_coefficients_95CI.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
  cat("Completed model reproduction. Conditional temporal validation passed: ",
      validation_metrics$passes_rule[validation_metrics$prediction_type == "conditional_country_RE"], "\n", sep = "")

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  summary_file <- file.path(root, "output", "analysis", "scenario_year_mean_HAMR_2007_2050.csv")
  figure_dir <- file.path(root, "output", "figures")
  source_dir <- file.path(root, "output", "source_data")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

  wide <- read.csv(summary_file, check.names = FALSE)
  wide <- wide[wide$prediction_type == "conditional_country_RE", ]
  policy_columns <- c("HAMR_anti", "HAMR_current", "HAMR_strength_5y", "HAMR_strength_10y")
  policy_labels <- c("HAMR-anti", "HAMR-current", "HAMR-strength 5y", "HAMR-strength 10y")
  ssp_labels <- c(SSP1 = "SSP1-2.6", SSP2 = "SSP2-4.5", SSP3 = "SSP3-7.0", SSP5 = "SSP5-8.5")

  long <- do.call(rbind, lapply(seq_along(policy_columns), function(i) {
    name <- policy_columns[i]
    data.frame(scenario = wide$scenario, year = as.integer(wide$year), policy = policy_labels[i],
               HAMR = wide[[name]], CI_low = wide[[paste0(name, "_CI_low")]],
               CI_high = wide[[paste0(name, "_CI_high")]])
  }))
  long$scenario <- factor(long$scenario, levels = names(ssp_labels), labels = unname(ssp_labels))
  long$policy <- factor(long$policy, levels = policy_labels)
  long <- long[order(long$policy, long$scenario, long$year), ]
  write.csv(long, file.path(source_dir, "policyColour_SSPshape_CI_source_data.csv"), row.names = FALSE)
  markers <- long[long$year >= 2022 & (long$year %% 2 == 0 | long$year == 2050), ]

  policy_colours <- c("HAMR-anti" = "#5B6472", "HAMR-current" = "#277DA1",
                      "HAMR-strength 5y" = "#149C7E", "HAMR-strength 10y" = "#D9852B")
  ssp_shapes <- c("SSP1-2.6" = 16, "SSP2-4.5" = 15, "SSP3-7.0" = 17, "SSP5-8.5" = 18)
  theme_figure <- ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.title = ggplot2::element_text(size = 8),
      axis.text = ggplot2::element_text(size = 7, colour = "black"),
      legend.title = ggplot2::element_text(size = 7.2, face = "bold"),
      legend.text = ggplot2::element_text(size = 6.8),
      legend.key.width = grid::unit(15, "pt"), legend.key.height = grid::unit(10, "pt"),
      plot.title = ggplot2::element_text(size = 9.5, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 7.2, colour = "#4D4D4D", hjust = 0),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 6.5, "pt")
    )

  figure <- ggplot2::ggplot(long, ggplot2::aes(x = year, y = HAMR, colour = policy,
    group = interaction(scenario, policy))) +
    ggplot2::geom_vline(xintercept = 2022, linewidth = 0.4, linetype = "dashed", colour = "#667085") +
    ggplot2::geom_line(data = long[long$policy != "HAMR-current", ], linewidth = 0.62, alpha = 0.74) +
    ggplot2::geom_point(data = markers[markers$policy != "HAMR-current", ],
      ggplot2::aes(shape = scenario), size = 1.55, stroke = 0.25, alpha = 0.74) +
    ggplot2::geom_line(data = long[long$policy == "HAMR-current", ], linewidth = 0.62, alpha = 0.74) +
    ggplot2::geom_point(data = markers[markers$policy == "HAMR-current", ],
      ggplot2::aes(shape = scenario), size = 1.55, stroke = 0.25, alpha = 0.74) +
    ggplot2::annotate("text", x = 2022.25, y = max(long$HAMR) + 1.5,
      label = "Projection", hjust = 0, vjust = 0, size = 2.5, colour = "#4D596A") +
    ggplot2::scale_colour_manual(values = policy_colours, name = "Policy strength") +
    ggplot2::scale_shape_manual(values = ssp_shapes, name = "SSP pathway") +
    ggplot2::scale_x_continuous(breaks = c(2007, 2015, 2022, 2030, 2040, 2050), expand = c(0.01, 0.015)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.075))) +
    ggplot2::labs(x = "Year", y = "European mean HAMR") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, override.aes = list(shape = NA, linewidth = 1.0)),
      shape = ggplot2::guide_legend(order = 2, override.aes = list(colour = "#343434", size = 2.3))
    ) + theme_figure

  ribbon_other <- ggplot2::geom_ribbon(
    data = long[long$policy != "HAMR-current" & is.finite(long$CI_low) & is.finite(long$CI_high), ],
    ggplot2::aes(x = year, ymin = CI_low, ymax = CI_high, fill = policy,
      group = interaction(scenario, policy)), inherit.aes = FALSE, alpha = 0.075, colour = NA)
  ribbon_current <- ggplot2::geom_ribbon(
    data = long[long$policy == "HAMR-current" & is.finite(long$CI_low) & is.finite(long$CI_high), ],
    ggplot2::aes(x = year, ymin = CI_low, ymax = CI_high, fill = policy,
      group = interaction(scenario, policy)), inherit.aes = FALSE, alpha = 0.11, colour = NA)
  figure <- figure
  figure$layers <- c(list(ribbon_other, ribbon_current), figure$layers)
  figure <- figure + ggplot2::scale_fill_manual(values = policy_colours, guide = "none")

  grDevices::cairo_pdf(file.path(figure_dir, "policyColour_SSPshape_CI.pdf"), width = 183 / 25.4,
                        height = 120 / 25.4, family = "Arial", onefile = TRUE)
  print(figure)
  grDevices::dev.off()
  cat("Created output/figures/policyColour_SSPshape_CI.pdf\n")

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  summary_file <- file.path(root, "output", "analysis", "scenario_year_mean_HAMR_2007_2050.csv")
  figure_dir <- file.path(root, "output", "figures")
  source_dir <- file.path(root, "output", "source_data")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

  wide <- read.csv(summary_file, check.names = FALSE)
  wide <- wide[wide$prediction_type == "conditional_country_RE" & wide$year %in% c(2021, 2050), ]
  outcomes <- c("HAMR_anti", "HAMR_current", "HAMR_strength_10y", "HAMR_strength_5y")
  labels <- c("HAMR-anti", "HAMR-current", "HAMR-strength 10y", "HAMR-strength 5y")
  source_outcomes <- c("HAMR_anti", "HAMR_current", "HAMR_strength_5y", "HAMR_strength_10y")
  ssps <- c("SSP1", "SSP2", "SSP3", "SSP5")
  delta_wide <- do.call(rbind, lapply(ssps, function(scenario) {
    data <- wide[wide$scenario == scenario, ]
    data <- data[match(c(2021, 2050), data$year), ]
    out <- data.frame(scenario = scenario)
    for (name in outcomes) out[[name]] <- data[[name]][2] - data[[name]][1]
    out
  }))
  names(delta_wide)[-1] <- paste0("delta_2050_minus_2021.", outcomes)
  source_columns <- c("scenario", paste0("delta_2050_minus_2021.", source_outcomes))
  write.csv(delta_wide[, source_columns], file.path(source_dir, "delta_ssp_source_data.csv"), row.names = FALSE)

  plot_data <- do.call(rbind, lapply(seq_along(outcomes), function(i) {
    data.frame(scenario = delta_wide$scenario, policy = labels[i], delta = delta_wide[[i + 1]])
  }))
  ssp_labels <- c(SSP1 = "SSP1-2.6", SSP2 = "SSP2-4.5", SSP3 = "SSP3-7.0", SSP5 = "SSP5-8.5")
  plot_data$scenario <- factor(plot_data$scenario, levels = names(ssp_labels), labels = unname(ssp_labels))
  plot_data$policy <- factor(plot_data$policy, levels = labels)
  policy_colours <- c("HAMR-anti" = "#5B6472", "HAMR-current" = "#277DA1",
                      "HAMR-strength 5y" = "#149C7E", "HAMR-strength 10y" = "#D9852B")
  dodge <- ggplot2::position_dodge(width = 0.86)
  y_limit <- max(ceiling((max(abs(plot_data$delta)) + 4) / 5) * 5, 25)

  figure <- ggplot2::ggplot(plot_data, ggplot2::aes(x = scenario, y = delta, fill = policy, group = policy)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.45, colour = "#333333") +
    ggplot2::geom_col(position = dodge, width = 0.86, colour = "black", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.1f", delta),
      y = ifelse(delta >= 0, delta + 0.9, delta - 0.9)), position = dodge, size = 2.5, colour = "#333333") +
    ggplot2::scale_fill_manual(values = policy_colours, limits = labels, breaks = labels, name = "Policy strength") +
    ggplot2::scale_y_continuous(limits = c(-y_limit, y_limit), breaks = seq(-y_limit, y_limit, by = 10),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::labs(x = "SSP pathway",
      y = "Change in European mean HAMR (2050 - 2021)") +
    ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.title = ggplot2::element_text(size = 8), axis.text = ggplot2::element_text(size = 7, colour = "black"),
      legend.position = "bottom", legend.justification = "center", legend.title = ggplot2::element_text(size = 7.2, face = "bold"),
      legend.text = ggplot2::element_text(size = 6.8), legend.key.width = grid::unit(13, "pt"),
      legend.key.height = grid::unit(9, "pt"), plot.title = ggplot2::element_text(size = 9.5, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 7.2, colour = "#4D4D4D", hjust = 0),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 6.5, "pt")
    )

  grDevices::cairo_pdf(file.path(figure_dir, "delta_ssp.pdf"), width = 183 / 25.4,
                        height = 110 / 25.4, family = "Arial", onefile = TRUE)
  print(figure)
  grDevices::dev.off()
  cat("Created output/figures/delta_ssp.pdf\n")
  figures <- c("policyColour_SSPshape_CI.pdf" = "Figure_4a.pdf", "delta_ssp.pdf" = "Figure_4b.pdf")
  target <- file.path(project, "output")
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  for (name in names(figures)) {
    if (!file.copy(file.path(root, "output", "figures", name), file.path(target, figures[[name]]), overwrite = TRUE)) {
      stop("Could not copy figure: ", name)
    }
  }
  cat("Created Figure_4a.pdf and Figure_4b.pdf in output/.\n")
}
run_figure()
