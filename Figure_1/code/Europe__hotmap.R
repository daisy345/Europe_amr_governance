# Recreate the four European heat maps from the CSV tables and bundled boundary data.
# Run with: Rscript code/Europe__hotmap.R
# The script writes PDFs to the sibling output/ directory.

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this file with Rscript so its location can be resolved.")
}
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "input")
output_dir <- file.path(project_dir, "output")
map_file <- file.path(input_dir, "geodata", "Europe_Cleaned.shp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(map_file)) stop("Missing bundled boundary file: ", map_file)

europe <- sf::st_read(map_file, quiet = TRUE)
if (!"CNTR_ID" %in% names(europe)) stop("Boundary data must contain CNTR_ID.")

maps <- list(
  list(file = "Human_AMR_score_hotmap.csv", output = "Human_score.pdf",
       legend = "Human AMR Score", limits = c(40, 100),
       breaks = seq(40, 100, by = 20), direction = 1, squish = FALSE),
  list(file = "Human_AMU_hotmap.csv", output = "Human_AMU.pdf",
       legend = "Human AMU", limits = c(0, 50),
       breaks = seq(0, 50, by = 10), direction = -1, squish = FALSE),
  list(file = "Animal_AMR_score_hotmap.csv", output = "Animal_AMR_score.pdf",
       legend = "animal AMR score", limits = c(60, 100),
       breaks = seq(60, 100, by = 20), direction = 1, squish = FALSE),
  list(file = "Animal_AMU_hotmap.csv", output = "Animal_AMU.pdf",
       legend = "animal AMU", limits = c(0, 200),
       breaks = seq(0, 200, by = 50), direction = -1, squish = TRUE)
)

for (spec in maps) {
  data_file <- file.path(input_dir, spec$file)
  if (!file.exists(data_file)) stop("Missing input table: ", data_file)
  dat <- utils::read.csv(data_file, check.names = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("", "NA"))
  required <- c("iso", "value")
  if (!all(required %in% names(dat))) {
    stop(spec$file, " must contain columns named iso and value.")
  }
  if (anyDuplicated(dat$iso)) stop(spec$file, " contains duplicated iso identifiers.")
  if (!is.numeric(dat$value)) stop(spec$file, " column value must be numeric.")
  unmatched <- setdiff(dat$iso, europe$CNTR_ID)
  if (length(unmatched)) {
    warning(spec$file, " has identifiers absent from boundary data: ",
            paste(unmatched, collapse = ", "))
  }

  map_data <- merge(europe, dat, by.x = "CNTR_ID", by.y = "iso",
                    all.x = TRUE, sort = FALSE)
  scale_args <- list(
    palette = "RdYlBu", direction = spec$direction,
    na.value = "#f0f0f0", name = spec$legend,
    limits = spec$limits, breaks = spec$breaks
  )
  if (spec$squish) scale_args$oob <- scales::squish

  p <- ggplot(map_data) +
    geom_sf(aes(fill = value), color = "black", linewidth = 0.1) +
    do.call(scale_fill_distiller, scale_args) +
    coord_sf(xlim = c(-25, 45), ylim = c(34, 72), expand = FALSE) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.key.width = grid::unit(0.8, "cm")
    )

  output_file <- file.path(output_dir, spec$output)
  ggsave(output_file, plot = p, width = 10, height = 8, device = grDevices::cairo_pdf)
  message("Wrote ", output_file)
}
