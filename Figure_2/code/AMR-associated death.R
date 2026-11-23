# Publication-style Joinpoint plot for deaths
# Input: input/AMR-associated death.txt

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
input_file <- file.path(script_dir, "..", "input", "AMR-associated death.txt")
output_pdf <- file.path(script_dir, "..", "output", "AMR-associated_death.pdf")
output_png <- file.path(script_dir, "..", "output", "AMR-associated_death.png")
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("year", "death_mean", "Model", "Flag")
if (!all(required %in% names(dat))) {
  stop("The Joinpoint data export must contain: ", paste(required, collapse = ", "))
}
dat <- dat[order(dat$year), ]
joinpoint_year <- dat$year[which(grepl("Joinpoint", dat$Flag, fixed = TRUE))]
if (length(joinpoint_year) != 1L || joinpoint_year != 2018) {
  stop("Expected one Joinpoint at 2018 in this death export.")
}

# 2007-2018: APC = 0.0445%, P = 0.407; 2018-2021: APC = -2.9653%, P < 0.001.
draw_plot <- function() {
  y_span <- range(c(dat$death_mean, dat$Model), na.rm = TRUE)
  y_pad <- diff(y_span) * 0.14
  ylim <- c(y_span[1] - y_pad * 0.35, y_span[2] + y_pad)

  par(
    mar = c(4.3, 5.2, 0.7, 0.7), las = 1,
    family = "sans", tcl = -0.25, mgp = c(3.0, 0.7, 0)
  )
  plot(
    dat$year, dat$death_mean, type = "n", xlim = range(dat$year), ylim = ylim,
    xlab = "Year", ylab = "Deaths", xaxt = "n", yaxt = "n", bty = "l"
  )
  axis(1, at = dat$year, labels = ifelse(dat$year %% 2 == 1, dat$year, ""), cex.axis = 0.84)
  y_ticks <- pretty(ylim, n = 5)
  axis(2, at = y_ticks, labels = format(y_ticks, big.mark = ",", trim = TRUE), cex.axis = 0.9)
  abline(v = joinpoint_year, lty = 2, lwd = 0.9, col = "grey55")
  lines(dat$year, dat$death_mean, col = "grey65", lwd = 1.1)
  points(dat$year, dat$death_mean, pch = 21, bg = "white", col = "grey35", cex = 1.05, lwd = 1)
  lines(dat$year, dat$Model, col = "#1F4E79", lwd = 2.3)
  points(joinpoint_year, dat$Model[dat$year == joinpoint_year], pch = 21,
         bg = "#1F4E79", col = "#1F4E79", cex = 1.15)
  legend(
    "bottomleft", inset = 0.018, bty = "n", cex = 0.86, bg = "white",
    legend = c("Observed deaths", "Joinpoint fitted trend"),
    col = c("grey45", "#1F4E79"), lwd = c(1.1, 2.3),
    pch = c(21, NA), pt.bg = c("white", NA), pt.cex = c(1, NA)
  )
  text(joinpoint_year, ylim[2] - 0.03 * diff(ylim), labels = "Joinpoint: 2018",
       pos = 4, offset = 0.30, cex = 0.73, col = "grey30")
  text(2007.2, 10270, labels = "APC = 0.04%", adj = c(0, 0), cex = 0.82, col = "#1F4E79")
  text(2018.15, 9630, labels = "APC = -2.97%*", adj = c(0, 0), cex = 0.82, col = "#1F4E79")
  mtext("* P < 0.05", side = 3, adj = 0, line = -0.1, cex = 0.72, col = "grey30")
}

pdf(output_pdf, width = 7.1, height = 5.0, useDingbats = FALSE)
draw_plot()
dev.off()

png(output_png, width = 7.1, height = 5.0, units = "in", res = 600)
draw_plot()
dev.off()
