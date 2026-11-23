# Publication-style Joinpoint plot for AAMU
# Input: input/animal AMU.txt

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
input_file <- file.path(script_dir, "..", "input", "animal AMU.txt")
output_pdf <- file.path(script_dir, "..", "output", "Animal_AMU.pdf")
output_png <- file.path(script_dir, "..", "output", "Animal_AMU.png")
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("year", "mean_value", "Model", "Flag")
if (!all(required %in% names(dat))) {
  stop("The Joinpoint data export must contain: ", paste(required, collapse = ", "))
}
dat <- dat[order(dat$year), ]
joinpoint_year <- dat$year[which(grepl("Joinpoint", dat$Flag, fixed = TRUE))]
if (length(joinpoint_year) != 2L || !all(joinpoint_year == c(2015, 2020))) {
  stop("Expected Joinpoints at 2015 and 2020 in this AAMU export.")
}

# Only the 2020-2022 APC is significant (P < 0.05).
apc_note <- c("APC = -1.08%", "APC = -4.43%", "APC = -11.15%*")

draw_plot <- function() {
  y_span <- range(c(dat$mean_value, dat$Model), na.rm = TRUE)
  y_pad <- diff(y_span) * 0.14
  ylim <- c(y_span[1] - y_pad * 0.35, y_span[2] + y_pad)

  par(
    mar = c(4.3, 4.7, 0.7, 0.7), las = 1,
    family = "sans", tcl = -0.25, mgp = c(2.6, 0.7, 0)
  )
  plot(
    dat$year, dat$mean_value, type = "n", xlim = range(dat$year), ylim = ylim,
    xlab = "Year", ylab = "AAMU", xaxt = "n", yaxt = "n", bty = "l"
  )
  axis(1, at = dat$year, labels = dat$year, cex.axis = 0.82)
  axis(2, at = pretty(ylim, n = 5), cex.axis = 0.9)
  abline(v = joinpoint_year, lty = 2, lwd = 0.9, col = "grey55")
  lines(dat$year, dat$mean_value, col = "grey65", lwd = 1.1)
  points(dat$year, dat$mean_value, pch = 21, bg = "white", col = "grey35", cex = 1.05, lwd = 1)
  lines(dat$year, dat$Model, col = "#1F4E79", lwd = 2.3)
  points(joinpoint_year, dat$Model[match(joinpoint_year, dat$year)], pch = 21,
         bg = "#1F4E79", col = "#1F4E79", cex = 1.15)
  legend(
    "bottomleft", inset = 0.018, bty = "n", cex = 0.86, bg = "white",
    legend = c("Observed AAMU", "Joinpoint fitted trend"),
    col = c("grey45", "#1F4E79"), lwd = c(1.1, 2.3),
    pch = c(21, NA), pt.bg = c("white", NA), pt.cex = c(1, NA)
  )
  text(2015, ylim[2] - 0.03 * diff(ylim), labels = "Joinpoint: 2015",
       pos = 4, offset = 0.28, cex = 0.72, col = "grey30")
  text(2020, ylim[2] - 0.10 * diff(ylim), labels = "Joinpoint: 2020",
       pos = 4, offset = 0.28, cex = 0.72, col = "grey30")
  text(2010.2, 113.5, labels = apc_note[1], adj = c(0, 0), cex = 0.80, col = "#1F4E79")
  text(2015.2, 96.5, labels = apc_note[2], adj = c(0, 0), cex = 0.80, col = "#1F4E79")
  text(2020.2, 74.5, labels = apc_note[3], adj = c(0, 0), cex = 0.80, col = "#1F4E79")
  mtext("* P < 0.05", side = 3, adj = 0, line = -0.1, cex = 0.72, col = "grey30")
}

pdf(output_pdf, width = 7.1, height = 5.0, useDingbats = FALSE)
draw_plot()
dev.off()

png(output_png, width = 7.1, height = 5.0, units = "in", res = 600)
draw_plot()
dev.off()
