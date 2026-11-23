# Publication-style Joinpoint plot for H-score
# Input: input/human AMR score.txt

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
input_file <- file.path(script_dir, "..", "input", "human AMR score.txt")
output_pdf <- file.path(script_dir, "..", "output", "Human_AMR_score.pdf")
output_png <- file.path(script_dir, "..", "output", "human_AMR_score.png")
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
names(dat)[match("t", names(dat))] <- "year"

required <- c("year", "amr", "Model", "Flag")
if (!all(required %in% names(dat))) {
  stop("The Joinpoint data export must contain: ", paste(required, collapse = ", "))
}
dat <- dat[order(dat$year), ]
joinpoint_year <- dat$year[which(grepl("Joinpoint", dat$Flag, fixed = TRUE))]
if (length(joinpoint_year) != 1L) stop("Expected exactly one joinpoint in the data export.")

# Compact in-panel statistics; full confidence intervals belong in the figure caption.
apc_note <- c("APC = -0.99%*", "APC = -0.00%")

draw_plot <- function() {
  y_span <- range(c(dat$amr, dat$Model), na.rm = TRUE)
  y_pad <- diff(y_span) * 0.14
  ylim <- c(y_span[1] - y_pad * 0.35, y_span[2] + y_pad)

  par(
    mar = c(4.3, 4.5, 0.7, 0.7), las = 1,
    family = "sans", tcl = -0.25, mgp = c(2.5, 0.7, 0)
  )
  plot(
    dat$year, dat$amr, type = "n", xlim = range(dat$year), ylim = ylim,
    xlab = "Year", ylab = "H-score", xaxt = "n", yaxt = "n", bty = "l"
  )
  axis(1, at = dat$year, labels = dat$year, cex.axis = 0.82)
  axis(2, at = pretty(ylim, n = 5), cex.axis = 0.9)
  abline(v = joinpoint_year, lty = 2, lwd = 0.9, col = "grey55")
  lines(dat$year, dat$amr, col = "grey65", lwd = 1.1)
  points(dat$year, dat$amr, pch = 21, bg = "white", col = "grey35", cex = 1.05, lwd = 1)
  lines(dat$year, dat$Model, col = "#1F4E79", lwd = 2.3)
  points(joinpoint_year, dat$Model[dat$year == joinpoint_year], pch = 21,
         bg = "#1F4E79", col = "#1F4E79", cex = 1.15)

  legend(
    "bottomleft", inset = 0.018, bty = "n", cex = 0.88,
    legend = c("Observed H-score", "Joinpoint fitted trend"),
    col = c("grey45", "#1F4E79"), lwd = c(1.1, 2.3),
    pch = c(21, NA), pt.bg = c("white", NA), pt.cex = c(1, NA)
  )
  text(joinpoint_year, ylim[2] - 0.02 * diff(ylim), labels = "Joinpoint: 2014",
       pos = 4, offset = 0.35, cex = 0.78, col = "grey30")
  text(2007.2, y_span[2] - 0.30, labels = apc_note[1], adj = c(0, 0),
       cex = 0.82, col = "#1F4E79")
  text(2015.0, y_span[1] + 0.42, labels = apc_note[2], adj = c(0, 0),
       cex = 0.82, col = "#1F4E79")
  mtext("* P < 0.05", side = 3, adj = 0, line = -0.1, cex = 0.72, col = "grey30")
}

pdf(output_pdf, width = 7.1, height = 5.0, useDingbats = FALSE)
draw_plot()
dev.off()

png(output_png, width = 7.1, height = 5.0, units = "in", res = 600)
draw_plot()
dev.off()
