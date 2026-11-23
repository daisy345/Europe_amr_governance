# Publication-style trend plot for Ascore (no Joinpoint selected)
# Input: input/animal AMR score.txt

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
input_file <- file.path(script_dir, "..", "input", "animal AMR score.txt")
output_pdf <- file.path(script_dir, "..", "output", "Animal_AMR_score.pdf")
output_png <- file.path(script_dir, "..", "output", "Animal_AMR_score.png")
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
names(dat)[match("t", names(dat))] <- "year"
required <- c("year", "amr", "Model", "Joinpoints")
if (!all(required %in% names(dat))) {
  stop("The Joinpoint data export must contain: ", paste(required, collapse = ", "))
}
dat <- dat[order(dat$year), ]
if (!all(dat$Joinpoints == 0)) stop("This script is intended for the zero-Joinpoint Ascore model.")

# APC = 0.5103% (95% CI 0.2830 to 0.7256), P < 0.001.
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
    xlab = "Year", ylab = "Ascore", xaxt = "n", yaxt = "n", bty = "l"
  )
  axis(1, at = dat$year, labels = ifelse(dat$year %% 2 == 1, dat$year, ""), cex.axis = 0.86)
  axis(2, at = pretty(ylim, n = 5), cex.axis = 0.9)
  lines(dat$year, dat$amr, col = "grey65", lwd = 1.1)
  points(dat$year, dat$amr, pch = 21, bg = "white", col = "grey35", cex = 1.05, lwd = 1)
  lines(dat$year, dat$Model, col = "#1F4E79", lwd = 2.3)
  legend(
    "bottomleft", inset = 0.018, bty = "n", cex = 0.88, bg = "white",
    legend = c("Observed Ascore", "Fitted linear trend"),
    col = c("grey45", "#1F4E79"), lwd = c(1.1, 2.3),
    pch = c(21, NA), pt.bg = c("white", NA), pt.cex = c(1, NA)
  )
  text(2010.1, y_span[2] - 0.10, labels = "APC = 0.51%*", adj = c(0, 0),
       cex = 0.84, col = "#1F4E79")
  mtext("* P < 0.05", side = 3, adj = 0, line = -0.1, cex = 0.72, col = "grey30")
}

pdf(output_pdf, width = 7.1, height = 5.0, useDingbats = FALSE)
draw_plot()
dev.off()

png(output_png, width = 7.1, height = 5.0, units = "in", res = 600)
draw_plot()
dev.off()
