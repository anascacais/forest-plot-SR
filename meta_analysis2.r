# Load libraries for meta-analysis
library(robumeta)
library(metafor)
library(dplyr)

# Load library for reading excel files
library(readxl)


data <- read_excel(
    "/Users/anasofiacc/Library/CloudStorage/OneDrive-UniversidadedeLisboa/PhD/PreEpiSeizures/Seizure forecast/Systematic review/meta-analysis.xlsx",
    sheet = "meta-analysis",
    na = "NA"
)
data <- data[1:(nrow(data) - 1), ]

### prepare data
dat <- subset(
    data, (!is.na(data$"AUC (SD)")) & (data$"AUC (SD)" != 0)
)

dat$"AUC (SD)" <- as.numeric(dat$"AUC (SD)")
dat$"AUC (mean)" <- as.numeric(dat$"AUC (mean)")
dat$"# Patients" <- as.numeric(dat$"# Patients")

dat <- mutate(dat,
    SE = dat$"AUC (SD)" / sqrt(dat$"# Patients")
)
dat <- mutate(dat,
    subgroup = dat$"Input data"
)
dat <- mutate(dat, vi = SE^2)
dat <- mutate(dat, yi = dat$"AUC (mean)")


### order data according to subgroup and yi
dat <- dat[order(dat$subgroup, dat$Year), decreasing = FALSE]

### fit random-effects model
res <- rma(yi, vi, data = dat)
res

############################################################################

### colors to be used in the plot
colp <- "#6b58a6"
coll <- "#a7a9ac"

### total number of studies
k <- nrow(dat)
n_subgroups <- length(unique(dat$subgroup))


### generate point sizes
psize <- weights(res)
psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))

### get the weights and format them as will be used in the forest plot
weights <- fmtx(weights(res), digits = 1)

### adjust the margins
# par(mar = c(2.7, 3.2, 2.3, 1.3), mgp = c(3, 0, 0), tcl = 0.15)
par(mar = c(5, 3.2, 5, 1.3), mgp = c(3, 0, 0), tcl = 0.15)

forest <- forest(dat$yi, dat$vi,
    xlim = c(min(dat$yi) - 2, max(dat$yi) + 1),
    ylim = c(-0.5, k + 3),
    alim = c(0, 1),
    cex = 0.88,
    pch = 18,
    psize = psize,
    efac = 0,
    refline = NA,
    xlab = "",
    ilab = cbind(weights),
    ilab.xpos = c(1.2), annosym = c(" (", " to ", ")"),
    rowadj = -.07,
    slab = paste(dat$Authors, dat$Year, sep = ", "),
)

### add the vertical reference line at 0
# segments(0, -1, 0, k + 1.6, col = coll)

### add the vertical reference line at the pooled estimate
segments(coef(res), 0, coef(res), k, col = colp, lty = "33", lwd = 0.8)

### redraw the CI lines and points in the chosen color
# segments(res$ci.lb, k:1, res$ci.ub, k:1, col = colp, lwd = 1.5)
# points(dat$yi, k:1, pch = 18, cex = psize * 1.15, col = "white")
# points(dat$yi, k:1, pch = 18, cex = psize, col = colp)

### add the summary polygon
addpoly(res,
    row = 0, mlab =
        paste(
            "Overall (I^2 =",
            round(res$I2, digits = 3), "%",
            fmtp(res$QMp, digits = 2, pname = "p", add0 = TRUE, sep = TRUE, equal = TRUE),
            ")"
        ),
    efac = 2, col = colp, border = colp, font = 2
)

### fit random-effects model in the three subgroups
res.s <- rma(
    yi, vi,
    subset = (subgroup == "surrogate measures of preictal state"), data = dat
)
res.r <- rma(
    yi, vi,
    subset = (subgroup == "cyclic distribution of events"), data = dat
)
res.a <- rma(
    yi, vi,
    subset = (subgroup == "both"), data = dat
)

### add summary polygons for the thr ee subgroups
addpoly(res.s,
    row = k - 4 + 0.5, mlab = paste(
        "Subgroup (I^2 =",
        round(res.s$I2, digits = 3), "%",
        fmtp(res.s$QMp, digits = 2, pname = "p", add0 = TRUE, sep = TRUE, equal = TRUE),
        ")"
    ),
    efac = 2, col = colp, border = colp, font = 2
)

addpoly(res.s,
    row = k - 4 - 5 + 0.5, mlab = paste(
        "Subgroup (I^2 =",
        round(res.s$I2, digits = 3), "%",
        fmtp(res.s$QMp, digits = 2, pname = "p", add0 = TRUE, sep = TRUE, equal = TRUE),
        ")"
    ),
    efac = 2, col = colp, border = colp, font = 2
)

### add the horizontal line at the top
# abline(h = k + 1.6, col = coll)

### redraw the x-axis in the chosen color
axis(side = 1, at = seq(0, 1, by = 0.5), col = coll, labels = FALSE)

### now we add a bunch of text; since some of the text falls outside of the
### plot region, we set xpd=NA so nothing gets clipped
par(xpd = NA)

### adjust cex as used in the forest plot and use a bold font
par(cex = forest$cex, font = 2)

### add headings
text(forest$xlim[1], k + 1.5, pos = 4, "Study or\nsubgroup")
text(
    c(
        forest$ilab.xpos[1], forest$xlim[2] - 0.25
    ),
    k + 1.5,
    c("Weight\n(%)", "AUC (95% CI)")
)

### use a non-bold font for the rest of the text
par(cex = forest$cex, font = 1)

### add the 100.0 for the sum of the weights
text(forest$ilab.xpos[1], 0, "100.0")

### fit meta-regression model to test for subgroup differences
res_sub <- rma(yi, vi, mods = ~subgroup, data = dat)
### add text for the test of subgroup differences
text(forest$xlim[1], -0.5, pos = 4, bquote(paste(
    "Heterogeneity between groups: ",
    .(fmtp(res_sub$QMp, digits = 2, pname = "p", add0 = TRUE, sep = TRUE, equal = TRUE))
)))
