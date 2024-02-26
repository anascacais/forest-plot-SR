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
weights <- weights(res)
print(weights)


### set up forest plot (with 2x2 table counts added; the 'rows' argument is
### used to specify in which rows the outcomes will be plotted)
xref_min <- -2
xref_max <- 2


forest(
    res,
    # xlim = c(
    #     min(par("usr")[1], par("usr")[2]) * 0.1,
    #     max(1.1, max(dat$yi) + 0.2)
    # ),
    xlim = c(xref_min, xref_max),
    at = seq(min(0, min(dat$yi) + 0.1), max(1, max(dat$yi) + 0.1), by = 0.5),
    cex = 1.00,
    ylim = c(-1, 25),
    # order = order(dat$subgroup, decreasing = FALSE),
    rows = c(2:8, 11:15, 18:21),
    mlab = mlabfun("RE Model for All Studies", res),
    header = "Author(s) and Year",
    slab = paste(Authors, Year, sep = ", ")
)

y_positions <- seq_along(weights)
text(
    x = max(par("usr")[1], par("usr")[2]) * 1.1,
    y = y_positions,
    labels = round(weights, 2),
    pos = 4
)


### Add line of "no-effect"
abline(v = 0.5, lty = 2) # , col = "red")

### set font expansion factor (as in forest() above) and use a bold font
op <- par(cex = 1.00, font = 2)

### switch to bold italic font
par(font = 4)
### add text for the subgroups
text(
    x = xref_min,
    y = c(22, 16, 9),
    labels = c(
        "Surrogate Measures of Preictal State",
        "Cyclic Distribution of Events",
        "Both"
    ), pos = 4
)


### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.s <- rma(yi, vi, subset = (subgroup == "surrogate measures of preictal state"), data = dat)
res.r <- rma(yi, vi, subset = (subgroup == "cyclic distribution of events"), data = dat)
res.a <- rma(yi, vi, subset = (subgroup == "both"), data = dat)

### add summary polygons for the three subgroups
addpoly(res.s, row = 17, mlab = mlabfun("RE Model for Subgroup", res.s))
addpoly(res.r, row = 10, mlab = mlabfun("RE Model for Subgroup", res.r))
addpoly(res.a, row = 1, mlab = mlabfun("RE Model for Subgroup", res.a))

### fit meta-regression model to test for subgroup differences
res <- rma(yi, vi, mods = ~subgroup, data = dat)

### add text for the test of subgroup differences
text(
    xref_min,
    -1.5,
    pos = 4, cex = 1.00,
    bquote(paste(
        "Test for Subgroup Differences: ",
        Q[M], " = ", .(fmtx(res$QM, digits = 2)),
        ", df = ", .(res$p - 1), ", ",
        .(fmtp(res$QMp, digits = 2, pname = "p", add0 = TRUE, sep = TRUE, equal = TRUE))
    ))
)
