# install.packages(c("robumeta", "metafor", "dplyr", "readxl", "vscDebugger"))

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


data2analyse <- subset(
    data, (!is.na(data$"AUC (SD)")) & (data$"AUC (SD)" != 0)
)

data2analyse$"AUC (SD)" <- as.numeric(data2analyse$"AUC (SD)")
data2analyse$"AUC (mean)" <- as.numeric(data2analyse$"AUC (mean)")
data2analyse$"# Patients" <- as.numeric(data2analyse$"# Patients")

data2analyse <- mutate(data2analyse,
    SE = data2analyse$"AUC (SD)" / sqrt(data2analyse$"# Patients")
)
data2analyse <- mutate(data2analyse,
    subgroup = data2analyse$"Input data"
)
data2analyse <- mutate(data2analyse, var = SE^2)


res1 <- rma(
    data2analyse$"AUC (mean)",
    var,
    data = data2analyse,
    subset = subgroup == "surrogate measures of preictal state"
)
res2 <- rma(
    data2analyse$"AUC (mean)",
    var,
    data = data2analyse,
    subset = subgroup == "cyclic distribution of events"
)

res <- rma(data2analyse$"AUC (mean)", var, data = data2analyse, mods = ~subgroup)
print(res)
forest(res)


# res.ee <- rma(data2analyse$"AUC (mean)", var, data = data2analyse, method = "EE", mods = ~ data2analyse$"Input data" - 1)
# res.re <- rma(data2analyse$"AUC (mean)", var, data = data2analyse, mods = ~ data2analyse$"Input data" - 1)

# w.ee.re <- cbind(
#     paste0(formatC(weights(res.ee), format = "f", digits = 1, width = 4), "%"),
#     paste0(formatC(weights(res.re), format = "f", digits = 1, width = 4), "%")
# )

# forest(data2analyse$"AUC (mean)", data2analyse$var,
#     xlim = c(-11, 5), ylim = c(-2.5, 16), header = TRUE, atransf = exp,
#     at = log(c(1 / 16, 1 / 4, 1, 4, 8)), digits = c(2L, 4L), ilab = w.ee.re, ilab.xpos = c(-6, -4)
# )

# abline(h = 0)
# addpoly(res.ee, row = -1)
# addpoly(res.re, row = -2)
# text(-6, 15, "EE Model", font = 2)
# text(-4, 15, "RE Model", font = 2)
# text(-5, 16, "Weights", font = 2)
# segments(-7, 15.5, -3, 15.5)
