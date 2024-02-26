# Load libraries for meta-analysis
library(robumeta)
library(metafor)
library(dplyr)

# Load library for reading excel files
library(readxl)
library(ggplot2)

library(tidyverse)
library(gt)
library(patchwork)


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

dat <- mutate(dat, key = paste(
    dat$Authors, dat$Year, dat$"Forecast horizon", dat$"Input data", dat$"Training and testing approach",
    sep = ", "
))
dat$key <- as.character(dat$key)

### order data according to subgroup and yi
dat <- dat[order(dat$subgroup, dat$Year), decreasing = FALSE]
subgroups <- sort(unique(dat$subgroup))
subgroup_counts <- table(dat$subgroup)
subgroup_counts_dict <- as.list(names = names(subgroup_counts), subgroup_counts)

res_dict <- setNames(lapply(subgroups, function(x) list(count = subgroup_counts[[x]])), subgroups)
dat <- transform(dat, ci_low = yi - 1.96 * sqrt(vi), ci_upp = yi + 1.96 * sqrt(vi))

res <- rma(yi, vi, data = dat)
weights <- fmtx(weights(res), digits = 2)
dat <- mutate(dat, weights = weights)

############################################################################
get_summary <- function(dat, y) {
    ### fit random-effects model
    res <- rma(dat$yi, dat$vi, data = dat)
    ci_midpoint <- (res$ci.ub - res$ci.lb) / 2 + res$ci.lb


    summary_polygon <- data.frame(
        x = c(res$ci.ub, ci_midpoint, res$ci.lb, ci_midpoint),
        y = c(y - 1, y - 1 + 0.5, y - 1, y - 1 - 0.5)
    )

    return(list(summary_polygon, res))
}

############################################################################
# get forest plot

layers <- list()
y_iter <- 1

entries <- data.frame(
    entry = numeric(),
    label = character(),
    weight = character(),
    metric = character(),
    stringsAsFactors = FALSE
)

for (sg in names(res_dict)) {
    # get subset of data for this subgroup and set the y-axis value as a sequence
    dat_subgroup <- subset(dat, subgroup == sg)
    dat_subgroup <- mutate(dat_subgroup, yn = seq(y_iter, y_iter + res_dict[[sg]]$count - 1))

    # get summary for this subgroup
    output <- get_summary(dat_subgroup, y_iter)
    summary_polygon <- output[[1]]
    res_aux <- output[[2]]

    # create label between the entry and the Study and respective weight
    new_entries <- data.frame(
        entry = dat_subgroup$yn,
        label = paste(dat_subgroup$Authors, dat_subgroup$Year, sep = ", "),
        weight = dat_subgroup$weights,
        metric = paste(fmtx(dat_subgroup$yi, digits = 2), " (", fmtx(dat_subgroup$ci_low, digits = 2), "-", fmtx(dat_subgroup$ci_upp, digits = 2), ")", sep = "")
    )
    entries <- rbind(entries, new_entries)

    new_entries <- data.frame(
        entry = c(y_iter + res_dict[[sg]]$count, y_iter + res_dict[[sg]]$count + 1),
        label = c(
            paste(
                "Subtotal", sg, "(I^2 =",
                round(res_aux$I2, digits = 3), "%, ",
                fmtp(res_aux$QMp, digits = 2, pname = "p", add0 = TRUE, sep = TRUE, equal = TRUE),
                ")"
            ), ""
        ),
        weight = c(sum(as.double(dat_subgroup$weights)), ""),
        metric = c("", "")
    )
    entries <- rbind(entries, new_entries)

    # Create layers for this subgroup, adding the mean metric, CI and summary polygon
    subgroup_layers <- list(
        geom_point(
            data = dat_subgroup, aes(x = yi, y = yn), shape = 15, size = 3
        ),
        geom_linerange(data = dat_subgroup, aes(xmin = ci_low, xmax = ci_upp, y = yn)),
        geom_polygon(
            data = summary_polygon, aes(x = x, y = y),
            fill = "gray", alpha = 0.3, inherit.aes = FALSE
        )
    )

    layers <- c(layers, subgroup_layers)
    # add empty entry between subgroups
    y_iter <- y_iter + res_dict[[sg]]$count + 2
}


# Combine all layers into a single plot
forest_plot <- ggplot() +
    theme_classic() +
    layers

# rename axes
forest_plot <- forest_plot +
    labs(x = "AUC", y = "")


forest_plot <- forest_plot +
    geom_vline(xintercept = coef(res), linetype = "dashed")


forest_plot <- forest_plot +
    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
    )
plot(forest_plot)

############################################################################
# get y-axis for the forest plot

entries$entry <- entries$entry + 1
new_entry <- data.frame(
    entry = 1,
    label = "Study",
    weight = "Weight (%)",
    metric = "AUC (95% CI)"
)
entries <- rbind(new_entry, entries)

p_labels <-
    entries |>
    ggplot(aes(y = rev(entry)))

# add the study as text (instead of as a label on the y axis)
p_labels <-
    p_labels +
    geom_text(aes(x = 0, label = label),
        hjust = 0,
        fontface = ifelse(grepl("Subtotal", entries$label) | grepl("Study", entries$label) | grepl("Overall", entries$label), "bold", "plain")
    )

# remove the background and edit the sizing so that this left size of the plot will match up neatly with the middle and right sides of the plot
p_labels <-
    p_labels +
    theme_void() +
    coord_cartesian(xlim = c(0, 4.5))

# plot(p_labels)


############################################################################
# get extra annotations

p_annot <- entries |>
    ggplot(aes(y = rev(entry)))

p_annot <- p_annot +
    geom_text(
        aes(x = 0, label = weight),
        hjust = 0,
        fontface = ifelse(entries$weight == "Weight (%)", "bold", "plain")
    )

p_annot <- p_annot +
    geom_text(
        aes(x = 1, label = metric),
        hjust = 0,
        fontface = ifelse(grepl("(95% CI)", entries$metric), "bold", "plain")
    )

# remove the background and edit the sizing so that this left size of the plot will match up neatly with the middle and right sides of the plot
p_annot <-
    p_annot +
    theme_void() +
    coord_cartesian(xlim = c(0, 2.5))

# plot(p_annot)

# ############################################################################

# layout <- c(
#     area(t = 0, l = 0, b = 30, r = 3), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
#     area(t = 1, l = 4, b = 30, r = 9), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
#     area(t = 0, l = 9, b = 30, r = 11) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
# )

p <- p_labels + forest_plot + p_annot # + plot_layout(design = layout)
# plot(p)
