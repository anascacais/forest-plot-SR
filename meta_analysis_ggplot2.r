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
dat <- dat[order(dat$Year, decreasing = TRUE), ]
subgroups <- sort(unique(dat$subgroup), decreasing = TRUE)
subgroup_counts <- table(dat$subgroup)

res_dict <- setNames(lapply(subgroups, function(x) list(count = subgroup_counts[[x]])), subgroups)
dat <- transform(dat, ci_low = yi - 1.96 * sqrt(vi), ci_upp = yi + 1.96 * sqrt(vi))

res <- rma(yi, vi, data = dat)
weights <- fmtx(weights(res), digits = 2)
dat <- mutate(dat, weights = weights)


############################################################################
get_summary <- function(dat, y) {
    ### fit random-effects model
    res <- rma(dat$yi, dat$vi, data = dat)

    summary_polygon <- data.frame(
        x = c(res$ci.ub, coef(res), res$ci.lb, coef(res)),
        y = c(y, y + 0.25, y, y - 0.25)
    )

    return(list(summary_polygon, res))
}

############################################################################
# get forest plot

layers <- list()
y_iter <- 2

entries <- data.frame(
    entry = numeric(),
    label = character(),
    weight = character(),
    metric = character(),
    stringsAsFactors = FALSE
)

# Entry for titles (last entry) + empty entry after each subgroup summary
n_entries <- length((dat$yi)) + length(subgroups) * 2 + 2

for (sg in rev(names(res_dict))) {
    # get subset of data for this subgroup and set the y-axis value as a sequence
    dat_subgroup <- subset(dat, subgroup == sg)
    dat_subgroup <- mutate(dat_subgroup, yn = seq(y_iter + 1, y_iter + res_dict[[sg]]$count))

    # get summary for this subgroup
    output <- get_summary(dat_subgroup, y_iter)
    summary_polygon <- output[[1]]
    res_aux <- output[[2]]

    # add subgroup label
    new_entries <- data.frame(
        entry = n_entries - y_iter - res_dict[[sg]]$count + 0.5,
        label = sg,
        weight = "",
        metric = ""
    )
    entries <- rbind(entries, new_entries)

    # create label between the entry and the Study and respective weight
    new_entries <- data.frame(
        entry = seq(n_entries - y_iter, n_entries - y_iter - res_dict[[sg]]$count + 1, by = -1),
        label = paste(dat_subgroup$Authors, dat_subgroup$Year, sep = ", "),
        weight = dat_subgroup$weights,
        metric = paste(fmtx(dat_subgroup$yi, digits = 2), " (", fmtx(dat_subgroup$ci_low, digits = 2), "-", fmtx(dat_subgroup$ci_upp, digits = 2), ")", sep = "")
    )
    entries <- rbind(entries, new_entries)

    new_entries <- data.frame(
        entry = n_entries - y_iter + 1,
        label =
            paste(
                "Subtotal (I^2 =",
                round(res_aux$I2, digits = 3), "%, ",
                fmtp(res_aux$QMp, digits = 2, pname = "p", add0 = TRUE, sep = TRUE, equal = TRUE),
                ")"
            ),
        weight = sum(as.double(dat_subgroup$weights)),
        metric = paste(fmtx(coef(res_aux), digits = 2), " (", fmtx(res_aux$ci.lb, digits = 2), "-", fmtx(res_aux$ci.ub, digits = 2), ")", sep = "")
    )
    entries <- rbind(entries, new_entries)

    # Create layers for this subgroup, adding the mean metric, CI and summary polygon

    subgroup_layers <- list(
        geom_polygon(
            data = summary_polygon, aes(x = x, y = y),
            # fill = "gray", alpha = 0.3,
            color = "black",
            fill = "transparent",
            inherit.aes = FALSE
        ),
        geom_point(
            data = dat_subgroup, aes(x = yi, y = yn), shape = 15, size = 3
        ),
        geom_linerange(data = dat_subgroup, aes(xmin = ci_low, xmax = ci_upp, y = yn))
    )

    layers <- c(layers, subgroup_layers)

    # add empty entry between subgroups
    y_iter <- y_iter + res_dict[[sg]]$count + 2
}


# Add summary for all studies
output <- get_summary(dat, 1)
summary_polygon <- output[[1]]
res <- output[[2]]
summary_layer <- list(
    geom_polygon(
        data = summary_polygon, aes(x = x, y = y),
        # fill = "gray", alpha = 0.3,
        color = "black",
        fill = "transparent",
        inherit.aes = FALSE
    )
)
layers <- c(layers, summary_layer)

new_entries <- data.frame(
    entry = n_entries,
    label = paste(
        "Overall (I^2 =",
        round(res$I2, digits = 3), "%, ",
        fmtp(res$QMp, digits = 2, pname = "p", add0 = TRUE, sep = TRUE, equal = TRUE),
        ")"
    ),
    weight = "100.00",
    metric = paste(fmtx(coef(res), digits = 2), " (", fmtx(res$ci.lb, digits = 2), "-", fmtx(res$ci.ub, digits = 2), ")", sep = "")
)
entries <- rbind(new_entries, entries)
# entries[entries$entry == n_entries, ] <- new_entries


# Combine all layers into a single plot
forest_plot <- ggplot() +
    theme_classic() +
    layers

# rename axes
forest_plot <- forest_plot +
    labs(x = "AUC", y = "")


# Add effect line
line_data <- data.frame(x = c(coef(res), coef(res)), y = c(0, 23.5))
forest_plot <- forest_plot +
    # geom_hline(yintercept = 1.5) +
    # geom_hline(yintercept = 10.5) +
    # geom_hline(yintercept = 17.5) +
    geom_hline(yintercept = 23.5) +
    geom_line(data = line_data, aes(x = x, y = y), linetype = "dashed")




forest_plot <- forest_plot +
    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
    ) +
    coord_cartesian(ylim = c(1, n_entries))

# plot(forest_plot)

new_entry <- data.frame(
    entry = 1,
    label = "Study",
    weight = "Weight (%)",
    metric = "AUC (95% CI)"
)
entries <- rbind(new_entry, entries)
entries <- entries[order(entries$entry, decreasing = FALSE), ]

############################################################################
# get y-axis for the forest plot

p_labels <-
    entries |>
    ggplot(aes(y = rev(entry)))
p_labels <- entries |>
    ggplot(aes(y = seq(max(entry), min(entry), by = -1)))

# add the study as text (instead of as a label on the y axis)
p_labels <-
    p_labels +
    # geom_hline(yintercept = 1.5) +
    # geom_hline(yintercept = 10.5) +
    # geom_hline(yintercept = 17.5) +
    geom_hline(yintercept = 23.5) +
    geom_text(aes(x = 0, label = label),
        hjust = 0,
        fontface = ifelse(grepl("Subtotal", entries$label) | grepl("Study", entries$label) | grepl("Overall", entries$label) | (entries$label %in% subgroups), "bold", "plain"),
    )

# p_labels <- ggplot(entries, aes(y = rev(entry))) +
#     # geom_point() + # Add points (just as an example)
#     geom_text(aes(x = 1, label = paste(label, entry)), y = rev(entries$entry), vjust = -0.5) # Add text annotations

# remove the background and edit the sizing so that this left size of the plot will match up neatly with the middle and right sides of the plot
p_labels <-
    p_labels +
    theme_void() +
    coord_cartesian(xlim = c(0, 4.5), ylim = c(1, n_entries))

plot(p_labels)


############################################################################
# get extra annotations


# p_annot <- entries |>
#     ggplot(aes(y = rev(entry)))

p_annot <- entries |>
    ggplot(aes(y = seq(max(entry), min(entry), by = -1)))


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
    # geom_hline(yintercept = 1.5) +
    # geom_hline(yintercept = 10.5) +
    # geom_hline(yintercept = 17.5) +
    geom_hline(yintercept = 23.5) +
    theme_void() +
    coord_cartesian(xlim = c(0, 2), ylim = c(1, n_entries))

# plot(p_annot)

# ############################################################################

layout <- c(
    area(t = 0, l = 0, b = 30, r = 4), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
    area(t = 0, l = 5, b = 30, r = 9), # starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
    area(t = 0, l = 10, b = 30, r = 12) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)

p <- p_labels + forest_plot + p_annot + plot_layout(design = layout)
plot(p)
ggsave("plot.pdf", plot = p)
