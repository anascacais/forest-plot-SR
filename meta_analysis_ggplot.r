library(readxl)
library(patchwork)

# Local
source("meta_analysis_aux.r")
source("data_aux.r")
source("graphic_aux.r")
source("filename.r")

subset_strategy <- NA # "Output approach"
subset_name <- NA

subgroup_strategy <- "Type of input data"

study_label_strategy <- "First author, year"

effect_size <- "AUC"


data <- read_excel(
    filename,
    sheet = "meta-analysis",
    na = "NA"
)

effect_size_mean <- paste(effect_size, "(mean)", sep = " ")
effect_size_sd <- paste(effect_size, "(SD)", sep = " ")


dat <- prepare_dataframe(
    data,
    subset_strategy,
    subset_name,
    subgroup_strategy,
    effect_size_mean,
    effect_size_sd
)


dat <- get_effect_size(dat, effect_size_mean, effect_size_sd)


# Get overall meta-analysis result
overall_res <- rma.uni(yi, vi, data = dat, method = "REML")
test_heterogeneity(data = dat, overall_res, mods = c("Output approach", "Input data source", "Forecast horizon", "Training and testing approach", "Year", "Data source"))

weights <- fmtx(weights(overall_res), digits = 2)
dat <- mutate(dat, weights = weights)

temp <- prepare_subgroup_analysis(dat, subgroup_strategy)
dat <- temp[[1]]
res_dict <- temp[[2]]
n_entries <- temp[[3]]

subgroup_analysis <- TRUE # !any(sapply(res_dict, function(x) x$count < 2))
if (!subgroup_analysis) {
    n_entries <- n_entries - length(names(res_dict))
}

############################################################################
# Prepare forest plot and entries for remaining plot info

layers <- list()
y_iter <- 1 + (1 * subgroup_analysis)

entries <- data.frame(
    entry = numeric(),
    label = character(),
    weight = character(),
    metric = character(),
    horizon = character(),
    size = character(),
    approach = character(),
    stringsAsFactors = FALSE
)


# Get soubgroup results
for (sg in rev(names(res_dict))) {
    # get subset of data for this subgroup and set the y-axis value as a sequence
    dat_subgroup <- subset(dat, subgroup == sg)
    dat_subgroup <- mutate(dat_subgroup, yn = seq(y_iter + 1, y_iter + res_dict[[sg]]$count))

    # add subgroup label
    new_entries <- data.frame(
        entry = n_entries - y_iter - res_dict[[sg]]$count + 0.5,
        label = sg
    )
    entries <- rbind.fill(entries, new_entries)

    # create label between the entry and the Study and respective weight
    new_entries <- data.frame(
        entry = seq(n_entries - y_iter, n_entries - y_iter - res_dict[[sg]]$count + 1, by = -1),
        label = dat_subgroup[[study_label_strategy]],
        weight = dat_subgroup$weights,
        metric = paste(fmtx(dat_subgroup$yi, digits = 2), " (", fmtx(dat_subgroup$ci_low, digits = 2), "-", fmtx(dat_subgroup$ci_upp, digits = 2), ")", sep = ""),
        horizon = dat_subgroup$"Forecast horizon",
        size = dat_subgroup$"# Patients",
        approach = ifelse(dat_subgroup$"Training and testing approach" == "prospective", "P", "R")
    )
    entries <- rbind.fill(entries, new_entries)
    layers <- draw_study_results(dat_subgroup, layers)


    # Add summary for this subgroup (if there are subgroups with less than 2 entries, subgroup analysis in not performed)
    if (subgroup_analysis) {
        temp <- make_summary(dat_subgroup, entries, layers, n_entries, subgroup_name = sg, pos = y_iter)
        entries <- temp[[1]]
        layers <- temp[[2]]
    }

    # add empty entry between subgroups
    y_iter <- y_iter + res_dict[[sg]]$count + 1 + (1 * subgroup_analysis)
}


# Add summary for all studies
temp <- make_summary(dat, entries, layers, n_entries, pos = 1, subtotal = FALSE)
entries <- temp[[1]]
layers <- temp[[2]]






# Add titles
new_entry <- data.frame(
    entry = 1,
    label = "First author, year",
    weight = "Weight (%)",
    metric = paste(effect_size, "(95% CI)", sep = " "),
    horizon = "Forecast\nhorizon",
    size = "Sample\nsize",
    approach = "R/P"
)
entries <- rbind.fill(new_entry, entries)
entries <- entries[order(entries$entry, decreasing = FALSE), ]
entries <- replace(entries, is.na(entries), "")


############################################################################
# Plot forest plot

forest_plot <- ggplot() +
    theme_classic() +
    layers

# rename axes
forest_plot <- forest_plot +
    labs(x = effect_size, y = 0)


# Add effect line
line_data <- data.frame(x = c(coef(overall_res), coef(overall_res)), y = c(0, n_entries - 1 + 0.5))
forest_plot <- forest_plot +
    geom_hline(yintercept = n_entries - 1 + 0.5, colour = "#696767") +
    scale_color_identity() +
    geom_line(data = line_data, aes(x = x, y = y), linetype = "dashed", colour = "#5698A3", linewidth = 1)

forest_plot <- forest_plot +
    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
    ) +
    coord_cartesian(ylim = c(1, n_entries))

############################################################################
# Plot study labels and extra info that will be on the left of the forest plot

p_labels <-
    entries |>
    ggplot(aes(y = rev(entry)))
p_labels <- entries |>
    ggplot(aes(y = seq(max(entry), min(entry), by = -1)))


# add the study and extra info as text
p_labels <-
    p_labels +
    geom_hline(yintercept = n_entries - 1 + 0.5, colour = "#696767") +
    scale_color_identity() +
    draw_labels(pos = 0, key = entries$label, hjust = 0, label = entries$label, subgroups = names(res_dict)) +
    draw_labels(pos = 1.5, key = entries$horizon, hjust = 0.5, label = entries$label) +
    draw_labels(pos = 2.5, key = entries$size, hjust = 0.5, label = entries$label) +
    draw_labels(pos = 3, key = entries$approach, hjust = 0.5, label = entries$label)

# remove the background and edit the sizing so that this left size of the plot will match up neatly with the middle and right sides of the plot
p_labels <-
    p_labels +
    theme_void() +
    coord_cartesian(xlim = c(0, 3), ylim = c(1, n_entries))


############################################################################
# Plot effect sizes and extra info that will be on the right of the forest plot

p_annot <- entries |>
    ggplot(aes(y = seq(max(entry), min(entry), by = -1)))


p_annot <- p_annot +
    draw_labels(pos = 0.5, key = entries$weight, hjust = 0.5, label = entries$label) +
    draw_labels(pos = 1.5, key = entries$metric, hjust = 0.5, label = entries$label)


# remove the background and edit the sizing so that this left size of the plot will match up neatly with the middle and right sides of the plot
p_annot <- p_annot +
    geom_hline(yintercept = n_entries - 1 + 0.5, colour = "#696767") +
    scale_color_identity() +
    theme_void() +
    coord_cartesian(xlim = c(0, 2), ylim = c(1, n_entries))

# ############################################################################

layout <- c(
    area(t = 0, l = 0, b = 30, r = 5), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
    area(t = 0, l = 6, b = 30, r = 9), # starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
    area(t = 0, l = 10, b = 30, r = 13) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)


p <- p_labels + forest_plot + p_annot + plot_layout(design = layout)
plot(p)
ggsave(
    paste("forest_", paste(subgroup_strategy, effect_size, sep = "_"), ".pdf", sep = ""),
    plot = p,
    height = n_entries * 0.5,
    # this keeps more or less the same distance between plot entries/rows
)
