library(ggplot2)
library(metafor)
library(dplyr)
library(plyr)

get_effect_size <- function(dat, effect_size_mean, effect_size_sd) {
    #' Get the effect size and its standard error
    #'
    #' @param dat dataframe with the data
    #' @param effect_size string with the name of the effect size


    dat <- mutate(dat,
        SE = dat[[effect_size_sd]] / sqrt(dat$"# Patients")
    )
    dat <- mutate(dat, vi = dat$SE^2)
    dat <- mutate(dat, yi = dat[[effect_size_mean]])

    colnames <- colnames(dat)
    dat <- transform(dat, ci_low = dat$yi - 1.96 * sqrt(dat$vi), ci_upp = dat$yi + 1.96 * sqrt(dat$vi))
    colnames(dat) <- c(colnames, c("ci_low", "ci_upp"))

    return(dat)
}


prepare_subgroup_analysis <- function(dat, subgroup_strategy) {
    # Order data according to subgroup and yi
    dat <- mutate(dat, subgroup = dat[[subgroup_strategy]])
    # dat <- dat[order(dat$"First author, year", decreasing = TRUE), ]

    subgroups <- sort(unique(dat$subgroup), decreasing = TRUE)
    subgroup_counts <- table(dat$subgroup)

    res_dict <- setNames(lapply(subgroups, function(x) list(count = subgroup_counts[[x]])), subgroups)

    # Entry for titles (last entry) + empty entry after each subgroup summary
    n_entries <- length((dat$yi)) + length(subgroups) * 2 + 2

    return(list(dat, res_dict, n_entries, subgroups))
}

get_summary <- function(dat, y) {
    ### fit random-effects model
    res <- rma.uni(dat$yi, dat$vi, data = dat, method = "REML")

    summary_polygon <- data.frame(
        x = c(res$ci.ub, coef(res), res$ci.lb, coef(res)),
        y = c(y, y + 0.25, y, y - 0.25)
    )

    return(list(summary_polygon, res))
}



draw_summary_polygon <- function(summary_polygon) {
    geom_polygon(
        data = summary_polygon, aes(x = x, y = y),
        # fill = "gray", alpha = 0.3,
        color = "black",
        fill = "transparent",
        inherit.aes = FALSE
    )
}


draw_study_results <- function(dat_subgroup) {
    list(
        geom_point(
            data = dat_subgroup, aes(x = yi, y = yn), shape = 15, size = 3
        ),
        geom_linerange(data = dat_subgroup, aes(xmin = ci_low, xmax = ci_upp, y = yn))
    )
}



get_studies_info <- function(dat_subgroup, study_label_strategy, pos) {
    summary_info <- data.frame(
        entry = pos,
        label = dat_subgroup[[study_label_strategy]],
        weight = dat_subgroup$weights,
        metric = paste(fmtx(dat_subgroup$yi, digits = 2), " (", fmtx(dat_subgroup$ci_low, digits = 2), "-", fmtx(dat_subgroup$ci_upp, digits = 2), ")", sep = ""),
        horizon = dat_subgroup$"Forecast horizon",
        approach = ifelse(dat_subgroup$"Training and testing approach" == "prospective", "(*)", "")
    )

    return(summary_info)
}

get_summary_info <- function(res, pos, weight, subtotal = TRUE) {
    label <- ifelse(subtotal, "Subtotal", "Overall")

    summary_info <- data.frame(
        entry = pos,
        label =
            paste(
                label,
                " (I^2 =",
                round(res$I2, digits = 3), "%, ",
                fmtp(res$QMp, digits = 2, pname = "p", add0 = TRUE, sep = TRUE, equal = TRUE),
                ")"
            ),
        weight = weight,
        metric = paste(fmtx(coef(res), digits = 2), " (", fmtx(res$ci.lb, digits = 2), "-", fmtx(res$ci.ub, digits = 2), ")", sep = "")
    )

    return(summary_info)
}


make_summary <- function(dat, entries, layers, n_entries, pos, subtotal = TRUE) {
    output <- get_summary(dat, pos)
    summary_polygon <- output[[1]]
    res <- output[[2]]

    new_entries <- get_summary_info(
        res,
        pos = ifelse(subtotal, n_entries - pos + 1, n_entries),
        weight = ifelse(subtotal, sum(as.double(dat$weights)), "100.00"),
    )
    entries <- rbind.fill(entries, new_entries)

    # Create layers for this subgroup, adding the mean metric, CI and summary polygon
    if (subtotal) {
        new_layers <- list(
            draw_summary_polygon(summary_polygon),
            draw_study_results(dat)
        )
    } else {
        new_layers <- list(
            draw_summary_polygon(summary_polygon)
        )
    }
    layers <- c(layers, new_layers)

    return(list(entries, layers))
}

# output <- get_summary(dat, 1)
# summary_polygon <- output[[1]]
# overall_res <- output[[2]]

# new_entries <- get_summary_info(
#     overall_res,
#     pos = n_entries,
#     weight = "100.00",
#     subtotal = FALSE
# )
# entries <- rbind.fill(new_entries, entries)

# summary_layer <- list(
#     draw_summary_polygon(summary_polygon)
# )
# layers <- c(layers, summary_layer)
