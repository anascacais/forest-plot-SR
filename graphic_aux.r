library(metafor)
library(ggplot2)
library(plyr)

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
        color = "#696767",
        fill = "transparent",
        inherit.aes = FALSE
    )
}


draw_study_results <- function(dat_subgroup, layers) {
    new_layers <- list(
        geom_point(
            data = dat_subgroup, aes(x = yi, y = yn), shape = 18, size = 5, colour = "#696767"
        ),
        geom_linerange(data = dat_subgroup, aes(xmin = ci_low, xmax = ci_upp, y = yn), colour = "#696767")
    )
    return(layers <- c(layers, new_layers))
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
                round(res$I2, digits = 2), "%, ",
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
        subtotal = subtotal
    )
    entries <- rbind.fill(entries, new_entries)

    # Create layers, adding the mean metric, CI and summary polygon
    new_layers <- list(
        draw_summary_polygon(summary_polygon)
    )
    layers <- c(layers, new_layers)

    print(paste(
        fmtx(coef(res), digits = 2),
        " (\\gls{CI}=", fmtx(res$ci.lb, digits = 2), "-", fmtx(res$ci.ub, digits = 2),
        ", I$^2$=", fmtx(res$I2, digits = 2), "\\%)",
        sep = ""
    ))

    return(list(entries, layers))
}


draw_labels <- function(pos, key, hjust, labels, subgroups = c(NA)) {
    if (anyNA(subgroups)) {
        fontface_strategy <- ifelse(grepl("First author, year", labels), "bold", "plain") # corresponds to pos = 1
    } else {
        fontface_strategy <- ifelse(grepl("Subtotal", labels) | grepl("First author, year", labels) | grepl("Overall", labels) | (labels %in% subgroups), "bold", "plain")
    }

    geom_text(
        aes(x = pos, label = key),
        hjust = hjust,
        fontface = fontface_strategy
    )
}
