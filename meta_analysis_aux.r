library(dplyr)


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


prepare_subgroup_analysis <- function(dat, subgroup_strategy, subgroup_analysis = TRUE) {
    # Order data according to subgroup and yi
    dat <- mutate(dat, subgroup = dat[[subgroup_strategy]])

    subgroups <- sort(unique(dat$subgroup), decreasing = TRUE)
    subgroup_counts <- table(dat$subgroup)

    res_dict <- setNames(lapply(subgroups, function(x) list(count = subgroup_counts[[x]])), subgroups)

    # if there are subgroups with less than 2 entries, subgroup analysis in not performed
    if (any(sapply(res_dict, function(x) x$count <= 2))) {
        subgroup_analysis <- FALSE
    }

    # Entry for titles (last entry) + empty entry after each subgroup summary
    n_entries <- length((dat$yi)) + length(names(res_dict)) * (1 + 1 * subgroup_analysis) + 2

    return(list(dat, res_dict, n_entries))
}
