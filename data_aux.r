library(dplyr)

prepare_dataframe <- function(data, subset_strategy, subset_name, effect_size_mean, effect_size_sd) {
    #' Get subset of data according to subset_strategy and studies that reported the SD of the variable of interest;
    #' convert columns to the correct type; and create a key for each study for the meta-analysis
    #'
    #' @param data dataframe with the data
    #' @param subset_strategy string with the name of the column to subset the data
    #' @param subset_name string with the value to subset the data
    #' @param effect_size_mean string with the name of the column with the mean of the effect size
    #' @param effect_size_sd string with the name of the column with the standard deviation of the effect size
    #' @param key_strategy string with the name of the column to create the key for the meta-analysis

    data <- subset(
        data, (data[subset_strategy] == subset_name)
    )

    data <- subset(
        data, (!is.na(data[[effect_size_sd]])) & (data[[effect_size_sd]] != 0)
    )

    # Get only required columns
    dat <- data[c(
        "First author, year",
        "# Patients",
        "Type of input data",
        "Input data",
        "Forecast horizon",
        "Training and testing approach",
        effect_size_sd,
        effect_size_mean
    )]


    dat <- mutate(dat, key = paste(
        dat$"First author, year", dat$"Forecast horizon", dat$"Input data", dat$"Training and testing approach",
        sep = ", "
    ))

    dat$key <- as.character(dat$key)
    dat[[effect_size_sd]] <- as.numeric(dat[[effect_size_sd]])
    dat[[effect_size_mean]] <- as.numeric(dat[[effect_size_mean]])
    dat$"# Patients" <- as.numeric(dat$"# Patients")

    return(dat)
}
