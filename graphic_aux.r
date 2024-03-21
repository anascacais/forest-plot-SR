library(ggplot2)

draw_labels <- function(pos, key, hjust, labels, subgroups = c(NA)) {
    if (anyNA(subgroups)) {
        fontface_strategy <- ifelse(grepl("Study", labels), "bold", "plain") # corresponds to pos = 1
    } else {
        fontface_strategy <- ifelse(grepl("Subtotal", labels) | grepl("Study", labels) | grepl("Overall", labels) | (entries$label %in% subgroups), "bold", "plain")
    }

    geom_text(
        aes(x = pos, label = key),
        hjust = hjust,
        fontface = fontface_strategy
    )
}
