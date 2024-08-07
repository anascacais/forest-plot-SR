# Forest plots with ggplot2

## Why forest plots?

Forest plots provide a graphical summary of multiple individual results, allowing for quick visual assessment. They are a flexible tool in evidence synthesis, allowing for interpretation of complex data in a coherent and comprehensive manner. Uses of forest plots include:

- Comparison of results and confidence intervals, making it easier to identify consistency or variability in results
- Estimation of overall result
- Stratification by specific factors

## What we are trying to understand?

In this case, I used forest plots in the context of a meta-analysis to summarize the current performance of automated algorithms for the forecast of seizure risk. The main objectives were to answer the following questions:

1. What is the benchmark performance of automated algorithms for forecast of seizure risk?
2. Which data are the most valuable biomarkers for seizures?
3. Which algorithm design factors provide more informative forecasts?

## Exploring data

As a newcomer to R, I started on a post by Katherine Hoffman (https://www.khstats.com/blog/forest-plots/), which used 𝑚𝑒𝑡𝑎𝑓𝑜𝑟, 𝑔𝑔𝑝𝑙𝑜𝑡2, and 𝑝𝑎𝑡𝑐ℎ𝑤𝑜𝑟𝑘. But, I was finding it challenging to customize the visualization (namely, when it came to **subgroup analysis** and adding **algorithm characteristics** to it).

So I adapted the original code to, without the need to modify the original data (spreadsheet), do the following:

- Perform subgroup analysis by providing only the name of the column with the wanted factor
- Expand visualization with algorithm characteristics (forecast horizon, sample size, number of seizures, and train/test approach)
- Customize colors and symbols

To

## Interpreting the results

1.

[Example](results/forest_Data%20source_BSS.png)
