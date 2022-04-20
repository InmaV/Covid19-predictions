# Covid19-predictions

Data_EU+EFTA+UK_2020.xlsx and Data_EU+EFTA+UK_2021.xlsx were downloaded from the World Health Organization webpage () and contain the cummulative COVID19 cases reported in Europe in 2020 and 2021, respectively. A minor preprocessing was done for Data_EU+EFTA+UK_2020.xlsx, in which countries that presented only one gap in their data were corrected by distributing the extra cases reported the following day to the previous. 

Patterns.m determines the daily reporting pattern for each country in a specified period of time. Predictions.m provides the COVID19 cases predictions computed by adjusting a Gompertz to their data, that can be outputted as csv files. Four different models can be used: model B (Baseline), using the raw data without corrections and the cumulative minimization; model I (Introduction of Patterns), using the pattern-weighted data and the cumulative minimization; model F (Fallback), using the raw data without corrections and both the minimization of the cumulative and new cases; and finally, model H (our Hallmark), using the pattern-weighted data and both the minimization of the cumulative and new cases. Predictions with less than a weekly average of 100 daily cases and/or close to a holiday are filtered out. The number of days introduced to the fitting function (Nfitting) can be varied to obtain their predictions and errors, collected in Analysis_cases_allN.xlsx and Errors_N.xlsx, respectively. These files can be employed in Instability_analysis.m to determine the unstable predictions (Analysis_stability.xlsx) and filter them out from Predictions.m.
