# Disruption of spatiotemporal dependence in dengue transmission by wMel Wolbachia in Yogyakarta, Indonesia

This repository holds all of the code necessary to recreate the analysis performed in the paper by the same name. The code has not been optimized, but should accurately return the results described in the paper.


## Organization of this Repository
This repository is organized as follows. Note: the raw data is not available in this repository. However, the data necessary to recreate most of the figures (excluding those relying on geolocated residence) is available in the `graphs` folder.

### `munge`

Contains all of the scripts that run the analyses on the raw data and return the results.

+ `overall-bootstrap-script.R` calls the following functions to run tau for the study area, naive to intervention designation. Significance is determined as the result of bootstrap resampling.
   + `lib/bootstrap-function_odds.R`
   + `lib/binary-matrix-function_odds.R`
   + Output: `data/overall-bootstrap-output.RData` 
   
+ `overall-permutation-script.R` calls the following functions to run tau for the study area, naive to intervention designation. Significance is determined as compared to the 95% confidence intervals on the permutation-based null distribution.
   + `lib/permutation-function.R`
   + `lib/binary-matrix-function_odds.R`
   +  Output: `data/overall-permutation-output.RData`
   
+ `cluster-specific-permutation-script.R` calls the following functions to run tau at the cluster level. Significance is determined as compared to the 95% confidence intervals on the permutation-based null distribution.
   + `lib/permutation-function.R`
   + `lib/binary-matrix-function_odds.R`
   + Output: `data/cluster-specific-permutation-output.RData`
   
+ `cluster-specific-spatial-permutation-sensitivity-script.R` and `cluster-specific-temporal-sensitivity-script.R` call on the following functions to run the sensitivity checks. Significance is determined as compared to the 95% confidence intervals on the permutation-based null distribution.
   + `lib/permutation-function.R`
   + `lib/binary-matrix-function_odds.R`
   + Output: `data/cluster-specific-permutation-output-temporal-sensitivity.RData`
   + Output: `data/cluster-specific-permutation-output-sensitivity-d50.RData`
   
### `analysis`
Turns the output from `munge` into the tables and figures presented in the article.

### `graphs`

Figures were generated with the following files:

+ **Figure 1.** Time series plot for illness onset among A) test-negative controls and B) virologically-confirmed dengue cases included in the primary analysis of the AWED trial in Yogyakarta, Indonesia from January 2018 until March 2020, by intervention arm. No dengue cases were enrolled in September 2018 and, in accordance with the trial protocol \citep{Anders2020update}, the test-negatives enrolled during that month were excluded from the analysis dataset. 
   + `analysis/04_time-series-plots.R`

+ **Figure 2.** Spatial distribution of A) enrolled dengue cases by serotype across Yogyakarta City, B) the cluster-aggregate test-positive fraction, i.e., the proportion of enrolled dengue cases among the total number of individuals enrolled in each cluster, and C) kernel smoothing estimates of the spatially-varying test-positive fraction. Each map includes participants enrolled from January 2018 through March 2020. The borders in each map represent the cluster boundaries for the AWED trial. Clusters are numbered with their administrative labels. Points represent the geolocated households of virologically confirmed dengue cases. Areas with darker shading are associated with a higher proportion of dengue cases among the AWED participants than areas with lighter shading. Smoothing bandwidth was selected by cross-validation.   
   + `analysis/00_overall-map.R`
   
+ **Figure 3.** Estimated odds ratio ($\tau (d_1,d_2)$) comparing the odds of a homotypic dengue case pair within $(d_1,d_2)$ versus the odds of a homotypic dengue case pair at any distance across the entire study area among participant pairs with illness onset occurring within 30 days with A) bootstrap 95\% confidence interval and B) against the 95\% CI on the permutation-based null rejection region.   
   + `analysis/02_overall-analysis.R`
   
+ **Figure 4.** Cluster-specific and pooled arm-level estimates of the $\tau$-statistic (points) and 95\% CIs on the null distribution (error bar) generated from 1,000 simulations, where the location at which a case occurs is randomly reassigned within each cluster. Each panel displays the estimated spatial dependence for homotypic case pairs with illness onset occurring within 30 days and resident within a given distance interval (meters) from each other . Statistically significant dependence is present when the point estimate falls outside of the 95\% CIs of the null distribution and, for improved visibility, is marked by the light blue points. The overall point estimate for each trial arm is found by taking the geometric mean of the cluster-level estimates and is then compared against the 95\% CIs of the null distribution of the permuted geometric mean.   
   + `analysis/07_cluster-permutation-analyses.R`
   
+ **Figure 5.** Residential locations of the enrolled serotyped dengue cases involved in homotypic pairs with residences within 300m and illness onset within 30 days, including pairs that cross cluster boundaries.   
   + `analysis/10_homotypic-vcd-maps.R`
   
+ **Figure 6.** Sensitivity analyses and comparison with the primary analysis (‘Full data’). Estimated geometric mean odds ratio, $\tau$, comparing the odds of a homotypic dengue case pair within the distance interval $(d_1, d_2)$ versus the odds of a homotypic dengue case pair at any distance across the entire study area for 1) the full dataset, 2) the dataset excluding those within 50m of a cluster border, 3) participant pairs with illness
onset occurring within 1 week of each other, and 4) participant pairs with illness onset occurring within 2 weeks of each other. The shaded area is the 95% CI of the permutation-based null distribution.
   + `analysis/09_sensitivity-comparisons.R`