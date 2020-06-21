# Histogram-Uncertainty

This work has been published by **Mehdi Rezvandehy** and **Clayton V. Deutsch** at Springer,
[Journal of Natural Resource Research](https://doi.org/10.1007/s11053-016-9322-3) https://doi.org/10.1007/s11053-016-9322-3.
**Python** implementation is presented here for calculating correct histogram uncertainty. If you have any question about the approach and implementation, please email rezvande@ualberta.ca.

**There is unavoidable uncertainty in the representative histogram required as an input for geostatistical modeling. This uncertainty should be quantified correctly and incorporated into final modeling because it affects resource/reserve estimation, investment and development decisions. This works confirms a correct approach of quantifying histogram uncertainty by comparing to reference uncertainty found by an automated scan-based approach. Similar patterns of a data configuration are found within a large image and the mean of the specified domain is directly calculated. The correct quantification of histogram uncertainty is defined. The spatial bootstrap provides prior uncertainty that does not consider the domain limits and is not conditioned to the data.  This uncertainty is updated in geostatistical modeling to consider the limits and data. The resulting uncertainty closely matches the scan based approach.**

## Introduction

Developing natural resourcea (mining and hydrocarbone) are risky. The industry aims to predict and mitigate risk. Exploration and production companies shows a decreasing rate of success for exploration and production megaprojects. One of the main reasons for this underperformance is due to use of evaluation methods that do not account for the full uncertainty, which leads to inaccurate production forecasts. In the past, the input statistics were held constant and relatively small fluctuations between realizations were used to characterize reservoir uncertainty. This approach underestimates uncertainty. Uncertainty is small because local fluctuations above and below the average cancel out between locationS. Thus, the uncertainty in input parameters is important and should be integrated into final resource modeling. There are a number of approaches that have been proposed to quantify uncertainty in the univariate distribution (histogram).

The bootstrap is the first simplest method of quantifying uncertainty in the histogram developed. This method uses MCS simulation to draw values from the data distribution to simulate different possible data sets; so, it can be easily applied to calculate the uncertainty in the mean and other statistical parameters. There are two critical assumptions for applying the bootstrap: 1- The distribution of the data should be representative of the whole domain, and 2- The data are independent. The bootstrap may be useful early in appraisal with widely spaced well data.

The spatial bootstrap was proposed in order to consider the spatial correlation of data. The spatial bootstrap in geostatistics applies unconditional simulation at the data locations according to spatial correlation of the data. This approach considers neither the conditioning data nor the area of interest. Increasing spatial correlation  leads to greater uncertainty because the data are more redundant. 

The last technique for assessing uncertainty in the mean is using kriging for estimation of the entire domain. The global kriging variance will decrease when the domain size increases due to the support effect. This technique is independent of data values and leads to relatively low uncertainty.

Khan and Deutsch, 2015 propose a simulation-based approach for quantifying histogram uncertainty. This method transfers the spatial bootstrap uncertainty to posterior uncertainty considering the conditioning data and domain size. They claim that this posterior uncertainty is the correct histogram uncertainty which is between the global kriging uncertainty as lower limit and the spatial bootstrap as an upper limit; this proposal has not been proved in practice. This paper work an experimental framework to evaluate the approaches of quantifying histogram uncertainty.

## Workflow for Histogram Uncertainty

The histogram uncertainty quantified by the spatial bootstrap is assumed to be the prior uncertainty. Considering this prior uncertainty in simulation conditioned to the original data and limited to the domain limits achieves more reliable histogram uncertainty. The result is called posterior histogram uncertainty. The procedure is summarized by:

1. Define a stationary covariance function.
2. Define the reference distribution.	
3. Perform the spatial bootstrap resampling as follows:		
	3-a. Construct the spatial data-to-data covariance matrix. 
	3-b. Compute the Cholesky decomposition of the correlation matrix as. 		
	3-c. Simulate a vector of uncorrelated standard normal deviate. 		
	3-d. Generate a vector of correlated Gaussian values. 	 	
	3-e. Transform the unconditional Gaussian values to original units.		

4. Conditional simulation; the spatial bootstrap realization used as a reference distribution for normal score transformation.	
5. Back transform the realization to original units with the spatial bootstrap reference distribution.	
6. Each realization is limited to the domain limits.	

Steps 3.c to 6 are repeated to achieve the posterior histogram uncertainty. This approach accounts for the conditioning data, the domain limits and the spatial correlation between data.

Figure 1 shows a 2D synthetic data set with in an area of $1000m\times1000m $ with nine data locations. The variable is effective porosity. The uncertainty in the mean of the distribution is calculated by the prior and posterior uncertainty. Sequential Gaussian simulation (SGS) is used for conditional simulation. An anisotropic variogram model is assumed with ranges of 800m for azimuth $0^{\circ}$ and 400m for azimuth $90^{\circ}$. Prior uncertainty in the mean is attained by averaging each realization of the spatial bootstrap. The variance of the mean is $2.11\times 10^{-2}$. Posterior uncertainty in the mean is attained by using each realization of the spatial bootstrap as a reference distribution for normal score transformation of data and back transforming the realizations from the Gaussian space to original units. The mean of each realization is computed. This process is repeated to attain the posterior uncertainty in the mean: the variance of the mean decreases from $2.11\times 10^{-2}$ (prior uncertainty) to $1.43\times 10^{-2}$ because of the conditioning data and domain limits. The posterior histogram uncertainty is claimed to be more accurate than the other techniques. The only way to check the approaches of quantifying histogram uncertainty is to design an experimental framework where the true uncertainty in the histogram is known.

<img src="./Images/fig1.png" alt="drawing" width="650"/>

**Figure 1**: *A synthetic 2D example of effective porosity. The uncertainty in the mean is calculated by prior and posterior uncertainty. The posterior approach leads to lower uncertainty in the mean of effective porosity: the variance of the mean decreases from $2.11\times 10^{-2}$ (prior uncertainty) to $1.43\times 10^{-2}$ because of the conditioning data and domain limits.*
