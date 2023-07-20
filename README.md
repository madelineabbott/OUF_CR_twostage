# Joint OU factor model and cumulative risk model
Code for fitting dynamic Ornstein-Uhlenbeck (OU) factor model and cumulative risk model using a weighted two-stage approach.  In the first stage, the OU factor model is fit using just data on the longitudinal outcomes.  Code for carrying out the first stage of estimation is available [here](https://github.com/madelineabbott/OUF).

This repository contains code for the second stage of estimation.  To fit the cumulative risk model given known values for parameters in the dynamic factor model (from stage 1 of estimation), run **stage2.R**.  This R file carries out the following steps:

1. Simulated data--for both the longitudinal and event outcomes--are imported from the [data](/data) directory. 
2. Factor scores are predicted at the measurement occasions using the dynamic factor model parameters.
3. Individuals within the dataset are bootstrapped.
4. Synthetic average values of the latent process are generated for each event interval, conditional on the factor scores predicted at measurement occasions, for each bootstrapped dataset.
5. A weighted Poisson regression model is fit to each bootstrapped dataset augmented with the synthetic average values of the latent process.
6. Point estimates and standard errors are pooled across the bootstrapped datasets.

Two sets of example data are provided in the [data](/data) directory.  Each set consists of a file containing the longitudinal outcome data and a file containing the cumulative event outcome data.  These sets of example data differ in the true OU process that underlies the observed longitudinal outcome; setting 1 corresponds to an easier setting for estimation (lower noise, higher autocorrelation) and setting 2 corresponds to a more difficult setting for estimation (higher noise, lower autocorrelation).

For questions, please contact mrabbott@umich.edu.
