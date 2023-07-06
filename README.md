# Joint OU factor model and cumulative risk model
Code for fitting dynamic Ornstein-Uhlenbeck (OU) factor model and cumulative risk model using a weighted two-stage approach.  In the first stage, the OU factor model is fit using just the longitudinal data.  Code for carrying out the first stage of estimation is available [here](https://github.com/madelineabbott/OUF).

This repository constains code for the second stage of estimation.  To fit the cumulative risk model given known parameter estimates for the factor model (from stage 1 of estimation), run **stage2.R**.  This R file carries out the following steps:

1. Factor scores are predicted at the measurement occasions using the factor model parameters.
2. Individuals within the dataset are bootstrapped.
3. Synthetic average values of the latent process are generated for each event interval, given predicted factor scores, for each bootstapped dataset.
4. A weighted Poisson regression model is fit to each bootstrapped dataset.
5. Point estimates and standard errors are pooled across the bootstrapped datasets.

Two example datasets are provided in in the [data](/data) directory.  These datasets differ in true OU process that underlies the observed longitudinal outcomes--setting 1 corresponds to an easier setting for estimation (lower noise, higher autocorrelation) and setting 2 corresponds to a more difficult setting for estimation (higher noise, lower autocorrelation).
