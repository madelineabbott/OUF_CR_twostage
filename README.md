# Joint OU factor model and cumulative risk model
Code for fitting dynamic Ornstein-Uhlenbeck (OU) factor model and cumulative risk model using a weighted two-stage approach.  In the first stage, the OU factor model is fit using just the longitudinal data.  Code for carrying out the first stage of estimation is available \href{here}{https://github.com/madelineabbott/OUF}.

The repository constains code for the second stage of estimation.  To fit the cumulative risk model given known parameter estimates for the factor model (from stage 1 of estimation), run **stage2.R**.  This R file carries out the following steps:

1. The data are bootstrapped
2. Synthetic average values of the latent process are generated for each event interval
3. A weighted Poisson regression model is fit to each bootstrapped dataset
4. Point estimates and standard errors are pooled across the bootstrapped dataset

Two example datasets are provided in in the data directory.  These datasets differ in the values of the parameters for the true OU process that underlies the observed longitudinal outcomes.
