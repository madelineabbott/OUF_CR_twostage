# Joint OU factor model and cumulative risk model
Code for fitting dynamic Ornstein-Uhlenbeck (OU) factor model and cumulative risk model using a weighted two-stage approach.  Two sets of example data are provided in the [data](/data) directory.  Each set consists of a file containing the longitudinal outcome data and a file containing the cumulative event outcome data.  These sets of example data differ in the true OU process that underlies the observed longitudinal outcome; setting 1 corresponds to a setting in which the true OU process has higher correlation and setting 2 corresponds to a setting in which the true OU process has lower correlation.

In **stage 1**, the OU factor model is fit using just data on the longitudinal outcomes.  To carry out stage 1, run **stage1.R**.  More details on stage 1 estimation can be found [here](https://github.com/madelineabbott/OUF).

In **stage2.R**, the cumulative risk model.  After carrying out stage 1 of estimation, stage 2 can be completed by running **stage2.R**.  This R file uses results from stage 1 to carry out the following steps:

1. Factor scores are predicted at the measurement occasions using the dynamic factor model parameters.
2. Individuals within the dataset are bootstrapped.
3. Synthetic average values of the latent process are generated for each event interval, conditional on the factor scores predicted at measurement occasions, for each bootstrapped dataset.
4. A weighted Poisson regression model is fit to each bootstrapped dataset augmented with the synthetic average values of the latent process.
5. Point estimates and standard errors are pooled across the bootstrapped datasets.

For questions, please contact mrabbott@umich.edu.
