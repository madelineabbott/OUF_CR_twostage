# Joint OU factor model and cumulative risk model
Code for fitting dynamic Ornstein-Uhlenbeck (OU) factor model and cumulative risk model using a weighted two-stage approach.

+ To fit the cumulative risk model given known parameter estimates from the dynamic factor model, run stage2.R.  This file bootstraps the data, generate sythetic average values of the latent process, and then combines the point estimates and returns pooled standard error estimates based on the von Hipple approach.
