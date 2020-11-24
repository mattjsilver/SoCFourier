# SoCFourier
Gambian DNA methylation - season of conception analysis using Fourier regression

Code for running Fourier regression models described in 
Silver et al. Environmentally sensitive hotspots in the methylome of the
early human embryo: https://www.biorxiv.org/content/10.1101/777508v2

Methylation and covariate data from the ENID (2yr) cohort is available at:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99863  
Data from the replication EMPHASIS (8-9yr) cohort will be made available once the results from the main EMPHASIS 
study are published.
<br/><br/>

The relationship between DNA methylation (outcome) and date of conception
(predictor) is modelled using Fourier regression with two Fourier terms:
sin(DoC) and cos(DoC) where DoC is date of conception measured in radians
(0=1st Jan; 2*pi=31st Dec). Adjustment covariates are specified by the user
in design matrix required as input.

Fourier model with 2 Fourier terms assumes a sinusoidal relationship between
methylation and DoC with a single methylation maximum and minimum within a
single year. The amplitude (distance between maximum and minimum)
and phase (time of year of methylation maximum/minimum) are determined by
the estimated coefficients for each Fourier term.

The significance of any 'seasonal' (date of conception) association is determined by likelihood ratio test (LRT) comparing the full model with both
Fourier terms against a nested, covariates-only model.
