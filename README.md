# non-Poisson-calculations
short version for explanation purposes
libraries gamlss and bigsplines are necessary for the negative binomial and splines.

The code is comprised of several functions and the main part which iterates along the passes. The data must be in a csv format as shown in the C14_C12_example_1 & _2. 
The format is three columns a passes label, 12C current and 14C counts with the headers "label", "raw_curr" and "raw_cts", respectively.
The label can be any integer number or string character that repeats for every cycle but is different for each pass.
The following functions are located in the beginning: 
mu_fitting: takes the counts and current and it is in charge of variating the smoothness parameter (lambda), number of knots (nknots_vect) to calculate a AIC value (loss function).
aic_lambda: actually fits the cubic spline to the 12C for a selected lambda, scale the 12C spline into 14C counts and optimize the proportional scale parameter A_scale using the ideal_cts_calc function.
ideal_cts_calc: Function the is used to optimise the proportional scale parameter by variating the A_scale and calculating the AIC value using the function my_pois_aic and my_pois_likelihood.
save_aic_data: "save" the current AIC value and current parameter values in a data frame.
chi2_transformation: performs the chi squared transformation from non-stationary to stationary counts.

Main code
iterates the passes. The code figures out the passes groups by checking the repeated strings in the column "labels".
Finds the optimum J and mu models using the function mu_fitting.
After finding the J and mu models for all the passes, it calculates the quassi Poisson phi and performs the non-stationary to stationary using the chi2_transformation function

