#*****
# Q1 *
#*****
# Forecast cumulative cases and cumulative deaths for MechBayes (SEIRD?) model, no intervention. Provide 90% prediction intervals
# Compare accuracy to true data over span

#*****
# Q2 *
#*****
# What is probability hospitalizations stay under deaths threshold for interval

#*****
# Q3 *
#*****
# If institute mask req, what is prob of staying below threshold. 
#  May need to incorporate masking
#  May want to minimize decrease to meet threshold
#  Use TA1 to translate intervention into transmission parameter, including uncertainty

#*****
# Q4 *
#*****
# Detection rate param, varies as Gauss RW on log-odds scale, likely captures other changes as well
# a) If increase detection by 20% (no other interventions), does that increase forecasted cumulative cases and deaths?
#    How does increase in detection affect the uncertainty of estimates
#    Characterize relationship b/w detection rate, forecasts, uncertainties
#    Does improving detection rates provide decision makers with more accurate forecasts or narrower prediction intervals
# b) Compute accuracy of forecasts (no mask policy) and determine if detection rate improves accuracy

#*****
# Q5 *
#*****
# Model modification, model space exploration, and model run_model_selection
# a) Convert SEIRHD to SEIRHD
#    Compute forecasts and compare accuracy to those of 1
# b) Further modify and do model space exploration and run_model_selection
#    Fit to previous month. Select from forecast accuracy
# c) Three-way structural comparison

#*****
# Q6 *
#*****
# Modify SEIRD to SEIRD+Renewal
# Explain and plot how models differ

#*****
# Q7 *
#*****
# Latest mask mandate can be imposed to ensure with 90% prob won't exceed thresh 
# Characterize relationship between extra cumulative deaths for each day delay