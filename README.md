# partial_residual_plots_lavaan
This project contains one function that creates partial residual plots for single paths of SEMs created with lavaan. It may or may not work for you.

# The function creates partial residual plots for single paths of SEMs created with lavaan
# It can handle multigroup and multilevel models and deal with funky survey designs fit with lavaan.survey
# because I am fancy, it can also convert data back to its original scale (if you did your sem right, you scaled your data beforehand ;-))
# and it can backtransform log and sqrt transformed y and x variables (it does so correctly if you log transformed them before the scaling, not after)

# requried input:
# - fit   fitted lavaan object (this needs to be fitted with lavaan, not lavaan.survey!)
# - y     variable to which the path goes from x
# - x                 variable from which the path goes to y
#
# voluntary input:
# - data              sometimes the function fails to grab the (scaled) data from the model output. You can feed it to the function manually
# - de_scale          set to TRUE if you want the plot to display data on its original units (before the scaling)
# - unscaled_data     you need to provide the unscaled data if you want the back-transformation to work
# - survey            set to TRUE if you want to account for the survey design with a lavaan.survey generated model fit. 
# - survey_fit        you need to provide the object fitted with lavaan.survey (in addition to the object fitted with lavaan in fit) if survey is set to TRUE
# - level             if you fit a multilevel SEM, specify at which level your path of interest takes place. It does work accross levels, but creates funky wave fitted lines. Not sure thats what you want, but hey, you do you
# - backtransform_x   if you log or sqrt transformed your variable (before the scaling!), you can backtransform it. At this point you can give "log", "log+1)", "sqrt", "sqrt+1" backtransformations for log(x), log(x+1), sqrt(x) and sqrt(x+1) transformed variables respectively
# - backtransform_y   same as backtransform_x, but for y

# created in R R-4.3.2 in RStudio 2023.12.0 Build 369
# packages required: ggplot2 (Version 3.5.1), cowplot (Version 1.1.2) and RColorBrewer (Version 1.1-3)
