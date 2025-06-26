# created 24.06.2025 by seraina.cappelli@gmail.com

# This function creates partial residual plots for single paths of SEMs created with lavaan
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



the_big_beautiful_residual_plot_function <- function(fit, y, x, 
                                           data = NULL,
                                           de_scale = FALSE, 
                                           unscaled_data = NULL, 
                                           survey = FALSE, 
                                           survey_fit = NULL, 
                                           level = NULL,
                                           backtransform_x = NULL,
                                           backtransform_y = NULL) {
  
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
  theme_set(theme_cowplot() +
              theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
              theme(axis.text = element_text(size = 8),
                    axis.title = element_text(size = 10),
                    legend.text = element_text(size = 8),
                    legend.title = element_text(size = 10),
                    title = element_text(size = 10),
                    text = element_text(size = 8)))
  
  # Get variable names as character strings
  x <- deparse(substitute(x))
  y <- deparse(substitute(y))
  
  
  # Extract data from model or use provided data
  model_data <- lavInspect(fit, "data")
  
  if (is.null(model_data)) {
    if (is.null(data)) {
      stop("No data found in model object. Please provide 'data' manually using the `data = ...` argument.")
    }
    df0 <- data
    df0$group <- NA
  } else {
    if (lavInspect(fit, "ngroups") < 1.5) {
      df0 <- data.frame(model_data)
      df0$group <- NA
    } else {
      dfs <- model_data
      df0 <- cbind.data.frame(do.call(rbind.data.frame, dfs),
                              rep(names(dfs), vapply(dfs, nrow, numeric(1))))
      names(df0) <- c(names(df0)[-length(names(df0))], "group")
    }
  }
  
  # Get parameter table from correct fit
  cofis <- if (survey) {
    if (is.null(survey_fit)) stop("Please provide survey_fit when survey = TRUE")
    parTable(survey_fit)
  } else {
    parTable(fit)
  }
  
  # check if it is a multilevel model
  if ("level" %in% names(cofis)){
    if (is.null(level)) {
      warning("Please provide level for which you want the partial plots")
    } else {
      cofis <- cofis[cofis$level == level,]
    }
  }
  
  # Variable checks
  if (!(y %in% names(df0))) {
    stop(paste("Variable", y, "is not present in the data"))
  }
  if (!(x %in% names(df0))) {
    stop(paste("Variable", x, "is not present in the data"))
  }
  if (!(y %in% cofis[cofis$op == "~" & cofis$rhs == x, ]$lhs)) {
    stop(paste("Path", y, "~", x, "is not present in the model"))
  }
  
  # Get other predictors of y (excluding x)
  vars <- unique(cofis[cofis$lhs == y & cofis$op == "~", ]$rhs)
  vars <- setdiff(vars, x)
  
  # Loop through groups
  df <- df0[0, ]
  df$y_resid <- c()
  fit_line <- data.frame(x = c(), y = c(), group = c())
  
  for (g in 1:lavInspect(fit, "ngroups")) {
    group_label <- lavInspect(fit, "group.label")[g]
    df_group <- df0[df0$group %in% group_label | is.na(df0$group), ]
    cofis_group <- cofis[cofis$group == g | is.na(cofis$group), ]
    
    # Drop missing values
    df_group <- df_group[complete.cases(df_group[, c(y, x, vars)]), ]
    
    # Compute partial residuals
    df_group$y_resid <- df_group[[y]]
    intercept <- cofis_group[cofis_group$lhs == y & cofis_group$op == "~1", "est"]
    df_group$y_resid <- df_group$y_resid - intercept
    for (i in vars) {
      effect <- cofis_group[cofis_group$lhs == y & cofis_group$rhs == i, "est"]
      df_group$y_resid <- df_group$y_resid - effect * df_group[[i]]
    }
    
    # Fitted line
    minx <- min(df_group[[x]], na.rm = TRUE)
    maxx <- max(df_group[[x]], na.rm = TRUE)
    fit_line_group <- data.frame(x = seq(minx, maxx, length.out = 50))
    slope <- cofis_group[cofis_group$lhs == y & cofis_group$rhs == x, "est"]
    fit_line_group$y <- intercept + slope * fit_line_group$x
    fit_line_group$group <- group_label
    
    # Collect results
    df <- rbind.data.frame(df, df_group)
    fit_line <- rbind.data.frame(fit_line, fit_line_group)
  }
  
  # Re-add x variable for plotting
  df$x <- df[[x]]
  
  # Unscale if requested
  if (de_scale) {
    if (is.null(unscaled_data)) stop("Please provide unscaled_data when de_scale = TRUE")
    mean_x <- mean(unscaled_data[[x]])
    sd_x <- sd(unscaled_data[[x]])
    mean_y <- mean(unscaled_data[[y]])
    sd_y <- sd(unscaled_data[[y]])
    
    df[[x]] <- sd_x * df[[x]] + mean_x
    df[["x"]] <- sd_x * df[["x"]] + mean_x
    df[[y]] <- sd_y * df[[y]] + mean_y
    fit_line[[x]] <- sd_x * fit_line[[x]] + mean_x
    fit_line[[y]] <- sd_y * fit_line[[y]] + mean_y
  }
  
  # Backtransform if requested
  if (!is.null(backtransform_x)){
    if(!(backtransform_x %in% c("log", "log+1", "sqrt", "sqrt+1"))){
      warning("at this point only 'log(x)', 'log(x+1)', 'sqrt(x)' and 'sqrt(x+1)' backtransformations are available\the displayed output is not backtransformed")
    } else {
      if(backtransform_x == "log") {
        df[[x]] <- exp(df[[x]])
        df[["x"]] <- exp(df[["x"]])
        fit_line[["x"]] <- exp(fit_line[["x"]])
      }
      if(backtransform_x == "log+1") {
        df[[x]] <- exp(df[[x]]) - 1
        df[["x"]] <- exp(df[["x"]]) - 1
        fit_line[["x"]] <- exp(fit_line[["x"]]) - 1
      }
      if(backtransform_x == "sqrt") {
        df[[x]] <- (df[[x]])^2
        df[["x"]] <- (df[["x"]])^2
        fit_line[["x"]] <- (fit_line[["x"]])^2
      }
      if(backtransform_x == "sqrt+1") {
        df[[x]] <- (df[[x]])^2 - 1 
        df[["x"]] <- (df[["x"]])^2 - 1 
        fit_line[["x"]] <- (fit_line[["x"]])^2 -1
      }
    }
  }
  
  
  if (!is.null(backtransform_y)){
    if(!(backtransform_y %in% c("log", "log+1", "sqrt", "sqrt+1"))){
      warning("at this point only 'log(y)', 'log(y+1)', 'sqrt(y)' and 'sqrt(y+1)' backtransformations are available\the displayed output is not backtransformed")
    } else {
      if(backtransform_y == "log") {
        df[[y]] <- exp(df[[y]])
        df[["y_resid"]] <- exp(df[["y_resid"]])
        fit_line[["y"]] <- exp(fit_line[["y"]])
      }
      if(backtransform_y == "log+1") {
        df[[y]] <- exp(df[[y]]) - 1
        df[["y_resid"]] <- exp( df[["y_resid"]]) -1
        fit_line[["y"]] <- exp(fit_line[["y"]]) - 1
      }
      if(backtransform_y == "sqrt") {
        df[[y]] <- (df[[y]])^2
        df[["y_resid"]] <- ( df[["y_resid"]])^2
        fit_line[["y"]] <- (fit_line[["y"]])^2
      }
      if(backtransform_y == "sqrt+1") {
        df[[y]] <- (df[[y]])^2 - 1  
        df[["y_resid"]] <- ( df[["y_resid"]])^2 -1
        fit_line[["y"]] <- (fit_line[["y"]])^2 -1
      }
    }
  }
  
  # Plot
  p <- ggplot(df, aes(y = y_resid, x = x, color = group)) +
    geom_point(alpha = 0.2) +
    geom_line(data = fit_line, aes(y = y, x = x), size = 1) +
    scale_color_brewer(palette = "Set2") +
    labs(y = paste("Partial residuals of", y),
         x = x,
         title = deparse(substitute(fit)))
  
  if (lavInspect(fit, "ngroups") < 1.5) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}
