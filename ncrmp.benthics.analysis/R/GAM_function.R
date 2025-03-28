GAM_function <- function(data, response, group_var = NULL, group_names = NULL, kvalue = 3){

  library(mgcv)
  library(rlang)
  
  if (!response %in% colnames(data)) {
    stop("Response variable not found in the dataset.")
  }

  #!!sym() allows us to generalize this function
  #it calls it as a column name
  
  # If group_var is provided, filter the data
  if (!is.null(group_var)) {
    if (!group_var %in% colnames(data)) {
      stop("Grouping variable not found in the dataset.")
    }
    
    data <- data %>%
      filter(!!sym(group_var) %in% group_names) %>%
      mutate(
        !!sym(response) := as.numeric(!!sym(response)),
        SE = as.numeric(SE)
      )
  }
  
  #Make sure that SE / standard error column exists
  if (!"SE" %in% colnames(data)) {
    stop("Column 'SE' not found in the dataset. Ensure it is included.")
  }

  #build and fit the model
  formula <- as.formula(paste(response, "~ s(YEAR, k = ", kvalue, ")"))
  model <- gam(formula, weights = 1/(SE^2), data = data)

  #outputs of model 
  
  #Adjusted r squared
  adj_r <- summary(model)$r.sq

  #p-value for the smooth term
  p_val <- summary(model)$s.table[1, "p-value"]

  #R-squared and p-value
  result <- list(adj_r_squared = adj_r, p_value = p_val)

  # Statistical significance check
  significance <- ifelse(p_val < 0.05, "Relationship is significant", "Relationship is not significant")

  # Shapiro-Wilk Test (Normality Check)
  shap_wilk_test <- shapiro.test(residuals(model))

  # Kolmogorov-Smirnov Test (Normality Check)
  kol_smir_test <- ks.test(residuals(model), "pnorm", mean = mean(residuals(model)), sd = sd(residuals(model)))

  # Residual Diagnostics Plot
  resid_hist <- hist(residuals(model), main = "Histogram of Residuals", xlab = "Residuals")

  # Q-Q Plot for Residuals
  qqnorm(residuals(model))
  qqline(residuals(model))

  # Homoscedasticity Check (Constant Variance)
  plot(model, residuals = TRUE, pch = 20, col = "blue")
  abline(h = 0, lty = 2, col = "red")

  # GAM Check
  gam_check_results <- gam.check(model)

  # Concurvity Check
  concurv <- concurvity(b = model, full = TRUE)


  # Return results
  return(list(
    model_summary = list(
      adj_r_squared = adj_r,
      p_value = p_val,
      significance = significance,
      shapiro_wilk = shap_wilk_test,
      kolmogorov_smirnov = kol_smir_test,
      concurvity = concurv,
      gam_check = gam_check_results
    )

  ))
}
