#' Simulate cancer event time data with discrete-event simulation
#'
#' @param n Integer for cohort size.
#' @param l_params Nested list of lists with entries params_Do, params_P, params_PC, and 
#' params_CD, each of the format list(distr, params), where "distr" is a string with the
#' distribution name and params is a named list of parameters for the distribution
#'
#' @return A data table with estimated prevalence and confidence intervals at each age.
#' @import data.table
#'
cancer_des <- function(n, l_params) {
  # Initialize matrix of patient trajectories
  m_patients <- data.table(pt_id = 1:n)
  setkey(m_patients, pt_id)
  
  # Simulate time to death from other causes, preclinical cancer, clinical cancer, and death from cancer
  for (event in names(l_params)) {
    m_patients[, paste0("time_", event) := query_distr("r", .N, l_params[[event]][["distr"]], l_params[[event]][["params"]])]
  }
  
  # Calculate time from birth to clinical cancer and death from cancer
  m_patients[, `:=` (time_C  = time_P + time_PC,
                     time_Dc = time_P + time_PC + time_CD)]
  
  # Calculate time to death and cause of death
  m_patients[, `:=` (time_D = pmin(time_Do, time_Dc),
                     Dc = (time_Dc < time_Do))]
  
  return(m_patients)
}

# Differential equation cohort model for cancer
cancer_cohort_ode <- function(t, v_x, params) {
  H    <- v_x[1] # Healthy
  P    <- v_x[2] # Preclinical cancer
  C    <- v_x[3] # Clinical cancer
  Do   <- v_x[4] # Death from other causes
  Dc   <- v_x[5] # Death from cancer
  CInc <- v_x[6] # Incidence of clinical cancer
  
  with(                                               # We can simplify code using "with"
    as.list(params),                                  # This argument to "with" lets us use the variable names
    {                                                  
      # Define differential equations
      dH  <- -(r_P + r_Do)*H
      dP  <- r_P*H - (r_PC + r_Do)*P
      dC  <- r_PC*P - (r_CD + r_Do)*C
      dDo <- r_Do*(H + P + C)
      dDc <- r_CD*C 
      dCInc <- r_PC*P # Incidence of clinical cancer
      dx  <- c(dH, dP, dC, dDo, dDc, dCInc) # Combine results into a single vector dx
      list(dx)                    # Return result as a list
    }                            
  )
}

#' General function to query a distribution's density, cumulative distribution 
#' function, quantile function, or random generation function following the 
#' format in the R stats package
#'
#' \code{query_distr} is a flexible method to generate values from an inputted 
#' distribution
#'
#' @param target Type of output required, consistent with the R stats package: 
#' \code{d} for density, \code{p} for cumulative probability, \code{q} for 
#' quantile, and \code{r} for randomly generated value
#' @param x For \code{target} = \code{d} or \code{p}, the vector of quantiles 
#' at which to evaluate; for \code{target} = \code{q}, the vector of 
#' probabilities; for \code{target} = \code{r}, the number of observations to 
#' generate
#' @param distr String with the distribution name
#' @param params List of distribution parameters named correspondingly to the 
#' distribution's family of functions
#' @param ... List of optional arguments to the target distribution function
#' 
#' @return Target value from the distribution
#' 
#' @export
query_distr <- function(target, x, distr, params, ...) {
  # Get distribution name for function
  if(tolower(distr) == "exponential") {
    distr_name <- "exp"
  } else if(tolower(distr) == "binomial") {
    distr_name <- "binom"
  } else if(tolower(distr) == "uniform") {
    distr_name <- "unif"
  } else {
    distr_name <- tolower(distr)
  }
  
  # Call function based on target and distribution name if function exists
  function_name <- paste0(target, distr_name)
  if(is.function(match.fun(function_name))) {
    val <- do.call(match.fun(function_name), c(list(x), params, ...))
    return(val)
  } else {
    stop(paste0("Function ", function_name, " not found"))
  }
}