#' Simulate cancer event time data with discrete-event simulation
#'
#' @param n Integer for cohort size.
#' @param l_params List of 4 vectors, params_Do, params_P, params_PC, and 
#' params_CD, each with named entries "min" and "max"
#'
#' @return A data table with estimated prevalence and confidence intervals at each age.
#' @import data.table
#'
cancer_des <- function(n, l_params) {
  # Initialize matrix of patient trajectories
  m_patients <- data.table(pt_id = 1:n)
  setkey(m_patients, pt_id)
  
  # Simulate time to death from other causes, preclinical cancer, clinical cancer, and death from cancer
  with(l_params, {
    m_patients[, `:=` (time_Do = runif(.N, min = params_Do["min"], max = params_Do["max"]),
                       time_P  = runif(.N, min = params_P["min"],  max = params_P["max"]),
                       time_PC = runif(.N, min = params_PC["min"], max = params_PC["max"]),
                       time_CD = runif(.N, min = params_CD["min"], max = params_CD["max"]))]
    
    # Calculate time from birth to clinical cancer and death from cancer
    m_patients[, `:=` (time_C  = time_P + time_PC,
                       time_Dc = time_P + time_PC + time_CD)]
    
    # Calculate time to death and cause of death
    m_patients[, `:=` (time_D = pmin(time_Do, time_Dc),
                       Dc = (time_Dc < time_Do))]
  })
  
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
      #Define differential equations
      dH  <- -(r_H_P + r_Do)*H
      dP  <- r_H_P*H - (r_P_C + r_Do)*P
      dC  <- r_P_C*P - (r_C_Dc + r_Do)*C
      dDo <- r_Do*(H + P + C)
      dDc <- r_C_Dc*C 
      dCInc <- r_P_C*P #Incidence of clinical cancer
      dx  <- c(dH, dP, dC, dDo, dDc, dCInc) #Combine results into a single vector dx
      list(dx)                    #Return result as a list
    }                             #NOTICE that here, we've assumed beta is constant
  )
}