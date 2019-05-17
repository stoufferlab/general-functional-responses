# Given a specified model and an ordered vector of parameter values to be passed to the likelihood function, assign the parameter values to the appropriate parameter names.
# Note: Most parameter values are exponentiated (to avoid negative values during fitting).

set_params<-function(params,
                     model=c('Holling.I',
                             'Holling.II',
                             'Beddington.DeAngelis',
                             'Crowley.Martin',
                             'Stouffer.Novak.I',
                             'Stouffer.Novak.II',
                             'Stouffer.Novak.III',
                             'Ratio',
                             'Arditi.Ginzburg',
                             'Hassell.Varley',
                             'Arditi.Akcakaya')
                     ){
  
  model <- match.arg(model)
  
  if(is.null(params)){
    stop("Must pass 'params' to set_params()")
  }
  
  if(model == "Holling.I"){
    assign('attack',       exp(params[1]), envir = .GlobalEnv)
    assign('handling',     0, envir = .GlobalEnv)
    assign('interference', 0, envir = .GlobalEnv)
    assign('phi_numer',    1, envir = .GlobalEnv)
    assign('phi_denom',    1, envir = .GlobalEnv)
    
  } else  if(model == "Holling.II"){
      assign('attack',       exp(params[1]), envir = .GlobalEnv)
      assign('handling',     exp(params[2]), envir = .GlobalEnv)
      assign('interference', 0, envir = .GlobalEnv)
      assign('phi_numer',    1, envir = .GlobalEnv)
      assign('phi_denom',    1, envir = .GlobalEnv)
      
  } else  if(model == "Beddington.DeAngelis"){
      assign('attack',       exp(params[1]), envir = .GlobalEnv)
      assign('handling',     exp(params[2]), envir = .GlobalEnv)
      assign('interference', exp(params[3]), envir = .GlobalEnv)
      assign('phi_numer',    1, envir = .GlobalEnv)
      assign('phi_denom',    1, envir = .GlobalEnv)
      
  } else  if(model == "Crowley.Martin"){
      assign('attack',       exp(params[1]), envir = .GlobalEnv)
      assign('handling',     exp(params[2]), envir = .GlobalEnv)
      assign('interference', exp(params[3]), envir = .GlobalEnv)
      assign('phi_numer',    1, envir = .GlobalEnv)
      assign('phi_denom',    0, envir = .GlobalEnv)
      
  } else  if(model == "Stouffer.Novak.I"){
      assign('attack',       exp(params[1]), envir = .GlobalEnv)
      assign('handling',     exp(params[2]), envir = .GlobalEnv)
      assign('interference', exp(params[3]), envir = .GlobalEnv)
      assign('phi_numer',    1, envir = .GlobalEnv)
      assign('phi_denom',    params[4], envir = .GlobalEnv)
      
  } else  if(model == "Stouffer.Novak.II"){
      assign('attack',       exp(params[1]), envir = .GlobalEnv)
      assign('handling',     exp(params[2]), envir = .GlobalEnv)
      assign('interference', exp(params[3]), envir = .GlobalEnv)
      assign('phi_numer',    params[4], envir = .GlobalEnv)
      assign('phi_denom',    1, envir = .GlobalEnv)
      
  } else  if(model == "Stouffer.Novak.III"){
      assign('attack',       exp(params[1]), envir = .GlobalEnv)
      assign('handling',     exp(params[2]), envir = .GlobalEnv)
      assign('interference', exp(params[3]), envir = .GlobalEnv)
      assign('phi_numer',    params[4], envir = .GlobalEnv)
      assign('phi_denom',    params[5], envir = .GlobalEnv)
      
  } else  if(model == "Ratio"){
      assign('attack',   exp(params[1]), envir = .GlobalEnv)
      assign('handling', 0, envir = .GlobalEnv)
      assign('exponent', 1, envir = .GlobalEnv)
      
  } else  if(model == "Arditi.Ginzburg"){
      assign('attack',   exp(params[1]), envir = .GlobalEnv)
      assign('handling', exp(params[2]), envir = .GlobalEnv)
      assign('exponent', 1, envir = .GlobalEnv)
      
  } else  if(model == "Hassell.Varley"){
      assign('attack',   exp(params[1]), envir = .GlobalEnv)
      assign('handling', 0, envir = .GlobalEnv)
      assign('exponent', exp(params[2]), envir = .GlobalEnv)
      
  } else  if(model == "Arditi.Akcakaya"){
      assign('attack',   exp(params[1]), envir = .GlobalEnv)
      assign('handling', exp(params[2]), envir = .GlobalEnv)
      assign('exponent', exp(params[3]), envir = .GlobalEnv)
      
  } else  stop("Model not correctly specified in set_params()")
}