#' Building a linear regression model using ordinary least squares method
#' 
#' @param formula object is a formula with dependent numeric variable on left and independent numeric variable on right.
#' @param data A numeric matrix.
#' @return Statistics for linear regression model.
#' 
#' 

linreg <- function(formula,data){
    
    #Independent X-values
    X <- model.matrix(formula, data)
    #Dependent y-values
    y_label <- all.vars(as.formula(formula))[1]
    y <- (data[[y_label]])
    
    #Obervations & parameters
    n <- nrow(X)
    p <- ncol(X)
    
    #Estimating the beta-coefficients
    coefficients_est <- solve(t(X) %*% X) %*% t(X) %*% y  
    
    #Estimating the fitted values of the model  
    fitted_vals <- X%*%coefficients_est
    
    #Estimating the residuals
    residuals_est <- y - X%*%coefficients_est
    
    degrees_of_freedom <- n - p
    
    #Estimation of residual variance (0,1851)
    res_var_est <- t(residuals_est)%*%residuals_est/degrees_of_freedom
    res_var_est <- drop(res_var_est) #Drops res_v_e from 1x1 to a scalar
    
    #Estimating the variance of the beta coeff.  
    var_of_coeff <- (solve(t(X)%*%X))*(res_var_est)
    
    #Extracting the variance diagonal of the cov. matrix (=Std.Error^2)
    v_coeff<-diag(var_of_coeff)
    
    #Finding the t-values
    find_t_vals <- 0
    
    for(i in 1:p){  
      
      find_t_vals[i] <-  coefficients_est[i]/sqrt(abs(v_coeff[i]))
      
    }
    
    p_vals <- pt(find_t_vals, degrees_of_freedom)
    
    return(coefficients_est)
    
  }

