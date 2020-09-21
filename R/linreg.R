#' Linear regression model
#'
#' \code{linreg} returns a linear regression model for formula and data inputs
#'
#' This is a function that builds a linear regression model using ordinary least squares method and RC class.
#'
#' @param formula A formula with dependent numeric variable on left and independent numeric variable on right.
#' @param data A dataset 
#' @export


library(ggplot2)
library(dplyr)

linreg<-setRefClass("linreg", 
        fields= list(
          formula="formula", 
          data="data.frame",
          call = "vector",
          coefficients_est = "matrix",
          fitted_vals = "matrix",
          residuals_est = "matrix",
          degrees_of_freedom = "numeric",
          res_var_est = "matrix",
          v_coeff = "vector",
          find_t_vals = "numeric",
          p_vals = "vector",
          std_res = "matrix"), 
        
        methods=list(
          initialize = function(formula,data) {
            call <<- c("linreg(formula = ",
                       Reduce(paste,deparse(formula)),
                       ", data = ",
                       deparse(substitute(data)),
                       ")")
            formula <<- formula
            data <<- data
            
            #Independent X-values
            X <- model.matrix(formula, data)
            #Dependent y-values
            y_label <- all.vars(formula)[1]
            y <- data[[y_label]]
            
            #Estimating the beta-coefficients
            coefficients_est <<- solve(t(X) %*% X) %*% t(X) %*% y  
            
            #Estimating the fitted values of the model  
            fitted_vals <<- X %*% coefficients_est
            
            #Estimating the residuals
            residuals_est <<- y - fitted_vals
            
            #Obervations & parameters
            n <- nrow(X)
            p <- ncol(X)
            degrees_of_freedom <<- n - p
            
            #Estimation of residual variance (0,1851, i.e (residuals stand. error)^2)
            res_var_est <<- (t(residuals_est) %*% residuals_est) / degrees_of_freedom
            
            
            #Estimating the variance of the beta coeff.  
            var_of_coeff <- (solve(t(X)%*%X))*(as.numeric(res_var_est))
            
            #Extracting the variance diagonal of the cov. matrix (=coeff. Std.Error^2)
            v_coeff<<-diag(var_of_coeff)
            
            #Finding the t-values
            
            find_t_vals <<- as.numeric(coefficients_est) / sqrt(v_coeff)

            
            p_vals <<- 2*pt(-abs(find_t_vals), degrees_of_freedom)
            
            #Estimating the standardized residuals
            
            hat_values <- X %*% solve(t(X) %*% X) %*% t(X)
            
            std_res <<- residuals_est/(as.numeric(res_var_est) * sqrt(1-diag(hat_values)))
            
          },
            
          print = function(){
            #print call
            cat("\nCall:\n")
            lapply(call, cat)
            #Names of the coefficients
            c_names <- row.names(coefficients_est)
            #Vals of coeff.
            c_vals <- as.numeric(coefficients_est)
            #print Coefficients:
            cat("\n\nCoefficients:\n")
            cat(paste(c_names,collapse = "  "),collapse="\n")
            cat(paste(c_vals,collapse = "  "),collapse="\n")
          },
          plot = function(){
            
            plotting_data <- data.frame(p_fitted <- fitted_vals, p_residuals <- residuals_est, p_std_res <- sqrt(abs(std_res)))
            
            residuals_v_fitted <- ggplot(data=plotting_data, aes(x=p_fitted, y=p_residuals))+
              geom_point(shape = 1) + 
              geom_smooth(method = "loess", formula= y~x) + 
              geom_hline(yintercept=0, col="red", linetype="dashed") +
              xlab("Fitted values")+
              ylab("Residuals") + 
              ggtitle("Residual vs Fitted") + 
              theme(plot.title = element_text(hjust = 0.5))
            
            
            scale_location <- ggplot(data=plotting_data, aes(x=p_fitted, y=p_std_res)) +
              geom_point(shape = 1) + 
              geom_smooth(method="loess",formula = y~x) + 
              xlab("Fitted values") + 
              ggtitle("Scale - Location") +
              ylab(expression(sqrt("|Standardized residuals|"))) + 
              theme(plot.title = element_text(hjust = 0.5))
            
            return(list(residuals_v_fitted, scale_location)) 
            
          },
          resid = function(){
            return(residuals_est)
          },
          pred = function(){
            return(fitted_vals)
          },
          coef = function(){
            
           coef_vector <- as.vector(t(coefficients_est))
            names(coef_vector) <- row.names(coefficients_est)
            
            return(coef_vector)
            
          },
          summary = function(){
            
            #print call
            cat("\nCall:\n")
            lapply(call, cat)
            
            #print Residuals:
            cat("\n\nCoefficients:\n")
            
            summary_vec <- cbind(coefficients_est,as.numeric(lapply(v_coeff, sqrt)), find_t_vals, p_vals)
            colnames(summary_vec) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
            
            printCoefmat(summary_vec, P.value=TRUE, has.Pvalue=TRUE)
            
            #Print the standard error
            s_res <- c("--- \nResidual standard error: ",sqrt(res_var_est)," on ",degrees_of_freedom," degrees of freedom")
            lapply(s_res, function(x) cat(x))
            cat("\n")
          }
        )
)


  
  
    
  
                      
                      
  
  


