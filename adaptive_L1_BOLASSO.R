library(gamlr)


# with great assistance from https://scholar.harvard.edu/files/emily/files/section-5-bootstrap-lasso-ridge.r_01.txt
# and matthew taddy's gamlr package (https://cran.r-project.org/web/packages/gamlr/gamlr.pdf)


# Bootstrapped gamma LASSO, helper function that will be used in the two other functions

# Inputs: x: the design matrix (covariates)
# y: a vector of the dependent variable
# data: the dataset
# nrep: number of bootstrap replicates
# free_indices: the indices revealing which coefficients will be unpenalized. see gamlr for more detail
# gamma_vec: the values of gamma we want to cross validate over
# ... other keyword parameters to be passed into the gamlr function.

# Output: a coefficient matrix with each row representing the value of the coefficients for one bootstrap replicate

gamma_BOLASSO <- function(x, y, data, nrep = 100, free_indices = NULL, gamma_vec = c(0,1,10), ...) {
  
  # initialize the coefficient matrix with the proper specifications. each row is one bootstrap replicate
  coef_mtx <- matrix(0, nrow = nrep, ncol = ncol(x) + 1)
  
  for (i in 1:nrep) {
    
    lowest_error <- Inf
    best_model <- NULL
    
    # sample some data from the sample randomnly
    bootstrap_id <- sample(seq(nrow(data)), nrow(data), replace = TRUE)
    bootstrap_data <- data[bootstrap_id,]

    for (j in 1:length(gamma_vec)) {
      
      # we run the gamlr regression and collect the best lambda that receives the highest score
      cv.obj <- cv.gamlr(x, y, nfold = 5, family = "gaussian", gamma = gamma_vec[j], free = free_indices, nlambda = 100)
      
      # we cross validate again, but this time over gamma
      if (cv.obj$cvm[which.min(cv.obj$seg.min)] < lowest_error) {
        best_model <- cv.obj
        lowest_error <- cv.obj$cvm[which.min(cv.obj$seg.min)]
      }

    }
    
    # select the coefficients from the bootstrap model
    coef_mtx[i,] <- as.numeric(coef(best_model, select = "min"))
  }
  return(coef_mtx)
}



# A two stage BOLASSO. In this case the BOLASSO acts as a variable selection method, and the variables which are selected are then run in a second-stage regression
# (please see the readme for more details)

# Inputs: x: the design matrix (covariates)
# y: a vector consisting of the dependent variable
# dep_var: a string revealing the name of the dependent variable
# data: the dataset
# nrep: number of bootstrap replicates
# free_indices: the indices revealing which coefficients will be unpenalized. see gamlr for more detail
# gamma_vec: the values of gamma we want to cross validate over
# thresh: a real number between 0 and 1 that we will use to include in our final (second-stage) regression
# finalmodel_family: the function specifying the error function and link function of the second stage/post-LASSSO regression to be passed into the glm function (https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/glm)
# default is Gaussian so default output will be a linear regression. 
# ... other keyword parameters to be passed into the gamlr function.

# Returns a glm regression object that is run on a subset of variables 

two_stage_gamma_BOLASSO <- function(x, y, dep_var, coef_mtx, data, nrep = 100, free_indices = NULL, gamma_vec = c(0,1,10), thresh = 0.9,  
                              finalModel_family = gaussian(link = identity), ...) {
  
  # we use the gamma BOLASSO helper function to extract the coefficient matrix
  # coef_mtx <- gamma_BOLASSO(x, y, data, nrep, free_indices, gamma_vec, ...)
  
  # return a list of variables' indices which have >= thresh of their bootstrapped values not equal to 0
  b <- apply(coef_mtx, MARGIN = 2, FUN = function(x) length(which(x != 0))) >= thresh
  
  # grab the names of the non-zero variables
  vars <- colnames(data[,b])
  # (and remove the dependent variable from this list)
  vars <- vars[which(vars != dep_var)]
  
  # we create an R regression formula
  vars_form <- paste(vars, collapse = " + ")
  reg_formula <- as.formula(paste(dep_var, " ~", vars_form))
  
  # the second stage post-LASSO regression with the non-zero variables
  regression <- glm(reg_formula, data=data, family = finalModel_family)
  
  regression
}




# The "ensemble" BOLASSO. This time, we do not use the LASSO as a variable selection method, but we average over all the bootstrap replicates to produce a final result

# Inputs: x: the design matrix (covariates)
# y: a vector of the dependent variables
# data: the dataset used
# nrep: number of bootstrap replicates
# free_indices: the indices revealing which coefficients will be unpenalized. see gamlr for more detail
# gamma_vec: the values of gamma we want to cross validate over
# ... other arguments that would be passed into glmnet (NOT for glm)

# output: a vector corresponding to an average of the NON-ZERO coefficients over nrep bootstrap replicates

ensemble_gamma_BOLASSO <- function(x, y, data, nrep = 100, free_indices = NULL, gamma_vec = c(0,1,10), ...) {
  
  # we call our BOLASSO helper function which will return the value of the coefficients for nrep bootstrap replicates
  coef_mtx <- gamma_BOLASSO(x, y, data, nrep, free_indices, gamma_vec, ...)
  
  # we collect the average of the 100 bootstrap estimates
  coef_vec <- apply(coef_mtx, MARGIN = 2, mean)
  
  # we extract the indices of the coefficients which are not zero
  indices <- which(coef_vec != 0)
  
  # we collect those non-zero coefficents
  coef_vec <- coef_vec[indices]
  
  # update the names of the vector
  names(coef_vec) <- colnames(data[, indices])
    
  coef_vec
}
