

temp <- lapply(c("randomForest", "gbm", "dismo", "quadprog", "quantreg"), library, character.only=T)



TT <- 1000

set.seed(10)

rtf <- function(a1, a2) { sin(a1+a2)/(a1+a2) }

df <- data.frame(x1=(runif(TT, min=1, max=6)), x2=(runif(TT, min=1, max=6)))

df$m <- rtf(df$x1, df$x2)

df$y <- df$m + rnorm(TT, sd=.1)

model_lm <- lm(y~x1+x2, data=df)

p_lm <- function(x1,x2) predict(model_lm, newdata=data.frame(x1=x1,x2=x2))

model_rf <- randomForest(y~x1+x2, data=df)

p_rf <- function(x1,x2) as.numeric(predict(model_rf, newdata=data.frame(x1=x1,x2=x2), type="response"))

model_gbm <- gbm.step(data=df, gbm.x = 1:2, gbm.y = 4, plot.main = F,
                      
                      family = "gaussian", tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)

p_boost <- function(x1,x2) predict(model_gbm, newdata=data.frame(x1=x1,x2=x2), n.trees=1200)


TT_out_of_sample <- 500

df_new <- data.frame(x1=(runif(TT_out_of_sample, min=1, max=6)), x2=(runif(TT_out_of_sample, min=1, max=6)))

df_new$m <- rtf(df_new$x1, df_new$x2)

df_new$y <- df_new$m+rnorm(TT_out_of_sample,sd=.1)

boost_fhat <-p_boost(df_new$x1, df_new$x2)

rf_fhat <- p_rf(df_new$x1, df_new$x2)

lm_fhat <- p_lm(df_new$x1, df_new$x2)


mat_fhat <- cbind(boost_fhat, rf_fhat, lm_fhat)



colnames(mat_fhat) <- c("Boosting", "Random Forest", "OLS")



FA_schemes <- c("simple", "ols", "robust", "variance based", "cls")



results <- list()



init_window = 30

for(i in 1:length(FA_schemes)){
  
  results[[i]] <- FA(df_new$y, mat_fhat, Averaging_scheme= FA_schemes[i], init_window = init_window )
  
  cat( paste(FA_schemes[i],": ", 100*as.numeric(format( results[[i]]$MSE, digits=4) ), "\n") )
  
}




mat_err <- apply(mat_fhat, 2, function(x) (df_new$y - x)^2 )

apply(mat_err[init_window:TT_out_of_sample, ], 2, function(x) 100*( mean(x) ) )








FA <- function(obs, mat_fhat, init_window=30,
               
               Averaging_scheme=c("simple", "ols", "robust", "variance based", "cls", "best")) {
  
  pckg = c("quantreg", "quadprog")
  
  temp <- unlist(lapply(pckg, require, character.only=T))
  
  if (!all(temp==1) ) {
    
    stop("This function relies on packages \"quadprog\" and \"quantreg\".

Use ?install.packages if they are not yet installed. \n")
    
  }
  
  mat_err <- apply(mat_fhat, 2, function(x) obs - x)
  
  TT <- NROW(mat_fhat)
  
  p <- NCOL(mat_fhat)
  
  pred <- NULL
  
  weights <- matrix(ncol=p, nrow = TT)
  
  ## Subroutine needed:
  
  sq_er <- function(obs, pred){ mean( (obs - pred)^2 )}
  
  
  
  if(length(Averaging_scheme) > 1) {stop("Pick only one of the following:

c(\"simple\", \"ols\", \"robust\", \"variance based\", \"cls\", \"best\")") }
  

  
  if(Averaging_scheme=="simple") {
    
    pred <- apply(mat_fhat, 1, mean)
    
    weights <- matrix( 1/p, nrow = TT, ncol = p)
    
    ## OLS weights
    
  } else if (Averaging_scheme=="ols") {
    
    weights <- matrix(ncol=(p+1), nrow = TT)
    
    for (i in init_window:TT) {
      
      weights[i,] <- lm(obs[1:(i-1)]~ mat_fhat[1:(i-1),])$coef
      
      pred[i] <- t(weights[i,])%*%c(1, mat_fhat[i,])
      
      ## Robust weights
      
    } } else if (Averaging_scheme=="robust") {
      
      weights <- matrix(ncol=(p+1), nrow = TT)
      
      for (i in init_window:TT) {
        
        weights[i,] <- rq(obs[1:(i-1)]~ mat_fhat[1:(i-1),])$coef
        
        pred[i] <- t(weights[i,])%*%c(1, mat_fhat[i,])
        

      } } else if (Averaging_scheme=="variance based") {
        
        for (i in init_window:TT) {
          
          temp = apply(mat_err[1:(i-1),]^2,2,mean)/sum(apply(mat_err[1:(i-1),]^2,2,mean))
          
          weights[i,] <- (1/temp)/sum(1/temp)
          
          pred[i] <- t(weights[i,])%*%c(mat_fhat[i,])
          

        }} else if (Averaging_scheme== "cls") {
          

          cls1 = function(y, predictions){
            
            Rinv <- solve(chol(t(predictions) %*% predictions))
            
            C <- cbind(rep(1,NCOL(predictions)), diag(NCOL(predictions)))
            
            b = c(1, rep(0,NCOL(predictions)))
            
            d = t(y) %*% predictions
            
            qp1 = solve.QP(Dmat= Rinv, factorized= TRUE, dvec= d, Amat= C, bvec = b, meq = 1)
            
            weights = qp1$sol
            
            yhatcls = t(weights%*%t(predictions))
            
            list(yhat= yhatcls, weights= weights)
            
          }
          
          ##
          
          for (i in init_window:TT) {
            
            weights[i,] <- cls1(obs[1:(i-1)], mat_fhat[1:(i-1),])$weights
            
            pred[i] <- t(weights[i,])%*%c(mat_fhat[i,])
            

          } } else if (Averaging_scheme== "best") {
            
            temp <- apply(mat_fhat[-c(1:init_window),], 2, sq_er, obs= obs[-c(1:init_window)])
            
            for (i in init_window:TT) {
              
              weights[i,] <- rep(0,p)
              
              weights[i, which.min(temp)] <- 1
              
              pred[i] <- t(weights[i,])%*%c(mat_fhat[i,])
              
            } }
  
  MSE <- sq_er(obs= obs[-c(1:init_window)], pred= pred[-c(1:init_window)])
  
  return(list(prediction= pred, weights = weights, MSE= MSE))
  
}

