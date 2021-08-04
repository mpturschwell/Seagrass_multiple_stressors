#----------------------------------------------
# Calculate LRR from multi-stressor models
#----------------------------------------------

#--------
# Inputs 
#--------
# X = dataframe of predicted output along gradients of Stressor 1 and Stressor 2
# S1 = length of values tested for stressor 1 (e.g. Temp)
# S2 = length of values tested for stressor 2 (e.g. Light)
# E.g. LRRX <- Log_response_ratio(X = X, S1 = length(tempsens), S2 = length(lightsens))

Log_response_ratio <- function(X, S1, S2){
              LRRX <- matrix(nrow = S1, ncol = S2)
              icontrol <- which(X == max(X, na.rm = TRUE), arr.ind = TRUE)[1]
              jcontrol <- which(X == max(X, na.rm = TRUE), arr.ind = TRUE)[2]
          for (i in 1:S1){
            for (j in 1:S2){
              LRRX[i,j] <- log((X[i,j] - X[icontrol,jcontrol])/
                               (X[i,jcontrol] + X[icontrol,j] - 2*X[icontrol,jcontrol]))
                            }
                       }
              LRRX[icontrol, jcontrol] <- 0 #control will be NaN so make it zero 
              return(LRRX)
}

#This one uses a fixed control scenario
Log_response_ratio_alt <- function(X, S1, S2, icontrol, jcontrol){
  LRRX <- matrix(nrow = S1, ncol = S2)
  for (i in 1:S1){
    for (j in 1:S2){
      LRRX[i,j] <- log(abs((X[i,j] - X[icontrol,jcontrol])/   # added in abs here to deal with negative logs which make NaN values
                         (X[i,jcontrol] + X[icontrol,j] - 2*X[icontrol,jcontrol])))
      }
  }
  LRRX[icontrol, jcontrol] <- 0 #control will be NaN so make it zero 
  return(LRRX)
}






