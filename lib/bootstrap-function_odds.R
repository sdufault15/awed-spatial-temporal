bootstrap_function_odds <- function(df, r1, r2, t, replace = TRUE){
  
  # Draw a sample with replacement
  temp <- df %>% 
    slice_sample(prop = 1,
                 replace = replace)
  
  ###############
  # MATRICES
  ###############
  # First, estimate pairwise distances in METERS
  dist.mat <- geodist(x = cbind(longitude = temp$longitude, latitude = temp$latitude))
  # Set up for difference in serotype and time (absolute time!)
  time.pairs <- expand.grid(temp$illness_onset, temp$illness_onset)
  time.dif <- as.Date(time.pairs[,1]) - as.Date(time.pairs[,2])
  time.mat <- abs(matrix(time.dif, nrow = nrow(temp)))
  
  # Serotype
  sero.pairs <- expand.grid(temp$serotype, temp$serotype)
  sero.related <- 1*(sero.pairs[,1] != "test-negative control")*(sero.pairs[,1] == sero.pairs[,2])
  sero.mat <- matrix(sero.related, nrow = nrow(temp))
  
  # Identify comparisons
  rownames(dist.mat) <- rownames(time.mat) <- rownames(sero.mat) <- temp$participant_id
  colnames(dist.mat) <- colnames(time.mat) <- colnames(sero.mat) <- temp$participant_id
  
  if (replace == FALSE){
    # Remove self-comparisons (observed data)
    diag(sero.mat) <- diag(time.mat) <- diag(dist.mat) <- NA 
  } else {
    # Remove self-comparisons (boostrapped data)  
    for (i in 1:nrow(sero.mat)){
      for (j in 1:ncol(sero.mat)){
        if (rownames(sero.mat)[i] == colnames(sero.mat)[j]){
          sero.mat[i,j] <- time.mat[i,j] <- dist.mat[i,j] <- NA
        } 
      }
    }
  }
  
  
  
  
  ###############
  # ESTIMATION
  ###############
  tau_overall <- map2_dfr(
    r1, # x 
    r2, # y
    ~bin.fun(r1 = .x, 
             r2 = .y,
             t = t,
             A = time.mat,
             B = dist.mat,
             C = sero.mat))
  return(tau_overall)
}
