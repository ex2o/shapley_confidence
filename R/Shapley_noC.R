library(speedglm)
library(RcppEigen)
library(compiler)
library(profvis)
library(foreach)
library(doParallel)
library(Rcpp)
#enableJIT(3)

shapley <- function(y_, X_) {
  
  cov_det_R <- function(corr_mat_, s1_, s2_, diff_mat_1_, diff_mat_2_, kappa) {
    output <- 0
    wii <- 0
    for (ii in s1_) {
      wii <- wii + 1 #which(s1_==ii); cat("wii: ", wii, "\n")
      wjj <- 0
      for (jj in s1_) {
        wjj <- wjj + 1 #which(s1_==jj); cat("wjj: ", wjj, "\n")
        dm_ij <- diff_mat_1_[wii, wjj]
        cm_ij <- corr_mat_[ii,jj]
        cm_ji <- corr_mat_[jj,ii]
        wll <- 0
        for (ll in s2_) {
          wll <- wll + 1 #which(s2_==ll); cat("wll: ", wll, "\n")
          cm_il <- corr_mat_[ii,ll]
          cm_jl <- corr_mat_[jj,ll]
          cm_li <- corr_mat_[ll,ii]
          cm_lj <- corr_mat_[ll,jj]
          wmm <- 0
          for (mm in s2_) {
            wmm <- wmm + 1 #which(s2_==mm); cat("wmm: ", wmm, "\n")
            dm_lm <- diff_mat_2_[wll,wmm]
            cm_im <- corr_mat_[ii,mm]
            cm_jm <- corr_mat_[jj,mm]
            cm_lm <- corr_mat_[ll,mm]
            cm_mi <- corr_mat_[mm,ii]
            cm_mj <- corr_mat_[mm,jj]
            cm_ml <- corr_mat_[mm,ll]
            
            output <- output +
              dm_ij * dm_lm * (
                0.5 * cm_ij * cm_lm *
                  ( cm_il^2 + cm_im^2 + cm_jl^2 + cm_jm^2 ) +
                  cm_il * cm_jm +
                  cm_im * cm_jl -
                  cm_ij * cm_il * cm_im -
                  cm_ji * cm_jl * cm_jm -
                  cm_li * cm_lj * cm_lm -
                  cm_mi * cm_mj * cm_ml ) * kappa
          }
        }
      }
    }
    
    return(output)
  }
    
  n_ <- nrow(X_)
  d_ <- ncol(X_)
  
  count_ <- 0
  R2_ <- c()
  weights_ <- c()
  sign_ <- c()
  belong_ <- c()
  covar_list_ <- list()
  for (ii in 0:(d_-1)) {
    comb_ <- combn(d_-1,ii)
    for (jj in 1:d_) {
      indep_ <- (1:d_)[-jj]
      if (ii == 0) {
        count_ <- count_ + 1
        weights_[count_] <- factorial(ii)*factorial(d_-ii-1)/factorial(d_)
        sign_[count_] <- 1
        belong_[count_] <- jj
        R2_[count_] <- summary(speedlm(y_~X_[,jj]))$r.squared
        covar_list_[[count_]] <- jj
      } else {
        for (kk in 1:ncol(comb_)) {
          count_ <- count_ + 1
          weights_[count_] <- factorial(ii)*factorial(d_-ii-1)/factorial(d_)
          sign_[count_] <- 1
          belong_[count_] <- jj
          R2_[count_] <- summary(speedlm(y_~X_[,c(jj,indep_[comb_[,kk]])]))$r.squared
          covar_list_[[count_]] <- c(jj,indep_[comb_[,kk]])
  
          count_ <- count_ + 1
          weights_[count_] <- factorial(ii)*factorial(d_-ii-1)/factorial(d_)
          sign_[count_] <- -1
          belong_[count_] <- jj
          R2_[count_] <- summary(speedlm(y_~X_[,c(indep_[comb_[,kk]])]))$r.squared
          covar_list_[[count_]] <- c(indep_[comb_[,kk]])
        }
      }
    }
  }
  
  shapley_ <- numeric(d_)
  for (ii in 1:d_) {
    shapley_[ii] <- sum(R2_[which(belong_==ii)]*weights_[which(belong_==ii)]*sign_[which(belong_==ii)])
  }
  
  shapley_
  sum(shapley_)
  summary(speedlm(y_~X_))$r.squared
  
  Z_ <- cbind(y_,X_)
  corr_Z_ <- cor(Z_)
  n_terms_ <- length(sign_)
  cov_mat_ <- matrix(0,n_terms_,n_terms_)
  diag(cov_mat_) <- 4*R2_*(1-R2_)^2/n_
  
  
  # Calculate kappa
  muZ <- apply(Z_, MARGIN = 2, FUN = mean)
  invCov <- solve(cov(Z_))
  mahalanobis_sq <- apply(Z_, MARGIN = 1,
      FUN = function(z){ t(z-muZ) %*% invCov %*% (z-muZ) })
  kappa <- sum( mahalanobis_sq^2 )/(n_*(d_+1)*(d_+3))
  
  ############################# Bottleneck 1
  # Use parellel processing on all cores (except one)
  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  system.time({
    ii_vec_ <- 1:n_terms_
    l_ <- foreach(ii = ii_vec_) %dopar% {
      s1_ <- c(1, 1+covar_list_[[ii]])
      a11 <- det(as.matrix(corr_Z_[s1_[-1],s1_[-1]]))
      b11 <- det(as.matrix(corr_Z_[s1_,s1_]))
      diff_mat_1_ <- det(as.matrix(corr_Z_[s1_,s1_]))*solve(corr_Z_[s1_,s1_])
      jj_vec_ <- ii_vec_[ii > ii_vec_]
      v_ <- vector(mode = "numeric", length = length(ii_vec_))
      for (jj in jj_vec_) {
        s2_ <- c(1,1+covar_list_[[jj]])
        diff_mat_2_ <- det(as.matrix(corr_Z_[s2_,s2_]))*solve(corr_Z_[s2_,s2_])
        a22 <- det(as.matrix(corr_Z_[s2_[-1],s2_[-1]]))
        b22 <- det(as.matrix(corr_Z_[s2_,s2_]))
  
        diff_mat_1_ <- det(as.matrix(corr_Z_[s1_,s1_]))*solve(corr_Z_[s1_,s1_])
        diff_mat_2_ <- det(as.matrix(corr_Z_[s2_,s2_]))*solve(corr_Z_[s2_,s2_])
        c1 <- cov_det_R(corr_Z_, s1_, s2_, diff_mat_1_, diff_mat_2_, kappa)
  
        diff_mat_1_ <- det(as.matrix(corr_Z_[s1_[-1],s1_[-1]]))*solve(corr_Z_[s1_[-1],s1_[-1]])
        diff_mat_2_ <- det(as.matrix(corr_Z_[s2_[-1],s2_[-1]]))*solve(corr_Z_[s2_[-1],s2_[-1]])
        c2 <- cov_det_R(corr_Z_,s1_[-1],s2_[-1],diff_mat_1_,diff_mat_2_, kappa)
  
        diff_mat_1_ <- det(as.matrix(corr_Z_[s1_[-1],s1_[-1]]))*solve(corr_Z_[s1_[-1],s1_[-1]])
        diff_mat_2_ <- det(as.matrix(corr_Z_[s2_,s2_]))*solve(corr_Z_[s2_,s2_])
        c3 <- cov_det_R(corr_Z_,s1_[-1],s2_,diff_mat_1_,diff_mat_2_, kappa)
  
        diff_mat_1_ <- det(as.matrix(corr_Z_[s1_, s1_]))*solve(corr_Z_[s1_,s1_])
        diff_mat_2_ <- det(as.matrix(corr_Z_[s2_[-1],s2_[-1]]))*solve(corr_Z_[s2_[-1],s2_[-1]])
        c4 <- cov_det_R(corr_Z_,s1_,s2_[-1],diff_mat_1_,diff_mat_2_, kappa)
  
        v_[jj] <- ( 1 / (a11 * a22) * c1 +
                      b11 * b22 / (a11^2 * a22^2) * c2 -
                      b11 / (a11^2 * a22) * c3 -
                      b22 / (a11 * a22^2) * c4 ) / n_
      }
      v_
    }
    # Unpack the list into the covariance matrix
    for ( ii in ii_vec_ ) {
      jj_vec_ <- ii_vec_[ii > ii_vec_]
      for ( jj in jj_vec_ ) {
        value_ <- l_[[ii]][jj]
        cov_mat_[ii,jj] <- value_
        cov_mat_[jj,ii] <- value_
      }
    }
    stopCluster(cl)
  })
  ########################### End of bottleneck 1
  
  ############################# Bottleneck 2
  shapley_var_ <- numeric(d_)
  for (ii in 1:d_) {
    bii_ = (belong_==ii)
    shapley_var_[ii] <- sum(diag(cov_mat_[bii_,bii_])*
                              weights_[bii_]^2)
    for (aa in 1:sum(belong_==ii)) {
      for (bb in 1:sum(belong_==ii)) {
        if (aa > bb) {
          shapley_var_[ii] <- shapley_var_[ii] +
            2*cov_mat_[which(belong_==ii),which(belong_==ii)][aa,bb]*
            sign_[which(belong_==ii)[aa]]*
            sign_[which(belong_==ii)[bb]]*
            weights_[which(belong_==ii)[aa]]*
            weights_[which(belong_==ii)[bb]]
        }
      }
    }
  }
  ########################### End of bottleneck 2
  shapley_
  sqrt(abs(shapley_var_))
  return(list(shapley = shapley_, var = shapley_var_))
}


# ######## UNUSED TESTING CODE
# summary(speedlm(y_~X_))$r.squared
# t(y_) %*% X_ %*% solve(t(X_) %*% X_) %*% t(X_) %*% y_ %*% solve(t(y_) %*% y_)
# summary(speedlm(y_~X_, method = "qr"))$r.squared
# diag(cov_mat_)
# cov_mat_[1:10,1:3]
