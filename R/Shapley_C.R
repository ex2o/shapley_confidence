library(speedglm)
library(RcppEigen)
library(compiler)
enableJIT(3)

d_ <- 4
n_ <- 1000
set.seed(0)
X_ <- matrix(rnorm(n_*d_),n_,d_)
y_ <- X_%*% (2*(0:(d_-1))) + rnorm(n_)

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

############################# Bottleneck 1
# library(profvis)
# library(foreach)
# library(doParallel)
# library(Rcpp)
# # Use parellel processing on all cores (except one)
# cores=detectCores()
# cl <- makeCluster(cores[1]-1)
# registerDoParallel(cl)
######## THIS IS NOT PARALLELIZED due to issues with Rcpp and foreach
v_ <- list()
system.time({
  ii_vec_ <- 1:n_terms_
  for (ii in ii_vec_) {
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
      c1 <- cov_det_Arma(corr_Z_, s1_, s2_, diff_mat_1_, diff_mat_2_)$output
      
      diff_mat_1_ <- det(as.matrix(corr_Z_[s1_[-1],s1_[-1]]))*solve(corr_Z_[s1_[-1],s1_[-1]])
      diff_mat_2_ <- det(as.matrix(corr_Z_[s2_[-1],s2_[-1]]))*solve(corr_Z_[s2_[-1],s2_[-1]])
      c2 <- cov_det_Arma(corr_Z_,s1_[-1],s2_[-1],diff_mat_1_,diff_mat_2_)$output
      
      diff_mat_1_ <- det(as.matrix(corr_Z_[s1_[-1],s1_[-1]]))*solve(corr_Z_[s1_[-1],s1_[-1]])
      diff_mat_2_ <- det(as.matrix(corr_Z_[s2_,s2_]))*solve(corr_Z_[s2_,s2_])
      c3 <- cov_det_Arma(corr_Z_,s1_[-1],s2_,diff_mat_1_,diff_mat_2_)$output
      
      diff_mat_1_ <- det(as.matrix(corr_Z_[s1_, s1_]))*solve(corr_Z_[s1_,s1_])
      diff_mat_2_ <- det(as.matrix(corr_Z_[s2_[-1],s2_[-1]]))*solve(corr_Z_[s2_[-1],s2_[-1]])
      c4 <- cov_det_Arma(corr_Z_,s1_,s2_[-1],diff_mat_1_,diff_mat_2_)$output
      
      value_ <- ( 1 / (a11 * a22) * c1 +
                    b11 * b22 / (a11^2 * a22^2) * c2 -
                    b11 / (a11^2 * a22) * c3 -
                    b22 / (a11 * a22^2) * c4 ) / n_
      cov_mat_[ii,jj] <- value_
      cov_mat_[jj,ii] <- value_
    }
  }
})
########################### End of bottleneck 1

############################# Bottleneck 2
shapley_var_ <- numeric(d_)
for (ii in 1:d_) {
  shapley_var_[ii] <- sum(diag(cov_mat_[which(belong_==ii),which(belong_==ii)])*
                            weights_[which(belong_==ii)]^2)
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
sqrt(shapley_var_)