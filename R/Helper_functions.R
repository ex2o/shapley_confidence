cov_det_R <- function(corr_mat_, s1_, s2_, diff_mat_1_, diff_mat_2_) {
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
                cm_mi * cm_mj * cm_ml )
        }
      }
    }
  }
  
  return(output)
}
