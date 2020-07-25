library(tidyverse)

result_path <- function(type, d_) {
  paste0("../Julia/results/",d_,"/sim",type,"/")
}

# Load results for all n, for given correlation c_
read_results <- function(type, d_, c_, bs = F, small = F) {
  c_bs <- c(0,0.1,0.2,0.3,0.6,0.9,0.99)
  c_all <- if (bs) {c_bs} else {c(seq(0,0.9,0.1), 0.99)}
  type  <- if (small) {paste0(type, "s")} else {type}
  type  <- if (bs) {paste0(type,"_b")} else {type}
  n_max <- if (bs) {2000} else {5000}
  n_all <- if (small) {seq(5,100,5)} else {seq(100, n_max, 100)}
  c_index <- pmatch(c_, c_all)
  path <- result_path(type, d_)
  columns <- c("coverage", "covrgCI_L", "covrgCI_U",
               "width_mu", "width_sd", "widthCI_L", "widthCI_U")
  result <- list()
  for (col in columns) {
    result[[col]] <- read.csv(paste0(path,col,".csv"), header = F)[,c_index]
  }
  result[["n"]] <- n_all
  return(tibble::as_tibble(result))
}

# Create tibble with results for all c, for chosen type A, B or C
read_all_results <- function(type, d_, small = F) {
  c_all <-  c(0,0.1,0.2,0.3,0.6,0.9,0.99)
  res <- data.frame()
  for (c_ in c_all) {
    res_c <- read_results(type = type, d_ = d_, c_ = c_, small = small)
    res_b <- read_results(type = type, d_ = d_, c_ = c_, small = small, bs = T)
    res_c$bootstrap <- "no"
    res_b$bootstrap <- "yes"
    rr <- 1:nrow(res_b)
    res_i <- rbind(res_c[rr,], res_b[rr,])
    res_i$c <- c_
    res <- rbind(res, res_i)
  }
  res <- if (small) {res[res$n <= 50,]} else {res} #res[res$n <= 1500,]
  return(tibble::as_tibble(res))
}


# Create tibble with results for all c, for chosen type A, B or C
read_all_without_bootstrap <- function(type, d_, small = F) {
  c_all <-  c(seq(0,0.9,0.1), 0.99)
  res <- data.frame()
  for (c_ in c_all) {
    res_i   <- read_results(type = type, d_ = d_, c_ = c_, small = small, bs = F)
    res_i$c <- c_
    res <- rbind(res, res_i)
  }
  res <- if (small) {res[res$n <= 50,]} else {res}
  return(tibble::as_tibble(res))
}


# Plot all grid of coverage and widths for all correlations except 0, with
# either large or small sample sizes
plot_cw_grid <- function(type, d_, small = F, c0 = F) {
  res <- read_all_results(type = type, d_ = d_, small = small)
  
  keep <- if (c0) {res$c == 0} else {res$c != 0}
  
  p1 <- ggplot(res[keep,], aes(n, coverage)) +
    geom_point(aes(colour = bootstrap, shape = bootstrap)) +
    geom_line(aes(colour = bootstrap, linetype = bootstrap)) +
    scale_color_manual(values = c("darkred", "steelblue")) + 
    geom_ribbon(aes(ymin = covrgCI_L, ymax = covrgCI_U,
                    fill = bootstrap, colour = bootstrap, 
                    linetype = bootstrap), alpha = 0.13) +
    facet_wrap(~c, scales = "free", ncol = 1) + 
    theme_set(theme_minimal()) +
    theme(legend.position = "none") +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x))
  
  p2 <- ggplot(res[keep,], aes(n, width_mu)) +
    geom_point(aes(colour = bootstrap, shape = bootstrap)) +
    geom_line(aes(colour = bootstrap, linetype = bootstrap)) +
    scale_color_manual(values = c("darkred", "steelblue")) + 
    geom_ribbon(aes(ymin = widthCI_L, ymax = widthCI_U,
                    fill = bootstrap, colour = bootstrap, 
                    linetype = bootstrap), alpha = 0.13) +
    facet_wrap(~c, scales = "free", ncol = 1) +
    ylab("mean width") + 
    theme_set(theme_minimal()) +
    #theme(legend.position = "none") +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x))
  
  size <- if (small) {"small"} else {"large"}
  title <- paste0("Coverage and Width Comparisons '",type,"' (",size," Samples)")
  gridExtra::grid.arrange(
    p1, p2, ncol = 2) #, top = grid::textGrob(tools::toTitleCase(title)))
}

make_plot_pdf2 <- function(type, small, c0 = F) {
  size <- if(small) {"small"} else {"large"}
  size <- if(c0) {paste0(size,"_c0")} else {size}
  pdf(paste0("MC_figures/",type,"_",size,".pdf"))
    plot_cw_grid(type,"22122019", small = small, c0 = c0)
  dev.off()
}

make_all_plots2 <- function(types = c("A","B","C"), sizes = c(T, F)) {
  for (type in types) {
    for (small in sizes) {
      make_plot_pdf2(type = type, small = small)
    }
  }
}


## Coverage vs c plots
cov_vs_c <- function(dat, ns = seq(15,50,5), main, legpos = "bottomleft") {
  require(ggplot2)
  require(reshape2)
  c1 <- data.frame(c =  c(seq(0,0.9,0.1), 0.99))
  for (n in ns) {
    c1 <- cbind(c1, dat[dat$n == n, c("coverage")])
  }
  names(c1) <- c("c", as.character(ns))
  c1 <- as_tibble(c1)
  c1
  matplot(c1[,1], c1[,-c(1)], type = c("b"), 
          pch = 1, col = rev(1:length(ns)), lty = 1:length(ns),
          xlab = "c", ylab = "estimated coverage")
  legend(legpos, legend = ns, title = "n",
         pch = 1, col = rev(1:length(ns)), lty = 1:length(ns)) # optional legend
}


## Coverage vs c plots
width_vs_c <- function(dat, ns = seq(15,50,5), main, legpos = "bottom") {
  require(ggplot2)
  require(reshape2)
  c1 <- data.frame(c =  c(seq(0,0.9,0.1), 0.99))
  for (n in ns) {
    c1 <- cbind(c1, dat[dat$n == n, c("width_mu")])
  }
  names(c1) <- c("c", as.character(ns))
  c1 <- as_tibble(c1)
  c1
  matplot(c1[,1], c1[,-c(1)], type = c("b"), 
          pch = 1, col = 1:length(ns), lty = 1:length(ns),
          xlab = "c", ylab = "mean width")
  legend(legpos, legend = ns, title = "n",
         pch = 1, col = 1:length(ns), lty = 1:length(ns)) # optional legend
}

