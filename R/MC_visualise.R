source("MC_visualise_helpers.R")

################################################
## BOOTSTRAP COVERAGE AND WIDTH COMPARISON PLOTS
################################################

# Look at all the plots one at a time
X11()
plot_cw_grid("A", "22122019")
plot_cw_grid("B", "22122019")
plot_cw_grid("C", "22122019")
plot_cw_grid("A", "22122019", small = T)
plot_cw_grid("B", "22122019", small = T)
plot_cw_grid("C", "22122019", small = T)

# These are the plots for correlation 0, treated separately
plot_cw_grid("A", "22122019", c0 = T)
plot_cw_grid("B", "22122019", c0 = T)
X11(); plot_cw_grid("C", "22122019", c0 = T)
plot_cw_grid("A", "22122019", small = T, c0 = T)
plot_cw_grid("B", "22122019", small = T, c0 = T)
plot_cw_grid("C", "22122019", small = T, c0 = T)

# Make all the plot pdfs, except for those with c = 0.
make_all_plots2()

# Make some plot pdfs with c = 0
type <- "C"; size <- "large"; small = F
pdf(paste0("MC_figures/",type,"_",size,"_c0.pdf"), height= 3)
  plot_cw_grid(type,"22122019", small = small, c0 = T)
dev.off()
type <- "C"; size <- "small"; small = T
pdf(paste0("MC_figures/",type,"_",size,"_c0.pdf"), height= 3)
  plot_cw_grid(type,"22122019", small = small, c0 = T)
dev.off()
type <- "B"; size <- "small"; small = T
pdf(paste0("MC_figures/",type,"_",size,"_c0.pdf"), height= 3)
  plot_cw_grid(type,"22122019", small = small, c0 = T)
dev.off()


####################################################
## GENERAL STATEMENTS, COV vs c and WIDTH vs c PLOTS
####################################################

# Read and bind data data
Al  <- read_all_without_bootstrap("A", "22122019", small = F)
As  <- read_all_without_bootstrap("A", "22122019", small = T)
Bl  <- read_all_without_bootstrap("B", "22122019", small = F)
Bs  <- read_all_without_bootstrap("B", "22122019", small = T)
Cl  <- read_all_without_bootstrap("C", "22122019", small = F)
Cs  <- read_all_without_bootstrap("C", "22122019", small = T)
A <- rbind(Al, As)
B <- rbind(Bl, Bs)
C <- rbind(Cl, Cs)

## Width vs c plot
pdf(file="MC_figures/width_vs_cor.pdf",width=10,height=5.8)
par(mfrow = c(1,2), cex = 1.2, mar=c(4,4,1,1)+0.1)
width_vs_c(A)
width_vs_c(A, ns = seq(1000,2000,200))
dev.off()

## Cov vs c plots
pdf(file="MC_figures/cov_vs_cor.pdf",width=10,height=5.5)
par(mfrow = c(1,1), mar=c(5,3,2,2)+0.1)
cov_vs_c(A)
dev.off()
cov_vs_c(B, ns = seq(100,2000,300), legpos = "bottomright")

### Observations on lower bound of confidence interval for coverage

#For $c \geq 0.3$ and $n > 100$ the lower bound of the confidence interval for the coverage never drops below $0.9$.
A[A$covrgCI_L  < 0.9 & A$c >= 0.3 & A$n > 100,]
B[B$covrgCI_L  < 0.9 & B$c >= 0.3 & B$n > 100,]
C[C$covrgCI_L  < 0.9 & C$c >= 0.3 & C$n > 100,]

# In all three studies the estimated coverage probability is greater than 0.85 for all sample sizes greater than 10 and all c
A[A$coverage  < 0.85 & A$n >= 10,]
B[B$coverage  < 0.85 & B$n >= 10,]
C[C$coverage  < 0.85 & C$n >= 10,]

# For $c \geq 0.2$ and $n > 15$ the lower bound of the confidence interval for coverage never drops below $0.86$. 
A[A$covrgCI_L < 0.86 & A$c >= 0.2 & A$n > 15,]
B[B$covrgCI_L < 0.86 & B$c >= 0.2 & B$n > 15,]
C[C$covrgCI_L < 0.86 & C$c >= 0.2 & C$n > 15,]

# For $c = 0.1$ and $n \in [10,100]$, the lower bound of the confidence interval for coverage never drops below $0.91$.
A[A$covrgCI_L < 0.91 & A$c <= 0.1 & A$n >= 10 & A$n <= 100,]
B[B$covrgCI_L < 0.91 & B$c <= 0.1 & B$n >= 10 & B$n <= 100,]
C[C$covrgCI_L < 0.91 & C$c <= 0.1 & C$n >= 10 & C$n <= 100,]

# For $c = 0$ and $n \geq 15$ the lower bound of the confidence interval for coverage never drops below $0.95$ in Studies A and B, while in Study C the lower bound is $0.88$ at the lowest.
A[A$covrgCI_L < 0.95  & A$c == 0 & A$n >= 15,]
B[B$covrgCI_L < 0.95  & B$c == 0 & B$n >= 15,]
C[C$covrgCI_L < 0.95  & C$c == 0 & C$n >= 15,]

####################################################
## GENERAL STATEMENTS COMPARING TO BOOTSTRAP
####################################################
Al  <- read_all_results("A", "22122019", small = F)
As  <- read_all_results("A", "22122019", small = T)
Bl  <- read_all_results("B", "22122019", small = F)
Bs  <- read_all_results("B", "22122019", small = T)
Cl  <- read_all_results("C", "22122019", small = F)
Cs  <- read_all_results("C", "22122019", small = T)
A <- rbind(Al, As)
B <- rbind(Bl, Bs)
C <- rbind(Cl, Cs)

A[A$bootstrap != "yes",]$width_mu <= A[A$bootstrap == "yes",]$width_mu + 0.001

A[A$bootstrap == "yes",][137:nrow(A[A$bootstrap == "yes",]),]

A_both <- dplyr::inner_join(A[A$bootstrap != "yes",], A[A$bootstrap == "yes",], by = c("n","c"))

Awb <- A_both[A_both$n < 50, c("width_mu.x", "width_mu.y")]

A_both[A_both$n < 50,c("width_mu.x", "width_mu.y","n","c")]

A_both <- dplyr::mutate(A_both, w = width_mu.y - width_mu.x < 0)

A_both[A_both$n < 50  & A_both$w == T,c("width_mu.x", "width_mu.y","n","c","w")]

ggplot2(Awb, aes(width_mu.y-width_mu.x))
plot(Awb$width_mu.y-Awb$width_mu.x, pch = A_both$n)
abline(a = 0, b = 0)

matplot(Awb, type = c("b"), 
        pch = 1, col = 1:2, lty = 1:2,
        xlab = "c", ylab = "mean width")


