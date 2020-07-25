# Shapley value confidence intervals for attributing variance explained (via R squared)
All data (in .csv and .rds format), figures (in .pdf format) and code (in the R and Julia languages), for the paper [Shapley value confidence intervals for attributing variance explained](https://arxiv.org/abs/2001.09593) by [Daniel Fryer](https://danielvfryer.com), [Inga Strumke](https://strumke.com) and [Hien Nguyen](https://hiendn.github.io/).

The [Monte Carlo](Julia/MC_produce_22122019.jl) simulation [results](Julia/results) and [Shapley values](Julia/Shapley.jl) were calculated in [Julia](Julia), because this worked out to be much faster than doing it in R, even when we tried with [parallel processing](Shapley_noC.R), and even when we [used C via Rcpp](R/Shapley_C.R).

The [figures](R/Figures) (including Monte Carlo and benchmark figures), and some [exploratory analysis](R/Real_estate_applications.R) of the real estate data, were produced in [R](R). Shapley values for the real estate data were [produced in Julia](Julia/Real_estate_application.jl).

Enjoy! Feel free to contact me if you have any questions or concerns.

