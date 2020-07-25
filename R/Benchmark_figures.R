##### The data here was entered manually from benchmark results seen in Julia
##### on Daniel Fryer's computer - intel i7-9700k processor with 16GB DDR4 RAM at 3000MHz.


library(tidyverse)

d <- c(3, 3, 3, 6)
n <- c(1000,5000,10000,1000)

# Bootstrap
memr_b <- c(88.21, 356.94, 692.9, 307.6)                 # GiB memory
allc_b <- c(425647829, 439476845, 455514240, 3430806647) # allocations
mean_b <- c(46.72, 154.085, 289.275, 256.309)            # seconds

# Asymptotic
memr_s <- c(2.1, 5.03, 8.7, 182.7)                       # GiB memory
allc_s <- c(26670093, 59502765, 100542815, 1622018399)   # allocations
mean_s <- c(1.711, 3.293, 5.448, 140.562)                # seconds


# Ratios of memory and mean time
df <- data.frame(n = n, memr = memr_b/memr_s, mean = mean_b/mean_s)[1:3,]
ggplot(df, aes(x = reorder(as.character(n),n), y = mean)) +
  geom_col()

df <- data.frame(n = n, mem_boot = memr_b, mem_shap = memr_s,
                 mean_shap = mean_s, mean_boot = mean_b,
                 order = c(1,2,3,4))[1:3,]



g1 <- ggplot(df, aes(x = reorder(as.character(n),n), y = mem_boot/mem_shap)) +
  geom_col() +
  ylab("Bootstrap memory use / Normal memory use") +
  xlab("Sample size (n)") +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        title = element_text(size = rel(2))) +
  ggtitle("Memory usage ratios") +
  geom_text(size = 7, aes(label = round(mem_boot/mem_shap, 3)),
            position = position_stack(vjust = 0.5))

g2 <- ggplot(df, aes(x = reorder(as.character(n),n), y = mean_boot/mean_shap)) +
  geom_col() +
  ylab("Bootstrap time / Normal time") +
  xlab("Sample size (n)") +
  ggtitle("Mean execution time ratios") +
  theme_bw() +
  theme(axis.text = element_text(size = rel(2)),
        title = element_text(size = rel(2))) +
  geom_text(size = 7, aes(label = round(mean_boot/mean_shap, 3)),
            position = position_stack(vjust = 0.5))
#X11()
#gridExtra::grid.arrange(g1, g2, ncol = 2, 
#      top = grid::textGrob("Benchmark metric ratios (Bootstrap / Normal)", 
#      gp=grid::gpar(fontsize=25,font=1)))


par()

pdf(file="MC_figures/benchmark_ratios.pdf",width=15,height=7.5)
gridExtra::grid.arrange(g1, g2, ncol = 2)
dev.off()


# memory vs sample size shap and boot together
df <- tidyr::gather(
  df, mem_boot, mem_shap, key = "type", value = "mem") %>% 
  dplyr::arrange(n)
ggplot(df, aes(x = reorder(paste0(type," (",n,")"),order), y = mem)) +
  geom_col() +
  xlab("sample size")
