source("Real_estate_helpers.R")
library(rgdal)
library(rgeos)
library(maptools)
library(broom)

path <- "../real_estate_dataset/"
data <- readRDS(paste0(path,"scrape_all_2019_to_2020_raw_20200714.Rds"))
targets <- read.csv(paste0(path,"suburbs_list_20200714.csv"))[,2:3]

trans_path <- paste0(path, "trans/")
data_clean <- clean_data(data, date = "20200718", name = "FULL", save = F, path = path)
data_trans <- transform_data(data_clean, date = "20200718", name = "FULL", save = F, path = trans_path)

nrow(data_trans) # n = 13,291

# VISUALISING WITH MAP DATA  --------------------------------------
data_clean <- readRDS(paste0(path, "data_clean_FULL_20200718.Rds"))
data_trans <- readRDS(paste0(path, "trans/data_trans_FULL_20200718.Rds"))

dists <- data_clean$CBD
data_clean$dist_group <- factor(ceiling(dists/10))
X11(); boxplot(land ~ dist_group, data = data_clean)
hist(data_clean$land)

sample_ind <- sample(1:nrow(data_clean), 1000)
X11(); ggplot(data_clean[sample_ind,]) +
  geom_point(aes(x = CBD, y = price, size = land, colour = room))
1203/953
1824/1130

hist(data_clean$CBD)

# looks like everyone rushed to buy houses at start of covid.
barplot(table(data_trans$date))
# Zipf-looking suburb barplot
barplot(sort(table(data_trans$suburb)))

top_5_suburbs <- rev(sort(table(data_trans$suburb)))[1:5]

# VIC Suburb/Locality Boundaries - PSMA Administrative Boundaries
# -- Source: https://data.gov.au/dataset/ds-dga-af33dd8c-0534-4e18-9245-fc64440f742e/details
path2 <- "../../../real_estate_dataset/"
spdf <- readOGR(dsn= paste0(path2,"/vic_locality_polygon_shp"))
mapdata <- tidy_map_data(spdf, data_trans, targets)

# Based on the below plot, we dropped hillside and bellfield earlier

X11(); plot_map_with_suburb_labels(mapdata)
# Note that CBD is bimodal
histograms(data_trans)

# We remove this suburb cocoroc because it does not plot well (and only has one house)
data_clean <- data_clean[data_clean$suburb != "COCOROC",]
# Now compare close to far
mapdata_close <- tidy_map_data(spdf, data_trans[data_clean$CBD <= 25,], targets)
mapdata_far <- tidy_map_data(spdf, data_trans[data_clean$CBD > 25,], targets)
pdf(file="suburb_map.pdf",width=6,height=5)
plot_map_with_near_and_far(mapdata_close, mapdata_far)
dev.off()

# Now compare 2019 to 2020
mapdata2019 <- tidy_map_data(spdf, data_trans[data_trans$year == "2019",], targets)
mapdata2020 <- tidy_map_data(spdf, data_trans[data_trans$year == "2020",], targets)
X11(); plot_map(mapdata2019)
X11(); plot_map(mapdata2020)

###################################################################
# SHAPLEY VALUES ARE CALCULATED IN Julia/Real_estate_application.jl
###################################################################
date_condition <- function(date) {
  (date >= "2019-02-01" & date <= "2019-04-01") |
  (date >= "2020-02-01" & date <= "2020-04-01")
}
date_table_df <- date_counts(data_trans$date)
#pdf(file="dates_barplot.pdf",width=5,height=4)
plot_date_table(date_table_df)
#dev.off()

data_st <- split_and_transform(data_clean, date_condition = date_condition, dist = 25)
sum(sapply(data_st, nrow))
sapply(data_st, nrow)

m1 <- lm(price ~ CBD+images+land+school+station+room, data = data_st$near_2019)
m2 <- lm(price ~ CBD+images+land+school+station+room, data = data_st$far_2019)
m3 <- lm(price ~ CBD+images+land+school+station+room, data = data_st$near_2020)
m4 <- lm(price ~ CBD+images+land+school+station+room, data = data_st$far_2020)
summary(m1)$r.sq
summary(m2)$r.sq
summary(m3)$r.sq
summary(m4)$r.sq

#plot(m1)

pfx <- "data_shapley_"
shp_fourway <- read.csv(paste0(path, pfx, "fourway_split_20200718.csv"))
shp_year <- read.csv(paste0(path, pfx, "year_split_20200718.csv"))

# The data for this was transformed in four groups: by year/dist

shp_fourway$year <- factor(shp_fourway$year)
shp_fourway$type <- factor(shp_fourway$type) # Type is distance near/far
View(shp_fourway)
pdf(file="realestate_shapley2.pdf",width=5,height=4)
ggplot(data = shp_fourway) +
  geom_point(aes(x = feature, y = Shapley, shape = year)) +
  geom_crossbar(
    fatten=1, alpha=0.3, width=0.3, linetype=0,
    aes(y=Shapley, ymin = CIL, ymax = CIU, x = feature,
        colour = year, fill = year)) +
  facet_grid(type ~ .)
dev.off()

# The data for this was transformed in four groups: by year/dist
pdf(file="realestate_shapley.pdf",width=5,height=4)
ggplot(data = shp_fourway) +
  geom_point(aes(x = feature, y = Shapley, shape = type)) +
  geom_crossbar(
    fatten=1, alpha=0.3, width=0.3, linetype=0,
    aes(y=Shapley, ymin = CIL, ymax = CIU, x = feature,
        colour = type, fill = type)) +
  facet_grid(year ~ .)
dev.off()

# The data for this was transformed in two groups: by year
shp_year$year <- factor(shp_year$year)
ggplot(data = shp_year) +
  geom_point(aes(x = feature, y = Shapley)) +
  geom_crossbar(
    fatten=1, alpha=0.3, width=0.3, linetype=0,
    aes(y=Shapley, ymin = CIL, ymax = CIU, x = feature,
        colour = year, fill = year))
