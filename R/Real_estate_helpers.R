library(tidyverse)
library(realscrape)
library(caret)


clean_data <- function(data, date = "TEST", name = "TEST", save = T,
                       path = "") {
  
  date_condition <- function(date) {
    (date >= "2019-01-01" & date <= "2019-07-18") |
    (date >= "2020-01-01" & date <= "2020-07-18")
  }
  
  data_clean <- data %>% 
    dplyr::select( address, price, type, bed, date, 
                   bath, car, CBD, images, land, 
                   school, station, suburb, postcode ) %>% 
    tidyr::drop_na() %>%
    dplyr::distinct() %>% 
    dplyr::mutate(room = bed + bath + car) %>% 
    dplyr::filter(land != -1 &
                  date_condition(date) &
                  address != "DUD PAGE" &
                  type == "house" &
                  !(suburb %in% c("BELLFIELD","HILLSIDE"))) %>% # distant suburbs
    dplyr::select(-type) 
  data_clean <- data_clean %>% 
    dplyr::filter(land < 9e5 & 
                    land > 10 &
                    price > 1e5 &
                    bed >= bath/2 &
                    bath < 6 &
                    bed < 8 &
                    CBD <= quantile(data_clean$CBD, probs = 0.98) &
                    land <= quantile(data_clean$land, probs = 0.97) &
                    land >= quantile(data_clean$land, probs = 0.03))
  data_clean$year <- substr(data_clean$date,1,4)
  if (save) {
    saveRDS(data_clean, paste0(path,"data_clean_", name,"_",date,".Rds"))
    write.csv(data_clean, paste0(path,"data_clean_", name,"_",date,".csv"))
  }
  return(data_clean)
}

transform_data <- function(data, date = "TEST", name = "TEST", save = T,
                           path  = "") {
  # Perform YeoJohnson, center and scaling transformations
  data_cont <- data %>% 
    select(price, CBD, images, land, school, station, room)
  pp <- preProcess(data_cont, method = c("center", "scale", "YeoJohnson"))
  data_trans <- predict(pp, data_cont)
  # Add suburb column
  data_trans$suburb <- data$suburb
  data_trans$date <- data$date
  data_trans$year <- substr(data_trans$date,1,4)
  # Save transformed data
  if (save) {
    saveRDS(data_trans, paste0(path,"data_trans_",name,"_",date,".Rds"))
    write.csv(data_trans, paste0(path,"data_trans_",name,"_",date,".csv"))
  }
  return(data_trans)
}

histograms <- function(data_trans) {
  hist(data_trans$images)
  hist(data_trans$price)
  hist(data_trans$land)
  hist(data_trans$CBD)
  hist(data_trans$school)
  hist(data_trans$station)
  hist(data_trans$room)
}


plot_map_with_suburb_labels <- function(mapdata) {
  centroids <- mapdata$centroids
  polygons <- mapdata$polygons
  ggplot() +
    scale_fill_gradient(low="gray", high="black") +
    geom_polygon(data = polygons, 
                 aes( x = long, y = lat, group = group, fill = n), 
                 color="white") +
    theme_void() +
    geom_point(data = centroids[centroids$SUBURB %in% polygons$SUBURB,], 
               aes(x = centroid.x, y = centroid.y)) +
    geom_text(data = centroids[centroids$SUBURB %in% polygons$SUBURB,],
               aes(x = centroid.x, y = centroid.y, label = SUBURB),
              check_overlap = TRUE, fontface = "bold")
}

plot_map <- function(mapdata) {
  centroids <- mapdata$centroids
  polygons <- mapdata$polygons
  ggplot() +
    scale_fill_gradient(low="gray", high="black") +
    geom_polygon(data = polygons, 
                 aes( x = long, y = lat, group = group, fill = n), 
                 color="white") +
    theme_void() +
    geom_point(data = centroids[centroids$SUBURB %in% polygons$SUBURB,], 
               aes(x = centroid.x, y = centroid.y))
}


### WARNING: RADIUS DOES NOT APPEAR TO BE CORRECT FOR THIS CIRCLE POLYGON
plot_map_with_circle <- function(
  mapdata, radius = 10, center = c(144.9631, -37.8136)) {
  require(swfscMisc)
  circle <- data.frame(circle.polygon(center[1], center[2], radius, units = "km"))
  centroids <- mapdata$centroids
  polygons <- mapdata$polygons
  ggplot() +
    scale_fill_gradient(low="gray", high="black") +
    geom_polygon(data = polygons, 
                 aes( x = long, y = lat, group = group, fill = n), 
                 color="white") +
    theme_void() +
    geom_point(data = centroids[centroids$SUBURB %in% polygons$SUBURB,], 
               aes(x = centroid.x, y = centroid.y)) +
    geom_path(data = circle, aes(x = x, y = y), size = 1.5, colour = "darkred",
              linetype = "dashed")  +
    geom_point(aes(x = center[1], y = center[2]), 
               colour = "darkred", size = 4, shape = 18)
}

plot_map_with_near_and_far <- function(mapdata_close, mapdata_far,
                                       center = c(144.9631, -37.8136)) {
  
  ggplot() +
    scale_fill_gradient(low="gray", high="black") +
    geom_polygon(data = mapdata_far$polygons, 
                 aes( x = long, y = lat, group = group, fill = n), 
                 color="darkcyan", linetype = "blank") +
    geom_polygon(data = mapdata_close$polygons, 
                 aes( x = long, y = lat, group = group, fill = n), 
                 color="violetred3", linetype = "blank") +
    geom_path(data = mapdata_close$polygons, 
              aes( x = long, y = lat, group = group),
              color="violetred3", alpha = 1, linetype = "dashed",
              size = 0.2) +
    geom_path(data = mapdata_far$polygons, 
              aes( x = long, y = lat, group = group),
              color="darkcyan", alpha = 1, 
              size = 0.2) +
    theme_void() +
    geom_point(aes(x = center[1], y = center[2]), 
               colour = "darkred", size = 2, shape = 18) +
    xlab("longitude") +
    ylab("latitude")
}


make_circle_line <- function(x,y,r, n = 100) {
  pts <- seq(0, 2 * pi, length.out = n)
  require(sp)
  xy <- cbind(x + r * sin(pts), y + r * cos(pts))
  sl <- SpatialLines(list(Lines(list(Line(xy)), "line")))
  return(sl)
}



tidy_map_data <- function(spdf, data, targets) {
  spdf_df <- tidy(spdf, region = "VIC_LOCA_2")
  suburb_table <- table(data$suburb)
  names(dimnames(suburb_table)) <- "SUBURB"
  suburb_table_df <- tidy(suburb_table)
  targets$suburb <- as.character(targets$suburb)
  polygons <- inner_join(suburb_table_df, spdf_df, by = c("SUBURB" = "id")) %>% 
    inner_join(targets, by = c("SUBURB" = "suburb"))
  centroid_coords <- gCentroid(spdf, byid=TRUE)
  centroids <- as_tibble(data.frame(SUBURB = as.character(spdf@data$VIC_LOCA_2), 
                          centroid = centroid_coords@coords,
                          stringsAsFactors = F))
  return(list(polygons = polygons, centroids = centroids))
}




# We need to split the data into groups and transform the independently
split_and_transform <- function(data_clean, date_condition, dist) {
  data_clean <- data_clean[date_condition(data_clean$date),]
  data_clean_near <- data_clean[data_clean$CBD <= dist,]
  data_clean_far <- data_clean[data_clean$CBD > dist,]
  data_clean_near_2019 <- data_clean_near[data_clean_near$year == "2019",]
  data_clean_far_2019 <- data_clean_far[data_clean_far$year == "2019",]
  data_clean_near_2020 <- data_clean_near[data_clean_near$year == "2020",]
  data_clean_far_2020 <- data_clean_far[data_clean_far$year == "2020",]
  data_clean_2019 <- data_clean[data_clean$year == "2019",]
  data_clean_2020 <- data_clean[data_clean$year == "2020",]
  data_trans_2019      <- transform_data(data_clean_2019, date = "20200718", name = "2019", save = T, path = trans_path)
  data_trans_2020      <- transform_data(data_clean_2020, date = "20200718", name = "2020", save = T, path = trans_path)
  data_trans_near_2019 <- transform_data(data_clean_near_2019, date = "20200718", name = "near_2019", save = T, path = trans_path)
  data_trans_far_2019  <- transform_data(data_clean_far_2019 , date = "20200718", name = "far_2019" , save = T, path = trans_path)
  data_trans_near_2020 <- transform_data(data_clean_near_2020, date = "20200718", name = "near_2020", save = T, path = trans_path)  
  data_trans_far_2020  <- transform_data(data_clean_far_2020 , date = "20200718", name = "far_2020" , save = T, path = trans_path)
  cat("near 2019, n = ", nrow(data_trans_near_2019), "\n")
  cat("far 2019, n = ",  nrow(data_trans_far_2019), "\n")
  cat("near 2020, n = ", nrow(data_trans_near_2020), "\n")
  cat("far 2020, n = ",  nrow(data_trans_far_2020), "\n")
  return(list(near_2019 = data_trans_near_2019, 
              far_2019  = data_trans_far_2019,
              near_2020 = data_trans_near_2020, 
              far_2020  = data_trans_far_2020))
}

date_counts <- function(dates) {
  all_dates <- data.frame( date = 
                             c(seq(as.Date("2019-01-01"), as.Date("2019-07-18"), by="days"),
                               seq(as.Date("2020-01-01"), as.Date("2020-07-18"), by="days")))
  all_dates$date <- as.character(all_dates$date)
  date_table <- table(dates[date_condition(dates)])
  names(dimnames(date_table)) <- "date"
  date_table_df <- tidy(date_table) %>% 
    right_join(all_dates, by = "date") %>% 
    mutate(year = substr(date, 1,4), date = substr(date, 6, 10))
  date_table_df$n[is.na(date_table_df$n)] <- 0 
  names(date_table_df) <- c("date", "sales", "year")
  return(date_table_df)
}

plot_date_table <- function(date_table_df) {
  keep_dates <- c("01-01", "02-01", "03-01", "04-01", "05-01", "06-01", "07-01", "07-18")
  dl <- unique(date_table_df$date)
  keep <- dl %in% keep_dates
  date_labels <- dl[keep]
  ggplot(data = date_table_df, aes(x = date, y = sales)) + 
    geom_bar(stat = "identity", width = 0.925, colour = "turquoise4") +
    facet_grid(year~.) +
    geom_vline(aes(xintercept = c("02-01")), 
               colour = "darkred", size = 0.5, linetype = "dashed") +
    geom_vline(aes(xintercept = c("04-01")), 
               colour = "darkred", size = 0.5, linetype = "dashed") +
    scale_x_discrete(breaks = date_labels) +
    theme_light()
}
