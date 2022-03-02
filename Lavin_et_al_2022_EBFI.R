# Data and analyses for:
# Warmer temperature decreases the maximum length
# of six species of marine fishes, crustacean and squid in New Zealand 

# Lavin CP, Gordo-Vilaseca C, Stephenson F, Shi Z, Costello MJ
# Environmental Biology of Fishes
# 2022

# Corresponding author: Charles P. Lavin, email: charles.p.lavin@nord.no

# Code for source functions "AssessRelativeImportancePredictors.R",
# "CustomizedBarPlot.R" and "ModelEvaluation.R" are available upon request.


#loading packages
library( tidyverse )
library( dplyr )
library( ggplot2 )
library( ggpubr )
library( ggpmisc )
library(magrittr)
library( grid )
library( gridExtra )
library( fitdistrplus )
library( Hmisc )
library( corrplot )
library( mgcv )
library( gratia )
library( raster )
library( sdmpredictors )
library( ggmap )
library( sdm )
library( spbabel )
library( maptools )
library( SoDA )
library(grid)
library(ggplotify)

####################################################

# Snapper
length_data = read.csv( "SNA_trawl_lengths_DATES.csv" )
length_data = length_data %>%
  mutate(lgth_SL = ((lgth-0.412)/1.152)) # FORK LENGTH to STANDARD LENGTH
# Conversion from Ferrell & Sumpton 1997

dim( length_data )
nrow( length_data ) #### 73572 lines
names( length_data )
table( length_data$trip_code ) 
summary( length_data$station_no ) #### 1 NA
summary( length_data$dlon_s ) #### 171 NA's
summary( length_data$dlat_s ) #### 171 NA's
summary( length_data$dlon_e ) #### 9651 NA's
summary( length_data$dlat_s ) #### 171 NA's
summary( length_data$min_gdepth ) #### 9234 NA's
summary( length_data$max_gdepth )#### 9206 NA's
summary( length_data$bot_temp ) #### 49427 NA's
table( length_data$species ) 
summary( length_data$lgth_SL ) #### 1 NA
summary( length_data$no_a ) #### 1 NA
summary( length_data$no_m )#### 325 NA's
summary( length_data$no_f ) #### 325 NA's

length_data_1 = drop_na( length_data )
nrow( length_data_1 ) #### 24077 lines remaining
length_data_2 = length_data[!is.na( length_data$lgth_SL ) & !is.na( length_data$dlon_s ) &
                              !is.na( length_data$dlat_s ) & !is.na( length_data$bot_temp ) &
                              !is.na( length_data$min_gdepth ) & !is.na( length_data$max_gdepth ),]
nrow( length_data_2 ) #### 24113 lines remaining
length_data <- length_data_2 

plot( length_data$dlon_s, length_data$dlat_s )

#Trip stations
summary( length_data$station_no ) 
length( unique( length_data$station_no ) )
table( length_data$station_no ) 
trip_station <- paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) )
length( unique( trip_station ) )
length_data$trip_station <- as.factor( paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) ) )
dim( length_data )
names( length_data )

#Mean depth for each station
length_data = length_data %>%
  add_column(mean_depth = ((length_data$min_gdepth + length_data$max_gdepth)/2))

#Mean length data for each station
mean_length_data <- aggregate( lgth_SL ~ trip_station, data = length_data, mean )
mean_length_data$dlon_s <- length_data$dlon_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$dlat_s <- length_data$dlat_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$bot_temp <- length_data$bot_temp[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$mean_depth <- length_data$mean_depth[match( mean_length_data$trip_station, length_data$trip_station)]
mean_length_data$date_s <- length_data$date_s[match( mean_length_data$trip_station,length_data$trip_station)]
dim( mean_length_data )
names( mean_length_data )
thedata = mean_length_data
plot( thedata$dlon_s, thedata$dlat_s )

#adding 95th quantile value of length for each Trip station
length_data = as_tibble(length_data)

quant_length_data = length_data %>%
  dplyr::select(date_s, trip_station, lgth_SL, dlon_s, dlat_s, bot_temp)

quant_data <-  data.frame(quant_length_data$date_s, quant_length_data$trip_station, quant_length_data$bot_temp, quant_length_data$dlon_s, quant_length_data$dlat_s, quant_length_data$lgth_SL )
names( quant_data ) <- c("Date", "Trip.Station",  "Situ.Temp", "Lon", "Lat", "Length" ) 

bXY = geoXY( latitude = quant_data$Lat, longitude = quant_data$Lon, 
             lat0 = mean( quant_data$Lat ), lon0 = mean( quant_data$Lon ), unit = 1000 )
quant_data = cbind( quant_data, bXY )

sort( unique( quant_data$Trip.Station ) )

#set stations as factor
quant_data$Trip.Station <- as.factor( paste0( as.factor( quant_data$Trip.Station )))

#here are the 95 quantile values for each trip station
quant_dat = quant_data %>%
  group_by(Trip.Station) %>%
  summarise(enframe(quantile(Length, c(0.95)), "quantile", "95.Length"))

quant_95=as_tibble(quant_dat$`95.Length`)


######## Load in environmental data (mean Temperature.bio, and mean oxygen at maximum depth)
######## Those data were extracted from https://www.bio-oracle.org/downloads-to-email.php

ox = read.asciigrid( "oxygen.asc" )
ox = raster( ox )
x1 <- crop( ox, extent( -180, 0, -90, 90 ) )
x2 <- crop( ox, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
ox <- merge( x1, x2 )
plot( ox, xlab="Longitude", ylab="Latitude" )
temp = read.asciigrid( "temp.asc" )
temp = raster( temp )
x1 <- crop( temp, extent( -180, 0, -90, 90 ) )
x2 <- crop( temp, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
temp <- merge( x1, x2 )
plot( temp, xlab="Longitude", ylab="Latitude" )

predictors <- stack( temp, ox )

######## Combine environmental data with monitoring data and check for outliers in the response and covariate variables
#### Combine environmental data with monitoring data
data_new <-  data.frame( raster::extract( predictors, thedata[3:4] ),
                         thedata$date_s, thedata$bot_temp, thedata$mean_depth, thedata$dlon_s, thedata$dlat_s, thedata$lgth_SL )

#adding 95 quantiles
data_new = cbind(data_new, quant_95)
names( data_new ) <- c( "Bio.Temp", "Oxygen", "Date", "Situ.Temp", "Mean.depth", "Lon", "Lat", "Length", "Max.length" ) 
nrow( data_new ) ## 1688 lines 


#### Check for outliers in the response variable
summary( data_new$Max.length )
hist( data_new$Max.length )
boxplot( data_new$Max.length )
sort( unique( data_new$Max.length ) )

#### Check for outliers in the bottom bio.temp
summary( data_new$Bio.Temp )
hist( data_new$Bio.Temp )
boxplot( data_new$Bio.Temp )
sort( unique( data_new$Bio.Temp ) )

#### Check for outliers in the bottom dissolved oxygen (DO) concentrations
summary( data_new$Oxygen ) ## 6 NA's
data_new <- data_new[!is.na( data_new$Oxygen ),] 
nrow( data_new ) ## 1555 lines remaining
summary( data_new$Oxygen )
hist( data_new$Oxygen )
boxplot( data_new$Oxygen )
sort( unique( data_new$Oxygen ) )

#### Check for outliers in the mean depth
summary( data_new$Mean.depth )
hist( data_new$Mean.depth )
boxplot( data_new$Mean.depth )
sort( unique( data_new$Mean.depth ) )

### The diagnostics suggest that minimum DO concentrations,
### maximum mean depth, three minimum temp, are outliers and
### therefore, they should be removed

#DO
nrow( data_new ) #1682
data_new <- data_new[!data_new$Oxygen == min( data_new$Oxygen ), ]
nrow( data_new ) #1678

#Depth
data_new <- data_new[!data_new$Mean.depth == max( data_new$Mean.depth ),]
nrow( data_new ) #1677

#Temp
data_new <- data_new[!data_new$Bio.Temp == min( data_new$Bio.Temp ), ]
data_new <- data_new[!data_new$Bio.Temp == min( data_new$Bio.Temp ), ]
data_new <- data_new[!data_new$Bio.Temp == min( data_new$Bio.Temp ), ]
nrow( data_new ) #1673 datapoints for snapper

######## Checking  multi-collinearity of environmental covariates
data_1  = data_new 
dim( data_1 )
data_1 <- na.omit( data_1 ) 
dim( data_1 )
bXY = geoXY( latitude = data_1$Lat, longitude = data_1$Lon, 
             lat0 = mean( data_1$Lat ), lon0 = mean( data_1$Lon ), unit = 1000 )
data_1 = cbind( data_1, bXY )

matrix_1 <- rcorr( as.matrix( data_1[,c( "Oxygen", "Bio.Temp", "X", "Y", "Mean.depth" )], type = "pearson" ) )
corrplot( matrix_1$r, type = "lower", tl.col = "black", method = "number", title = "Snapper", 
          mar = c( 0, 0, 1, 0 ), p.mat = matrix_1$P, sig.level = 0.05) #### "X" indicates that the correlation is not significant
#Depth and Oxygen collinear (-0.81)


######## Check the distribution of max length data 
species_length = data_1$Max.length
summary( species_length )
hist( species_length ) #### It would be reasonable to fit a Gaussian GAM for the species
species_g = fitdist( species_length, distr = "gamma", method = "mle", lower = c( 0, 0 ) )
species_n = fitdist( species_length, distr = "norm", method = "mle", lower = c( 0, 0 ) )
species_l = fitdist( species_length, distr = "lnorm", method = "mle", lower = c( 0, 0 ) )
par( mfrow = c( 2, 2 ) )
plot.legend <- c( "Gamma distribution", "Lognormal distribution", "Normal distribution" )
denscomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
cdfcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
qqcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
ppcomp( list( species_g, species_n, species_l), legendtext = plot.legend )
gofstat( list( species_g, species_n, species_l ) ) 

#Converting dissolved oxygen concentration from mol.m^-3 to mg/L
data_1 = data_1 %>% mutate(oxygen_mgL = Oxygen*31.9988/1000)

#Full GAM
gam_max_length.sn.1 = gam( Max.length ~ te( X, Y ) + s( Bio.Temp, k = 4, bs = "ts" ) + s( oxygen_mgL, k = 4, bs = "ts" ), 
                      Gamma( link = "log" ), method = 'REML', data = data_1 )
summary( gam_max_length.sn.1 ) #### R-sq.(adj) =  0.281; Deviance explained = 20.8%
#draw ( gam_max_length.sn.1 )


#removing oxygen (not significant), re-running GAM
gam_max_length.sn = gam( Max.length ~ te( X, Y ) + s( Bio.Temp, k = 4, bs = "ts" ), 
                         Gamma( link = "log" ), method = 'REML', data = data_1 )
summary( gam_max_length.sn ) #### R-sq.(adj) =  0.281; Deviance explained = 20.8%
#draw ( gam_max_length.sn )

# AIC
logLik.gam(gam_max_length.sn) # -6228.759 (df=16.72021)
logLik.gam(gam_max_length.sn.1) # -6228.759 (df=16.7219)

snapper = data_1

######## Estimate the indices of relative importance of the predictors
#### Estimate the indices of relative importance of the predictors for the GAM fitted to max length data
model_list = list( gam_max_length.sn )

source( "AssessRelativeImportancePredictors.R" )
RIeval_mean_length = relImportantVar( model_list[[1]], data_1, nRep = 10 )
RIeval_mean_length_2 = RIeval_mean_length[,c( 1 : 2, 4 )]
source( "CustomizedBarPlot.R" )
RIeval_mean_length_2 
rownames( RIeval_mean_length_2 ) = c( "Temperature", "Eastings", "Northings")
RIeval_mean_length_2 <- as.data.frame( RIeval_mean_length_2 )
RIeval_mean_length_2 <- RIeval_mean_length_2[order( - RIeval_mean_length_2$meanInd ),]

tiff( file = "SNA_max_length_GAM_IndicesOfRelativeImportance.tif", width = 3.5, height = 2.5, units = "in", 
      res = 600, compression = "lzw" ) 
par( cex = 0.5 )
custBarChart( RIeval_mean_length_2, confGiven = TRUE, col.bar = 'gray', col.se = 1, ylab = '', xlab = '' )
mtext( side = 1, text = "Predictor", cex = 0.7, line = 1.5, padj = 0.5 )
mtext( side = 2, text = "Index of relative importance", cex = 0.7, line = 1.5, padj = -1 )
dev.off()
save( "RIeval_mean_length", "RIeval_mean_length_2", file = "SNA_max_length_GAM_IndicesOfRelativeImportance.RData" )

######## Evaluate GAMs using Leave Group Out Cross Validation
#### Evaluate the GAM fitted to max length data
model_list = list( gam_max_length.sn )

source( "ModelEvaluation.R" )
Evaluation_GAM_max_length_SNA = evaluateModel( model_list, nIter = 10, confInt = TRUE ) 

save( "Evaluation_GAM_max_length_SNA", file = "Evaluation_GAM_max_length_SNA.RData" )

####################################################

#Hoki
length_data = read.csv( "HOK_trawl_lengths_DATES.csv" )
length_data = length_data %>%
  mutate(lgth_SL = (lgth/1.068)) # TOTAL LENGTH to STANDARD LENGTH
# Conversion from Kloster et al. 2011

dim( length_data )
nrow( length_data ) #### 271168 lines
names( length_data )
table( length_data$trip_code ) 
summary( length_data$station_no ) #### 1 NA
summary( length_data$dlon_s ) #### 10 NA's
summary( length_data$dlat_s ) #### 10 NA's
summary( length_data$dlon_e ) #### 125 NA's
summary( length_data$dlat_s ) #### 10 NA's
summary( length_data$min_gdepth ) #### 9728 NA's
summary( length_data$max_gdepth ) #### 9728 NA's
summary( length_data$bot_temp ) #### 59834 NA's
table( length_data$species ) 
summary( length_data$lgth_SL ) #### 1 NA
summary( length_data$no_a ) #### 3 NA's
summary( length_data$no_m ) #### 84 NA's
summary( length_data$no_f ) #### 83 NA's

length_data_1 = drop_na( length_data )
nrow( length_data_1 ) #### 204409 lines remaining
length_data_2 = length_data[!is.na( length_data$lgth_SL ) & !is.na( length_data$dlon_s ) &
                              !is.na( length_data$dlat_s ) & !is.na( length_data$bot_temp ) &
                              !is.na( length_data$min_gdepth ) & !is.na( length_data$max_gdepth ),]
nrow( length_data_2 ) #### 204503 lines remaining
length_data <- length_data_2 

plot( length_data$dlon_s, length_data$dlat_s )

#Trip station
summary( length_data$station_no ) 
length( unique( length_data$station_no ) )
table( length_data$station_no ) 
trip_station <- paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) )
length( unique( trip_station ) )
length_data$trip_station <- as.factor( paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) ) )
dim( length_data )
names( length_data )

#Mean depth for each station
length_data = length_data %>%
  add_column(mean_depth = ((length_data$min_gdepth + length_data$max_gdepth)/2))

#Mean length data for each station
mean_length_data <- aggregate( lgth_SL ~ trip_station, data = length_data, mean )
mean_length_data$dlon_s <- length_data$dlon_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$dlat_s <- length_data$dlat_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$bot_temp <- length_data$bot_temp[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$mean_depth <- length_data$mean_depth[match( mean_length_data$trip_station, length_data$trip_station)]
mean_length_data$date_s <- length_data$date_s[match( mean_length_data$trip_station,length_data$trip_station)]
dim( mean_length_data )
names( mean_length_data )
thedata = mean_length_data
plot( thedata$dlon_s, thedata$dlat_s )

#adding 95th quantile value of length for each Trip station
length_data = as_tibble(length_data)

quant_length_data = length_data %>%
  dplyr::select(date_s, trip_station, lgth_SL, dlon_s, dlat_s, bot_temp)

quant_data <-  data.frame(quant_length_data$date_s, quant_length_data$trip_station, quant_length_data$bot_temp, quant_length_data$dlon_s, quant_length_data$dlat_s, quant_length_data$lgth_SL )
names( quant_data ) <- c("Date", "Trip.Station",  "Situ.Temp", "Lon", "Lat", "Length" ) 

bXY = geoXY( latitude = quant_data$Lat, longitude = quant_data$Lon, 
             lat0 = mean( quant_data$Lat ), lon0 = mean( quant_data$Lon ), unit = 1000 )
quant_data = cbind( quant_data, bXY )

sort( unique( quant_data$Trip.Station ) )

#set stations as factor
quant_data$Trip.Station <- as.factor( paste0( as.factor( quant_data$Trip.Station )))

#here are the 95 quantile values for each trip station
quant_dat = quant_data %>%
  group_by(Trip.Station) %>%
  summarise(enframe(quantile(Length, c(0.95)), "quantile", "95.Length"))

quant_95=as_tibble(quant_dat$`95.Length`)


######## Load in environmental data (mean temperature.bio, and oxygen at maximum depth)
######## Those data were extracted from https://www.bio-oracle.org/downloads-to-email.php

ox = read.asciigrid( "oxygen.asc" )
ox = raster( ox )
x1 <- crop( ox, extent( -180, 0, -90, 90 ) )
x2 <- crop( ox, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
ox <- merge( x1, x2 )
plot( ox )
temp = read.asciigrid( "temp.asc" )
temp = raster( temp )
x1 <- crop( temp, extent( -180, 0, -90, 90 ) )
x2 <- crop( temp, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
temp <- merge( x1, x2 )
plot( temp )

predictors <- stack( temp, ox )

######## Combine environmental data with monitoring data and check for outliers in the response and covariate variables
#### Combine environmental data with monitoring data
data_new <-  data.frame( raster::extract( predictors, thedata[3:4] ),
                         thedata$date_s, thedata$bot_temp, thedata$mean_depth, thedata$dlon_s, thedata$dlat_s, thedata$lgth_SL )

#adding 95 quantiles
data_new = cbind(data_new, quant_95)
names( data_new ) <- c( "Bio.Temp", "Oxygen", "Date", "Situ.Temp", "Mean.depth", "Lon", "Lat", "Length", "Max.length" ) 
nrow( data_new ) ## 9422 lines 


#### Check for outliers in the response variable
summary( data_new$Max.length )
hist( data_new$Max.length )
boxplot( data_new$Max.length )
sort( unique( data_new$Max.length ) )

#### Check for outliers in the bottom temperatures
summary( data_new$Bio.Temp )
hist( data_new$Bio.Temp )
boxplot( data_new$Bio.Temp )
sort( unique( data_new$Bio.Temp ) )

#### Check for outliers in the bottom dissolved oxygen (DO) concentrations
summary( data_new$Oxygen )
hist( data_new$Oxygen )
boxplot( data_new$Oxygen )
sort( unique( data_new$Oxygen ) )

#### Check for outliers in the mean depth
summary( data_new$Mean.depth )
hist( data_new$Mean.depth )
boxplot( data_new$Mean.depth )
sort( unique( data_new$Mean.depth ) )

## The diagnostics suggest that maximum depth has 2 outliers,
#therefore they should be removed
nrow( data_new ) #9422
data_new <- data_new[!data_new$Mean.depth == max( data_new$Mean.depth), ]
nrow( data_new ) #9421
data_new <- data_new[!data_new$Mean.depth == max( data_new$Mean.depth), ]
nrow( data_new ) #9420 datapoints for hoki

######## Checking  multi-collinearity of environmental covariates
data_1  = data_new 
dim( data_1 )
data_1 <- na.omit( data_1 ) 
dim( data_1 ) #9420
bXY = geoXY( latitude = data_1$Lat, longitude = data_1$Lon, 
             lat0 = mean( data_1$Lat ), lon0 = mean( data_1$Lon ), unit = 1000 )
data_1 = cbind( data_1, bXY )
matrix_1 <- rcorr( as.matrix( data_1[,c( "Oxygen", "Bio.Temp", "X", "Y", "Mean.depth" )], type = "pearson" ) )
corrplot( matrix_1$r, type = "lower", tl.col = "black", method = "number", title = "Hoki", 
          mar = c( 0, 0, 1, 0 ), p.mat = matrix_1$P, sig.level = 0.05) #### "X" indicates that the correlation is not significant
# Depth collinear with bio.temperature (-0.79)

######## Check the distribution of max length data 
species_length = data_1$Max.length
summary( species_length )
hist( species_length ) #### It would be reasonable to fit a Gaussian GAM for the species
species_g = fitdist( species_length, distr = "gamma", method = "mle", lower = c( 0, 0 ) )
species_n = fitdist( species_length, distr = "norm", method = "mle", lower = c( 0, 0 ) )
species_l = fitdist( species_length, distr = "lnorm", method = "mle", lower = c( 0, 0 ) )
par( mfrow = c( 2, 2 ) )
plot.legend <- c( "Gamma distribution", "Lognormal distribution", "Normal distribution" )
denscomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
cdfcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
qqcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
ppcomp( list( species_g, species_n, species_l), legendtext = plot.legend )
gofstat( list( species_g, species_n, species_l ) ) 


#Converting dissolved oxygen concentration from mol.m^-3 to mg/L
data_1 = data_1 %>% mutate(oxygen_mgL = Oxygen*31.9988/1000)

#GAM
gam_max_length.h = gam( Max.length ~ te( X, Y ) + s( Bio.Temp, k = 4, bs = "ts" )+ s( oxygen_mgL, k = 4, bs = "ts" ), 
                        Gamma( link = "log" ), method = 'REML', data = data_1 )
summary( gam_max_length.h ) #R-sq.(adj) =  0.461   Deviance explained = 41.9%

hoki = data_1

######## Estimate the indices of relative importance of the predictors
#### Estimate the indices of relative importance of the predictors for the GAM fitted to max length data
model_list = list( gam_max_length.h )

source( "AssessRelativeImportancePredictors.R" )
RIeval_mean_length = relImportantVar( model_list[[1]], data_1, nRep = 10 )
RIeval_mean_length_2 = RIeval_mean_length[,c( 1 : 2, 4 )]
source( "CustomizedBarPlot.R" )
RIeval_mean_length_2 
rownames( RIeval_mean_length_2 ) = c( "Temperature", "Eastings", "Northings",  "Oxygen" )
RIeval_mean_length_2 <- as.data.frame( RIeval_mean_length_2 )
RIeval_mean_length_2 <- RIeval_mean_length_2[order( - RIeval_mean_length_2$meanInd ),]

tiff( file = "HOK_max_length_GAM_IndicesOfRelativeImportance.tif", width = 3.5, height = 2.5, units = "in", 
      res = 600, compression = "lzw" ) 
par( cex = 0.5 )
custBarChart( RIeval_mean_length_2, confGiven = TRUE, col.bar = 'gray', col.se = 1, ylab = '', xlab = '' )
mtext( side = 1, text = "Predictor", cex = 0.7, line = 1.5, padj = 0.5 )
mtext( side = 2, text = "Index of relative importance", cex = 0.7, line = 1.5, padj = -1 )
dev.off()

save( "RIeval_mean_length", "RIeval_mean_length_2", file = "HOK_max_length_GAM_IndicesOfRelativeImportance.RData" )

######## Evaluate GAMs using Leave Group Out Cross Validation
#### Evaluate the GAM fitted to max length data
model_list = list( gam_max_length.h )

source( "ModelEvaluation.R" )
Evaluation_GAM_max_length_HOK = evaluateModel( model_list, nIter = 10, confInt = TRUE ) 

save( "Evaluation_GAM_max_length_HOK", file = "Evaluation_GAM_max_length_HOK.RData" )

####################################################

#Southern blue whiting
length_data = read.csv( "SBW_trawl_lengths_DATES.csv" )
length_data = length_data %>%
  mutate(lgth_SL = (lgth*0.941)) # FORK LENGTH to STANDARD LENGTH
# Conversion from Cohen et al. 1990

dim( length_data )
nrow( length_data ) #### 20886 lines
names( length_data )
table( length_data$trip_code ) 
summary( length_data$station_no ) #### 1 NA
summary( length_data$dlon_s ) #### 1 NA
summary( length_data$dlat_s ) #### 1 NA
summary( length_data$dlon_e ) #### 1 NA
summary( length_data$dlat_s ) #### 1 NA
summary( length_data$min_gdepth ) #### 680 NA's
summary( length_data$max_gdepth ) #### 680 NA's
summary( length_data$bot_temp ) #### 4393 NA's
table( length_data$species ) 
summary( length_data$lgth_SL ) #### 1 NA
summary( length_data$no_a ) #### 1 NA
summary( length_data$no_m ) #### 1 NA
summary( length_data$no_f ) #### 1 NA

length_data_1 = drop_na( length_data )
nrow( length_data_1 ) #### 16025 lines remaining
length_data_2 = length_data[!is.na( length_data$lgth_SL ) & !is.na( length_data$dlon_s ) &
                              !is.na( length_data$dlat_s ) & !is.na( length_data$bot_temp ) &
                              !is.na( length_data$min_gdepth ) & !is.na( length_data$max_gdepth ),]
nrow( length_data_2 ) #### 16025 lines remaining
length_data <- length_data_2 

plot( length_data$dlon_s, length_data$dlat_s )

#Trip station
summary( length_data$station_no ) 
length( unique( length_data$station_no ) )
table( length_data$station_no ) 
trip_station <- paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) )
length( unique( trip_station ) )
length_data$trip_station <- as.factor( paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) ) )
dim( length_data )
names( length_data )

#Mean depth of station
length_data = length_data %>%
  add_column(mean_depth = ((length_data$min_gdepth + length_data$max_gdepth)/2))

#Mean length at station
mean_length_data <- aggregate( lgth_SL ~ trip_station, data = length_data, mean )
mean_length_data$dlon_s <- length_data$dlon_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$dlat_s <- length_data$dlat_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$bot_temp <- length_data$bot_temp[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$mean_depth <- length_data$mean_depth[match( mean_length_data$trip_station, length_data$trip_station)]
mean_length_data$date_s <- length_data$date_s[match( mean_length_data$trip_station,length_data$trip_station)]
dim( mean_length_data )
names( mean_length_data )
thedata = mean_length_data
plot( thedata$dlon_s, thedata$dlat_s )

#adding 95th quantile value of length for each Trip station
length_data = as_tibble(length_data)

quant_length_data = length_data %>%
  dplyr::select(date_s, trip_station, lgth_SL, dlon_s, dlat_s, bot_temp)

quant_data <-  data.frame(quant_length_data$date_s, quant_length_data$trip_station, quant_length_data$bot_temp, quant_length_data$dlon_s, quant_length_data$dlat_s, quant_length_data$lgth_SL )
names( quant_data ) <- c("Date", "Trip.Station",  "Situ.Temp", "Lon", "Lat", "Length" ) 

bXY = geoXY( latitude = quant_data$Lat, longitude = quant_data$Lon, 
             lat0 = mean( quant_data$Lat ), lon0 = mean( quant_data$Lon ), unit = 1000 )
quant_data = cbind( quant_data, bXY )

sort( unique( quant_data$Trip.Station ) )

#set stations as factor
quant_data$Trip.Station <- as.factor( paste0( as.factor( quant_data$Trip.Station )))

#here are the 95 quantile values for each trip station
quant_dat = quant_data %>%
  group_by(Trip.Station) %>%
  summarise(enframe(quantile(Length, c(0.95)), "quantile", "95.Length"))

quant_95=as_tibble(quant_dat$`95.Length`)

######## Load in environmental data (mean temperature, and oxygen at maximum depth)
######## Those data were extracted from https://www.bio-oracle.org/downloads-to-email.php

ox = read.asciigrid( "oxygen.asc" )
ox = raster( ox )
x1 <- crop( ox, extent( -180, 0, -90, 90 ) )
x2 <- crop( ox, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
ox <- merge( x1, x2 )
plot( ox )
temp = read.asciigrid( "temp.asc" )
temp = raster( temp )
x1 <- crop( temp, extent( -180, 0, -90, 90 ) )
x2 <- crop( temp, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
temp <- merge( x1, x2 )
plot( temp )

predictors <- stack( temp, ox )

######## Combine environmental data with monitoring data and check for outliers in the response and covariate variables
#### Combine environmental data with monitoring data
data_new <-  data.frame( raster::extract( predictors, thedata[3:4] ),
                         thedata$date_s, thedata$bot_temp, thedata$mean_depth, thedata$dlon_s, thedata$dlat_s, thedata$lgth_SL )

#adding 95 quantiles
data_new = cbind(data_new, quant_95)
names( data_new ) <- c( "Bio.Temp", "Oxygen", "Date", "Situ.Temp", "Mean.depth", "Lon", "Lat", "Length", "Max.length" ) 
nrow( data_new ) ## 1361 


#### Check for outliers in the response variable
summary( data_new$Max.length )
hist( data_new$Max.length )
boxplot( data_new$Max.length )
sort( unique( data_new$Max.length ) )

#### Check for outliers in the bottom temperatures
summary( data_new$Bio.Temp )
hist( data_new$Bio.Temp )
boxplot( data_new$Bio.Temp )
sort( unique( data_new$Bio.Temp ) )

#### Check for outliers in the bottom dissolved oxygen (DO) concentrations
summary( data_new$Oxygen )
hist( data_new$Oxygen )
boxplot( data_new$Oxygen )
sort( unique( data_new$Oxygen ) )

#### Check for outliers in the mean depth
summary( data_new$Mean.depth )
hist( data_new$Mean.depth )
boxplot( data_new$Mean.depth )
sort( unique( data_new$Mean.depth ) )

## The diagnostics suggest that the minimum DO is an outlier, plus maximum of mean depth,
## and should, therefore, be removed
data_new <- data_new[!data_new$Oxygen == min( data_new$Oxygen ), ] 
nrow( data_new ) ## 1360
data_new <- data_new[!data_new$Mean.depth == max( data_new$Mean.depth), ]
nrow( data_new ) ## 1359
data_new <- data_new[!data_new$Mean.depth == max( data_new$Mean.depth), ]
nrow( data_new ) ## 1358

data_1  = data_new 
dim( data_1 ) ## 1358
data_1 <- na.omit( data_1 ) 
dim( data_1 ) ## 1358

boxplot( data_1$Bio.Temp )
data_1 <- data_1 %>% filter(Bio.Temp > 5, Bio.Temp < 8, Oxygen > 230)
nrow( data_1 ) #1312
boxplot(data_1$Oxygen)
data_1 <- data_1[!data_1$Oxygen == min( data_1$Oxygen ), ]
data_1 <- data_1[!data_1$Oxygen == min( data_1$Oxygen ), ]
boxplot(data_1$Mean.depth)
data_1 <- data_1[!data_1$Mean.depth == min( data_1$Mean.depth ), ]
data_1 <- data_1[!data_1$Mean.depth == min( data_1$Mean.depth ), ]
nrow( data_1 ) #1308 datapoints for SBW

######## Check collinearity 
bXY = geoXY( latitude = data_1$Lat, longitude = data_1$Lon, 
             lat0 = mean( data_1$Lat ), lon0 = mean( data_1$Lon ), unit = 1000 )
data_1 = cbind( data_1, bXY )
matrix_1 <- rcorr( as.matrix( data_1[,c( "Oxygen", "Bio.Temp", "X", "Y", "Mean.depth" )], type = "pearson" ) )
( matrix_1$r ) 
corrplot( matrix_1$r, type = "lower", tl.col = "black", method = "number", title = "Southern blue whiting", 
          mar = c( 0, 0, 1, 0 ), p.mat = matrix_1$P, sig.level = 0.05) #### "X" indicates that the correlation is not significant

######## Check the distribution of max length data 
species_length = data_1$Max.length
summary( species_length )
hist( species_length ) #### It would be reasonable to fit a Gaussian GAM for the species
species_g = fitdist( species_length, distr = "gamma", method = "mle", lower = c( 0, 0 ) )
species_n = fitdist( species_length, distr = "norm", method = "mle", lower = c( 0, 0 ) )
species_l = fitdist( species_length, distr = "lnorm", method = "mle", lower = c( 0, 0 ) )
par( mfrow = c( 2, 2 ) )
plot.legend <- c( "Gamma distribution", "Lognormal distribution", "Normal distribution" )
denscomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
cdfcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
qqcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
ppcomp( list( species_g, species_n, species_l), legendtext = plot.legend )
gofstat( list( species_g, species_n, species_l ) ) 


#Converting dissolved oxygen concentration from mol.m^-3 to mg/L
data_1 = data_1 %>% mutate(oxygen_mgL = Oxygen*31.9988/1000)

#GAM
gam_max_length.sb = gam( Max.length ~ te( X, Y ) + s( Bio.Temp, k = 4, bs = "ts" ) + s( oxygen_mgL, k = 4, bs = "ts" ), 
                         Gamma( link = "log" ), method = 'REML', data = data_1 )
summary( gam_max_length.sb ) #R-sq.(adj) =  0.229   Deviance explained = 21.6%

sbw = data_1

######## Estimate the indices of relative importance of the predictors
#### Estimate the indices of relative importance of the predictors for the GAM fitted to mean length data
model_list = list( gam_max_length.sb )

source( "AssessRelativeImportancePredictors.R" )
RIeval_mean_length = relImportantVar( model_list[[1]], data_1, nRep = 10 )
RIeval_mean_length_2 = RIeval_mean_length[,c( 1 : 2, 4 )]
source( "CustomizedBarPlot.R" )
RIeval_mean_length_2 
rownames( RIeval_mean_length_2 ) = c( "Temperature", "Eastings", "Northings", "Oxygen" )
RIeval_mean_length_2 <- as.data.frame( RIeval_mean_length_2 )
RIeval_mean_length_2 <- RIeval_mean_length_2[order( - RIeval_mean_length_2$meanInd ),]

tiff( file = "SBW_max_length_GAM_IndicesOfRelativeImportance.tif", width = 3.5, height = 2.5, units = "in", 
      res = 600, compression = "lzw" ) 
par( cex = 0.5 )
custBarChart( RIeval_mean_length_2, confGiven = TRUE, col.bar = 'gray', col.se = 1, ylab = '', xlab = '' )
mtext( side = 1, text = "Predictor", cex = 0.7, line = 1.5, padj = 0.5 )
mtext( side = 2, text = "Index of relative importance", cex = 0.7, line = 1.5, padj = -1 )
dev.off()

save( "RIeval_mean_length", "RIeval_mean_length_2", file = "SBW_max_length_GAM_IndicesOfRelativeImportance.RData" )

######## Evaluate GAMs using Leave Group Out Cross Validation
#### Evaluate the GAM fitted to max length data
model_list = list( gam_max_length.sb )

source( "ModelEvaluation.R" )
Evaluation_GAM_max_length_SBW = evaluateModel( model_list, nIter = 10, confInt = TRUE ) 

save( "Evaluation_GAM_max_length_SBW", file = "Evaluation_GAM_max_length_SBW.RData" )

####################################################

#Orange roughy
length_data = read.csv( "ORH_trawl_lengths_DATES.csv" ) #SL

dim( length_data )
nrow( length_data ) #### 140566 lines
names( length_data )
table( length_data$trip_code ) 
summary( length_data$station_no ) #### 1 NA
summary( length_data$dlon_s ) #### 28 NA
summary( length_data$dlat_s ) #### 28 NA
summary( length_data$dlon_e ) #### 180 NA
summary( length_data$dlat_s ) #### 28 NA
summary( length_data$min_gdepth ) #### 4454 NA's
summary( length_data$max_gdepth ) #### 4479 NA's
summary( length_data$bot_temp ) #### 49470 NA's
table( length_data$species ) 
summary( length_data$lgth ) #### 1 NA
summary( length_data$no_a ) #### 1 NA
summary( length_data$no_m ) #### 17 NA
summary( length_data$no_f ) #### 8 NA

length_data_1 = drop_na( length_data )
nrow( length_data_1 ) #### 90172 lines remaining
length_data_2 = length_data[!is.na( length_data$lgth ) & !is.na( length_data$dlon_s ) &
                              !is.na( length_data$dlat_s ) & !is.na( length_data$bot_temp ) &
                              !is.na( length_data$min_gdepth ) & !is.na( length_data$max_gdepth ),]
nrow( length_data_2 ) #### 90219 lines remaining
length_data <- length_data_2 

plot( length_data$dlon_s, length_data$dlat_s )

#Trip station
summary( length_data$station_no ) 
length( unique( length_data$station_no ) )
table( length_data$station_no ) 
trip_station <- paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) )
length( unique( trip_station ) )
length_data$trip_station <- as.factor( paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) ) )
dim( length_data )
names( length_data )

#Mean depth of station
length_data = length_data %>%
  add_column(mean_depth = ((length_data$min_gdepth + length_data$max_gdepth)/2))

#Mean length at station
mean_length_data <- aggregate( lgth ~ trip_station, data = length_data, mean )
mean_length_data$dlon_s <- length_data$dlon_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$dlat_s <- length_data$dlat_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$bot_temp <- length_data$bot_temp[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$mean_depth <- length_data$mean_depth[match( mean_length_data$trip_station, length_data$trip_station)]
mean_length_data$date_s <- length_data$date_s[match( mean_length_data$trip_station,length_data$trip_station)]
dim( mean_length_data )
names( mean_length_data )
thedata = mean_length_data
plot( thedata$dlon_s, thedata$dlat_s )

#adding 95th quantile value of length for each Trip station
length_data = as_tibble(length_data)

quant_length_data = length_data %>%
  dplyr::select(date_s, trip_station, lgth, dlon_s, dlat_s, bot_temp)

quant_data <-  data.frame(quant_length_data$date_s, quant_length_data$trip_station, quant_length_data$bot_temp, quant_length_data$dlon_s, quant_length_data$dlat_s, quant_length_data$lgth )
names( quant_data ) <- c("Date", "Trip.Station",  "Situ.Temp", "Lon", "Lat", "Length" ) 

bXY = geoXY( latitude = quant_data$Lat, longitude = quant_data$Lon, 
             lat0 = mean( quant_data$Lat ), lon0 = mean( quant_data$Lon ), unit = 1000 )
quant_data = cbind( quant_data, bXY )

sort( unique( quant_data$Trip.Station ) )


#set stations as factor
quant_data$Trip.Station <- as.factor( paste0( as.factor( quant_data$Trip.Station )))

#here are the 95 quantile values for each trip station
quant_dat = quant_data %>%
  group_by(Trip.Station) %>%
  summarise(enframe(quantile(Length, c(0.95)), "quantile", "95.Length"))

quant_95=as_tibble(quant_dat$`95.Length`)

######## Load in environmental data (mean temperature, and oxygen at maximum depth)
######## Those data were extracted from https://www.bio-oracle.org/downloads-to-email.php

ox = read.asciigrid( "oxygen.asc" )
ox = raster( ox )
x1 <- crop( ox, extent( -180, 0, -90, 90 ) )
x2 <- crop( ox, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
ox <- merge( x1, x2 )
plot( ox )
temp = read.asciigrid( "temp.asc" )
temp = raster( temp )
x1 <- crop( temp, extent( -180, 0, -90, 90 ) )
x2 <- crop( temp, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
temp <- merge( x1, x2 )
plot( temp )

predictors <- stack( temp, ox )

######## Combine environmental data with monitoring data and check for outliers in the response and covariate variables
#### Combine environmental data with monitoring data
data_new <-  data.frame( raster::extract( predictors, thedata[3:4] ),
                         thedata$date_s ,thedata$bot_temp, thedata$mean_depth, thedata$dlon_s, thedata$dlat_s, thedata$lgth )

#adding 95 quantiles
data_new = cbind(data_new, quant_95)
names( data_new ) <- c( "Bio.Temp", "Oxygen", "Date", "Situ.Temp", "Mean.depth", "Lon", "Lat", "Length", "Max.length" ) 
nrow( data_new ) ## 6767 lines 


#### Check for outliers in the response variable
summary( data_new$Max.length )
hist( data_new$Max.length )
boxplot( data_new$Max.length )
sort( unique( data_new$Max.length ) )

#### Check for outliers in the bottom temperatures
summary( data_new$Bio.Temp )
hist( data_new$Bio.Temp )
boxplot( data_new$Bio.Temp )
sort( unique( data_new$Bio.Temp ) )

#### Check for outliers in the bottom dissolved oxygen (DO) concentrations
summary( data_new$Oxygen )
hist( data_new$Oxygen )
boxplot( data_new$Oxygen )
sort( unique( data_new$Oxygen ) )

#### Check for outliers in the mean depth
summary( data_new$Mean.depth )
hist( data_new$Mean.depth )
boxplot( data_new$Mean.depth )
sort( unique( data_new$Mean.depth ) )

## The diagnostics suggest that the maximum temp and maximum oxygen have several outliers, plus maximum of mean depth,  and should, therefore,
## be removed
nrow( data_new ) ## 6767
data_new <- data_new[!data_new$Mean.depth == max( data_new$Mean.depth), ]
nrow( data_new ) ## 6766
data_new <- data_new[!data_new$Oxygen == max( data_new$Oxygen), ]
nrow( data_new ) ## 6765
data_new <- data_new[!data_new$Oxygen == max( data_new$Oxygen), ]
nrow( data_new ) ## 6764
data_new <- data_new[!data_new$Oxygen == max( data_new$Oxygen), ]
nrow( data_new ) ## 6763
boxplot( data_new$Bio.Temp )
nrow( data_new ) #6763 datapoints for orange roughy

######## Check collinearity 
data_1  = data_new 
dim( data_1 )
data_1 <- na.omit( data_1 ) 
dim( data_1 )
bXY = geoXY( latitude = data_1$Lat, longitude = data_1$Lon, 
             lat0 = mean( data_1$Lat ), lon0 = mean( data_1$Lon ), unit = 1000 )
data_1 = cbind( data_1, bXY )

matrix_1 <- rcorr( as.matrix( data_1[,c( "Oxygen", "Bio.Temp", "X", "Y", "Mean.depth" )], type = "pearson" ) )
( matrix_1$r ) 
corrplot( matrix_1$r, type = "lower", tl.col = "black", method = "number", title = "Orange roughy", 
          mar = c( 0, 0, 1, 0 ), p.mat = matrix_1$P, sig.level = 0.05) #### "X" indicates that the correlation is not significant

######## Check the distribution of max length data 
species_length = data_1$Max.length
summary( species_length )
hist( species_length ) #### It would be reasonable to fit a Gaussian GAM for the species
species_g = fitdist( species_length, distr = "gamma", method = "mle", lower = c( 0, 0 ) )
species_n = fitdist( species_length, distr = "norm", method = "mle", lower = c( 0, 0 ) )
species_l = fitdist( species_length, distr = "lnorm", method = "mle", lower = c( 0, 0 ) )
par( mfrow = c( 2, 2 ) )
plot.legend <- c( "Gamma distribution", "Lognormal distribution", "Normal distribution" )
denscomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
cdfcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
qqcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
ppcomp( list( species_g, species_n, species_l), legendtext = plot.legend )
gofstat( list( species_g, species_n, species_l ) ) 


#Converting oxygen
data_1 = data_1 %>% mutate(oxygen_mgL = Oxygen*31.9988/1000)

# GAM
gam_max_length.or = gam( Max.length ~ te( X, Y ) + s( Bio.Temp, k = 4, bs = "ts" ) + s( oxygen_mgL, k = 4, bs = "ts" ), 
                         Gamma( link = "log" ), method = 'REML', data = data_1 )
summary( gam_max_length.or ) #R-sq.(adj) =  0.131   Deviance explained = 9.94%


roughy = data_1


######## Estimate the indices of relative importance of the predictors
#### Estimate the indices of relative importance of the predictors for the GAM fitted to max length data
model_list = list( gam_max_length.or )

source( "AssessRelativeImportancePredictors.R" )
RIeval_mean_length = relImportantVar( model_list[[1]], data_1, nRep = 10 )
RIeval_mean_length_2 = RIeval_mean_length[,c( 1 : 2, 4 )]
source( "CustomizedBarPlot.R" )
RIeval_mean_length_2 
rownames( RIeval_mean_length_2 ) = c( "Temperature", "Eastings", "Northings", "Oxygen" )
RIeval_mean_length_2 <- as.data.frame( RIeval_mean_length_2 )
RIeval_mean_length_2 <- RIeval_mean_length_2[order( - RIeval_mean_length_2$meanInd ),]

tiff( file = "ORH_max_length_GAM_IndicesOfRelativeImportance.tif", width = 3.5, height = 2.5, units = "in", 
      res = 600, compression = "lzw" ) 
par( cex = 0.5 )
custBarChart( RIeval_mean_length_2, confGiven = TRUE, col.bar = 'gray', col.se = 1, ylab = '', xlab = '' )
mtext( side = 1, text = "Predictor", cex = 0.7, line = 1.5, padj = 0.5 )
mtext( side = 2, text = "Index of relative importance", cex = 0.7, line = 1.5, padj = -1 )
dev.off()

save( "RIeval_mean_length", "RIeval_mean_length_2", file = "ORH_max_length_GAM_IndicesOfRelativeImportance.RData" )

######## Evaluate GAMs using Leave Group Out Cross Validation
#### Evaluate the GAM fitted to max length data
model_list = list( gam_max_length.or )

source( "ModelEvaluation.R" )
Evaluation_GAM_max_length_ORH = evaluateModel( model_list, nIter = 10, confInt = TRUE ) 

save( "Evaluation_GAM_max_length_ORH", file = "Evaluation_GAM_max_length_ORH.RData" )

####################################################

#New Zealand arrow squid
length_data = read.csv( "NOS_trawl_lengths_DATES.csv" ) #Mantle length

dim( length_data )
nrow( length_data ) #### 23910 lines
names( length_data )
table( length_data$trip_code ) 
summary( length_data$station_no ) #### 1 NA
summary( length_data$dlon_s ) #### 1 NA
summary( length_data$dlat_s ) #### 1 NA
summary( length_data$dlon_e ) #### 16 NA
summary( length_data$dlat_s ) #### 1 NA
summary( length_data$min_gdepth ) #### 426 NA's
summary( length_data$max_gdepth ) #### 426 NA's
summary( length_data$bot_temp ) #### 3807 NA's
table( length_data$species ) 
summary( length_data$lgth ) #### 1 NA
summary( length_data$no_a ) #### 1 NA
summary( length_data$no_m ) #### 1 NA
summary( length_data$no_f ) #### 1 NA

length_data_1 = drop_na( length_data )
nrow( length_data_1 ) #### 19848 lines remaining
length_data_2 = length_data[!is.na( length_data$lgth ) & !is.na( length_data$dlon_s ) &
                              !is.na( length_data$dlat_s ) & !is.na( length_data$bot_temp ) &
                              !is.na( length_data$min_gdepth ) & !is.na( length_data$max_gdepth ),]
nrow( length_data_2 ) #### 19849 lines remaining
length_data <- length_data_2 

plot( length_data$dlon_s, length_data$dlat_s )

#Trip station
summary( length_data$station_no ) 
length( unique( length_data$station_no ) )
table( length_data$station_no ) 
trip_station <- paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) )
length( unique( trip_station ) )
length_data$trip_station <- as.factor( paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) ) )
dim( length_data )
names( length_data )

#Mean depth of station
length_data = length_data %>%
  add_column(mean_depth = ((length_data$min_gdepth + length_data$max_gdepth)/2))

#Mean length at station
mean_length_data <- aggregate( lgth ~ trip_station, data = length_data, mean )
mean_length_data$dlon_s <- length_data$dlon_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$dlat_s <- length_data$dlat_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$bot_temp <- length_data$bot_temp[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$mean_depth <- length_data$mean_depth[match( mean_length_data$trip_station, length_data$trip_station)]
mean_length_data$date_s <- length_data$date_s[match( mean_length_data$trip_station,length_data$trip_station)]
dim( mean_length_data )
names( mean_length_data )
thedata = mean_length_data
plot( thedata$dlon_s, thedata$dlat_s )

#adding 95th quantile value of length for each Trip station
length_data = as_tibble(length_data)

quant_length_data = length_data %>%
  dplyr::select(date_s, trip_station, lgth, dlon_s, dlat_s, bot_temp)

quant_data <-  data.frame(quant_length_data$date_s,quant_length_data$trip_station, quant_length_data$bot_temp, quant_length_data$dlon_s, quant_length_data$dlat_s, quant_length_data$lgth )
names( quant_data ) <- c("Date", "Trip.Station",  "Situ.Temp", "Lon", "Lat", "Length" ) 

bXY = geoXY( latitude = quant_data$Lat, longitude = quant_data$Lon, 
             lat0 = mean( quant_data$Lat ), lon0 = mean( quant_data$Lon ), unit = 1000 )
quant_data = cbind( quant_data, bXY )

sort( unique( quant_data$Trip.Station ) )

#set stations as factor
quant_data$Trip.Station <- as.factor( paste0( as.factor( quant_data$Trip.Station )))

#here are the 95 quantile values for each trip station
quant_dat = quant_data %>%
  group_by(Trip.Station) %>%
  summarise(enframe(quantile(Length, c(0.95)), "quantile", "95.Length"))

quant_95=as_tibble(quant_dat$`95.Length`)

######## Load in environmental data (mean temperature, and oxygen at maximum depth)
######## Those data were extracted from https://www.bio-oracle.org/downloads-to-email.php

ox = read.asciigrid( "oxygen.asc" )
ox = raster( ox )
x1 <- crop( ox, extent( -180, 0, -90, 90 ) )
x2 <- crop( ox, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
ox <- merge( x1, x2 )
plot( ox )
temp = read.asciigrid( "temp.asc" )
temp = raster( temp )
x1 <- crop( temp, extent( -180, 0, -90, 90 ) )
x2 <- crop( temp, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
temp <- merge( x1, x2 )
plot( temp )

predictors <- stack( temp, ox )

######## Combine environmental data with monitoring data and check for outliers in the response and covariate variables
#### Combine environmental data with monitoring data
data_new <-  data.frame( raster::extract( predictors, thedata[3:4] ),
                         thedata$date_s, thedata$bot_temp, thedata$mean_depth, thedata$dlon_s, thedata$dlat_s, thedata$lgth )

#adding 95 quantiles
data_new = cbind(data_new, quant_95)
names( data_new ) <- c( "Bio.Temp", "Oxygen", "Date", "Situ.Temp", "Mean.depth", "Lon", "Lat", "Length", "Max.length" ) 
nrow( data_new ) ## 2919 lines 


#### Check for outliers in the response variable
summary( data_new$Max.length )
hist( data_new$Max.length )
boxplot( data_new$Max.length )
sort( unique( data_new$Max.length ) )

#### Check for outliers in the bottom temperatures
summary( data_new$Bio.Temp )
hist( data_new$Bio.Temp )
boxplot( data_new$Bio.Temp )
sort( unique( data_new$Bio.Temp ) )

#### Check for outliers in the bottom dissolved oxygen (DO) concentrations
summary( data_new$Oxygen )
hist( data_new$Oxygen )
boxplot( data_new$Oxygen )
sort( unique( data_new$Oxygen ) )

#### Check for outliers in the mean depth
summary( data_new$Mean.depth )
hist( data_new$Mean.depth )
boxplot( data_new$Mean.depth )
sort( unique( data_new$Mean.depth ) )

## The diagnostics suggest that the maximum of mean length (x4), plus min of oxygen,
# and max depth are outliers, and should, therefore, be removed
nrow( data_new ) ## 2919
data_new <- data_new[!data_new$Max.length == max( data_new$Max.length ), ] 
data_new <- data_new[!data_new$Max.length == max( data_new$Max.length ), ] 
data_new <- data_new[!data_new$Max.length == max( data_new$Max.length ), ] 
nrow( data_new ) ## 2916
boxplot(data_new$Max.length)
data_new <- data_new[!data_new$Max.length == max( data_new$Max.length ), ] 
nrow( data_new ) ## 2915
boxplot(data_new$Max.length)

data_new <- data_new[!data_new$Oxygen == min( data_new$Oxygen ), ] 
nrow( data_new ) # 2913
data_new <- data_new[!data_new$Mean.depth == max( data_new$Mean.depth ), ] 
nrow( data_new ) ## 2912
boxplot(data_new$Mean.depth)

######## Check collinearity 
data_1  = data_new 
dim( data_1 )
data_1 <- na.omit( data_1 ) 
dim( data_1 )
bXY = geoXY( latitude = data_1$Lat, longitude = data_1$Lon, 
             lat0 = mean( data_1$Lat ), lon0 = mean( data_1$Lon ), unit = 1000 )
data_1 = cbind( data_1, bXY )
matrix_1 <- rcorr( as.matrix( data_1[,c( "Oxygen", "Bio.Temp", "X", "Y", "Mean.depth" )], type = "pearson" ) )
( matrix_1$r ) 
corrplot( matrix_1$r, type = "lower", tl.col = "black", method = "number", title = "New Zealand arrow squid", 
          mar = c( 0, 0, 1, 0 ), p.mat = matrix_1$P, sig.level = 0.05) #### "X" indicates that the correlation is not significant

######## Check the distribution of max length data 
species_length = data_1$Max.length
summary( species_length )
hist( species_length ) #### It would be reasonable to fit a Gaussian GAM for the species
species_g = fitdist( species_length, distr = "gamma", method = "mle", lower = c( 0, 0 ) )
species_n = fitdist( species_length, distr = "norm", method = "mle", lower = c( 0, 0 ) )
species_l = fitdist( species_length, distr = "lnorm", method = "mle", lower = c( 0, 0 ) )
par( mfrow = c( 2, 2 ) )
plot.legend <- c( "Gamma distribution", "Lognormal distribution", "Normal distribution" )
denscomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
cdfcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
qqcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
ppcomp( list( species_g, species_n, species_l), legendtext = plot.legend )
gofstat( list( species_g, species_n, species_l ) ) 

# converting oxygen
data_1 = data_1 %>% mutate(oxygen_mgL = Oxygen*31.9988/1000)

#GAM
gam_max_length.nos.1 = gam( Max.length ~ te( X, Y ) + s( Bio.Temp, k = 4, bs = "ts" ) + s( oxygen_mgL, k = 4, bs = "ts" ), 
                          Gamma( link = "log" ), method = 'REML', data = data_1 )
summary( gam_max_length.nos.1 ) # R-sq.(adj) =  0.2   Deviance explained = 20.1%
#draw ( gam_max_length.nos.1 )

#Re-run GAM, withour oxygen_mgL (not significant)
gam_max_length.nos = gam( Max.length ~ te( X, Y ) + s( Bio.Temp, k = 4, bs = "ts" ), 
                          Gamma( link = "log" ), method = 'REML', data = data_1 )
summary( gam_max_length.nos ) # R-sq.(adj) =  0.2   Deviance explained = 20.1%

# AIC
logLik.gam(gam_max_length.nos) # -8856.571 (df=20.6297)
logLik.gam(gam_max_length.nos.1) # -8856.527 (df=21.17027)

arrow = data_1

######## Estimate the indices of relative importance of the predictors
#### Estimate the indices of relative importance of the predictors for the GAM fitted to max length data
model_list = list( gam_max_length.nos )

source( "AssessRelativeImportancePredictors.R" )
RIeval_mean_length = relImportantVar( model_list[[1]], data_1, nRep = 10 )
RIeval_mean_length_2 = RIeval_mean_length[,c( 1 : 2, 4 )]
source( "CustomizedBarPlot.R" )
RIeval_mean_length_2 
rownames( RIeval_mean_length_2 ) = c( "Temperature", "Eastings", "Northings" )
RIeval_mean_length_2 <- as.data.frame( RIeval_mean_length_2 )
RIeval_mean_length_2 <- RIeval_mean_length_2[order( - RIeval_mean_length_2$meanInd ),]

tiff( file = "NOS_max_length_GAM_IndicesOfRelativeImportance.tif", width = 3.5, height = 2.5, units = "in", 
      res = 600, compression = "lzw" ) 
par( cex = 0.5 )
custBarChart( RIeval_mean_length_2, confGiven = TRUE, col.bar = 'gray', col.se = 1, ylab = '', xlab = '' )
mtext( side = 1, text = "Predictor", cex = 0.7, line = 1.5, padj = 0.5 )
mtext( side = 2, text = "Index of relative importance", cex = 0.7, line = 1.5, padj = -1 )
dev.off()

save( "RIeval_mean_length", "RIeval_mean_length_2", file = "NOS_max_length_GAM_IndicesOfRelativeImportance.RData" )

######## Evaluate GAMs using Leave Group Out Cross Validation
#### Evaluate the GAM fitted to max length data
model_list = list( gam_max_length.nos )

source( "ModelEvaluation.R" )
Evaluation_GAM_max_length_NOS = evaluateModel( model_list, nIter = 10, confInt = TRUE ) 

save( "Evaluation_GAM_max_length_NOS", file = "Evaluation_GAM_max_length_NOS.RData" )

####################################################

#Scampi
length_data = read.csv( "SCI_trawl_lengths_DATES.csv") #Carapace length

dim( length_data )
nrow( length_data ) #### 29430 lines
names( length_data )
table( length_data$trip_code ) 
summary( length_data$station_no ) #### 1 NA
summary( length_data$dlon_s ) #### 31 NA
summary( length_data$dlat_s ) #### 31 NA
summary( length_data$dlon_e ) #### 134 NA
summary( length_data$dlat_s ) #### 31 NA
summary( length_data$min_gdepth ) #### 4725 NA's
summary( length_data$max_gdepth ) #### 4736 NA's
summary( length_data$bot_temp ) #### 9781 NA's
table( length_data$species ) 
summary( length_data$lgth ) #### 1 NA
summary( length_data$no_a ) #### 1 NA
summary( length_data$no_m ) #### 416 NA
summary( length_data$no_f ) #### 904 NA

length_data_1 = drop_na( length_data )
nrow( length_data_1 ) #### 16564 lines remaining
length_data_2 = length_data[!is.na( length_data$lgth ) & !is.na( length_data$dlon_s ) &
                              !is.na( length_data$dlat_s ) & !is.na( length_data$bot_temp ) &
                              !is.na( length_data$min_gdepth ) & !is.na( length_data$max_gdepth ),]
nrow( length_data_2 ) #### 16983 lines remaining
length_data <- length_data_2 

plot( length_data$dlon_s, length_data$dlat_s )

#Trip station
summary( length_data$station_no ) 
length( unique( length_data$station_no ) )
table( length_data$station_no ) 
trip_station <- paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) )
length( unique( trip_station ) )
length_data$trip_station <- as.factor( paste0( as.factor( length_data$trip_code ), as.factor( length_data$station_no ) ) )
dim( length_data )
names( length_data )

#Mean depth of station
length_data = length_data %>%
  add_column(mean_depth = ((length_data$min_gdepth + length_data$max_gdepth)/2))

#Mean length at station
mean_length_data <- aggregate( lgth ~ trip_station, data = length_data, mean )
mean_length_data$dlon_s <- length_data$dlon_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$dlat_s <- length_data$dlat_s[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$bot_temp <- length_data$bot_temp[match( mean_length_data$trip_station, length_data$trip_station )]
mean_length_data$mean_depth <- length_data$mean_depth[match( mean_length_data$trip_station, length_data$trip_station)]
mean_length_data$date_s <- length_data$date_s[match( mean_length_data$trip_station,length_data$trip_station)]
dim( mean_length_data )
names( mean_length_data )
thedata = mean_length_data
plot( thedata$dlon_s, thedata$dlat_s )

#adding 95th quantile value of length for each Trip station
length_data = as_tibble(length_data)

quant_length_data = length_data %>%
  dplyr::select(date_s, trip_station, lgth, dlon_s, dlat_s, bot_temp)

quant_data <-  data.frame(quant_length_data$date_s, quant_length_data$trip_station, quant_length_data$bot_temp, quant_length_data$dlon_s, quant_length_data$dlat_s, quant_length_data$lgth )
names( quant_data ) <- c("Date", "Trip.Station",  "Situ.Temp", "Lon", "Lat", "Length" ) 

bXY = geoXY( latitude = quant_data$Lat, longitude = quant_data$Lon, 
             lat0 = mean( quant_data$Lat ), lon0 = mean( quant_data$Lon ), unit = 1000 )
quant_data = cbind( quant_data, bXY )

sort( unique( quant_data$Trip.Station ) )

#set stations as factor
quant_data$Trip.Station <- as.factor( paste0( as.factor( quant_data$Trip.Station )))

#here are the 95 quantile values for each trip station
quant_dat = quant_data %>%
  group_by(Trip.Station) %>%
  summarise(enframe(quantile(Length, c(0.95)), "quantile", "95.Length"))

quant_95=as_tibble(quant_dat$`95.Length`)

######## Load in environmental data (mean temperature, and oxygen at maximum depth)
######## Those data were extracted from https://www.bio-oracle.org/downloads-to-email.php

ox = read.asciigrid( "oxygen.asc" )
ox = raster( ox )
x1 <- crop( ox, extent( -180, 0, -90, 90 ) )
x2 <- crop( ox, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
ox <- merge( x1, x2 )
plot( ox )
temp = read.asciigrid( "temp.asc" )
temp = raster( temp )
x1 <- crop( temp, extent( -180, 0, -90, 90 ) )
x2 <- crop( temp, extent( 0, 180, -90, 90 ) )   
extent( x1 ) <- c( 180, 360, -90, 90 )
temp <- merge( x1, x2 )
plot( temp )

predictors <- stack( temp, ox )

######## Combine environmental data with monitoring data and check for outliers in the response and covariate variables
#### Combine environmental data with monitoring data
data_new <-  data.frame( raster::extract( predictors, thedata[3:4] ),
                         thedata$date_s, thedata$bot_temp, thedata$mean_depth, thedata$dlon_s, thedata$dlat_s, thedata$lgth )

#adding 95 quantiles
data_new = cbind(data_new, quant_95)
names( data_new ) <- c( "Bio.Temp", "Oxygen", "Date", "Situ.Temp", "Mean.depth", "Lon", "Lat", "Length", "Max.length" ) 
nrow( data_new ) ## 1407 lines 


#### Check for outliers in the response variable
summary( data_new$Max.length )
hist( data_new$Max.length )
# Two size classes -> take all measurements > 20
data_new = data_new %>% filter(Max.length>20)
hist(data_new$Max.length)
boxplot( data_new$Max.length )
sort( unique( data_new$Max.length ) )

#### Check for outliers in the bottom temperatures
summary( data_new$Bio.Temp )
hist( data_new$Bio.Temp )
boxplot( data_new$Bio.Temp )
sort( unique( data_new$Bio.Temp ) )

#### Check for outliers in the bottom dissolved oxygen (DO) concentrations
summary( data_new$Oxygen )
hist( data_new$Oxygen )
boxplot( data_new$Oxygen )
sort( unique( data_new$Oxygen ) )

#### Check for outliers in the mean depth
summary( data_new$Mean.depth )
hist( data_new$Mean.depth )
boxplot( data_new$Mean.depth )
sort( unique( data_new$Mean.depth ) )

## The diagnostics suggest that the minimum and max (x2) of temperature,
# plus min oxygen (x2) and max depth are outliers, and should, therefore, be removed
nrow( data_new ) ## 564
data_new <- data_new[!data_new$Bio.Temp == min( data_new$Bio.Temp ), ]
data_new <- data_new[!data_new$Bio.Temp == max( data_new$Bio.Temp ), ]
data_new <- data_new[!data_new$Bio.Temp == max( data_new$Bio.Temp ), ]
nrow( data_new ) ## 561
data_new <- data_new[!data_new$Oxygen == min( data_new$Oxygen ), ]
data_new <- data_new[!data_new$Oxygen == min( data_new$Oxygen ), ] 
data_new <- data_new[!data_new$Oxygen == max( data_new$Oxygen ), ]
data_new <- data_new[!data_new$Oxygen == max( data_new$Oxygen ), ] 
nrow( data_new ) ## 556
data_new <- data_new[!data_new$Mean.depth == max( data_new$Mean.depth ), ] 
nrow( data_new ) ## 555

######## Check collinearity 
data_1  = data_new 
dim( data_1 )
data_1 <- na.omit( data_1 ) 
dim( data_1 )
bXY = geoXY( latitude = data_1$Lat, longitude = data_1$Lon, 
             lat0 = mean( data_1$Lat ), lon0 = mean( data_1$Lon ), unit = 1000 )
data_1 = cbind( data_1, bXY )
matrix_1 <- rcorr( as.matrix( data_1[,c( "Oxygen", "Bio.Temp", "X", "Y", "Mean.depth" )], type = "pearson" ) )
corrplot( matrix_1$r, type = "lower", tl.col = "black", method = "number", title = "New Zealand scampi", 
          mar = c( 0, 0, 1, 0 ), p.mat = matrix_1$P, sig.level = 0.05) #### "X" indicates that the correlation is not significant

#Converting length data from mm to cm
data_1 = data_1 %>% mutate(length_cm = Length*0.1)
data_1 = data_1 %>% mutate(max.length_cm = Max.length*0.1)


######## Check the distribution of max length data 
species_length = data_1$Max.length
summary( species_length )
hist( species_length ) #### It would be reasonable to fit a Gaussian GAM for the species
species_g = fitdist( species_length, distr = "gamma", method = "mle", lower = c( 0, 0 ) )
species_n = fitdist( species_length, distr = "norm", method = "mle", lower = c( 0, 0 ) )
species_l = fitdist( species_length, distr = "lnorm", method = "mle", lower = c( 0, 0 ) )
par( mfrow = c( 2, 2 ) )
plot.legend <- c( "Gamma distribution", "Lognormal distribution", "Normal distribution" )
denscomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
cdfcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
qqcomp( list( species_g, species_n, species_l ), legendtext = plot.legend )
ppcomp( list( species_g, species_n, species_l), legendtext = plot.legend )
gofstat( list( species_g, species_n, species_l ) ) 


#Converting oxygen...
data_1 = data_1 %>% mutate(oxygen_mgL = Oxygen*31.9988/1000)

#GAM 
gam_max_length.sci.1 = gam( max.length_cm ~ te( X,Y ) + s( Bio.Temp, k = 4, bs = "ts" ) + s( oxygen_mgL, k = 4, bs = "ts" ), 
                          Gamma( link = "log" ), method = 'REML', data = data_1 )
summary( gam_max_length.sci.1 ) #R-sq.(adj) =  0.276   Deviance explained = 28.1%


#Re-run GAM, removing oxygen_mgL (non-significant & collinear with Y)
gam_max_length.sci.2 = gam( max.length_cm ~ te( X,Y ) + s( Bio.Temp, k = 4, bs = "ts" ), 
                          Gamma( link = "log" ), method = 'REML', data = data_1 )
summary( gam_max_length.sci.2 ) #R-sq.(adj) =  0.277   Deviance explained = 28.1%

#Re-run GAM, removing Y (collinear with oxygen)
gam_max_length.sci = gam( max.length_cm ~ s( X ) + s( Bio.Temp, k = 4, bs = "ts" ) + s( oxygen_mgL, k = 4, bs = "ts" ), 
                            Gamma( link = "log" ), method = 'REML', data = data_1 )
summary( gam_max_length.sci ) #R-sq.(adj) =  0.253   Deviance explained = 25.1%

logLik.gam(gam_max_length.sci) # -324.1929 (df=12.3302)
logLik.gam(gam_max_length.sci.1) # -312.9506 (df=17.1447)
logLik.gam(gam_max_length.sci.2) # -312.7011 (df=17.04145)

scampi = data_1

######## Estimate the indices of relative importance of the predictors
#### Estimate the indices of relative importance of the predictors for the GAM fitted to max length data
model_list = list( gam_max_length.sci )

source( "AssessRelativeImportancePredictors.R" )
RIeval_mean_length = relImportantVar( model_list[[1]], data_1, nRep = 10 )
RIeval_mean_length_2 = RIeval_mean_length[,c( 1 : 2, 4 )]
source( "CustomizedBarPlot.R" )
RIeval_mean_length_2 
rownames( RIeval_mean_length_2 ) = c( "Temperature", "Eastings", "Oxygen" )
RIeval_mean_length_2 <- as.data.frame( RIeval_mean_length_2 )
RIeval_mean_length_2 <- RIeval_mean_length_2[order( - RIeval_mean_length_2$meanInd ),]

tiff( file = "SCI_max_length_GAM_IndicesOfRelativeImportance.tif", width = 3.5, height = 2.5, units = "in", 
      res = 600, compression = "lzw" ) 
par( cex = 0.5 )
custBarChart( RIeval_mean_length_2, confGiven = TRUE, col.bar = 'gray', col.se = 1, ylab = '', xlab = '' )
mtext( side = 1, text = "Predictor", cex = 0.7, line = 1.5, padj = 0.5 )
mtext( side = 2, text = "Index of relative importance", cex = 0.7, line = 1.5, padj = -1 )
dev.off()

save( "RIeval_mean_length", "RIeval_mean_length_2", file = "SCI_max_length_GAM_IndicesOfRelativeImportance.RData" )

######## Evaluate GAMs using Leave Group Out Cross Validation
#### Evaluate the GAM fitted to max length data
model_list = list( gam_max_length.sci )

source( "ModelEvaluation.R" )
Evaluation_GAM_max_length_SCI = evaluateModel( model_list, nIter = 10, confInt = TRUE ) 

save( "Evaluation_GAM_max_length_SCI", file = "Evaluation_GAM_max_length_SCI.RData" )

all_spp_filtered_dates = rbind(hoki, snapper, sbw, roughy, arrow, scampi )

