##1. Descriptive analysis of the dependent and independent variables
library(sp)
library(spdep)
library(ggplot2)
library(dplyr)
library(sf)
library(readr)
library(car)
np2<-read.csv("COVID_Model2_prev.csv",header=TRUE, as.is=TRUE)
colnames(np2)<-c("DISTRICT","COVID_cases2", "Tot_pop21","COVID_Inc2",  
                "Sex_ratio21","Pop_den21",
                "DM_1718","BrA_1920",
                "COPD_1920","Health_fac_1718","avg_temp1822",
                "avg_precip1822","p_wo_safe_dr_w", "HH_crowding",
                "p_unemployed_F","p_unemployed_M","Smoking_p_men","Smoking_rate_female",
                "p_wo_hw_fac", 
                 "Acc_to_healthfac","Obesity_M","Obesity_F",
                "p_graduate","Gini_coeff","Median_age", "HTN_1920")

head(np2)
np2[is.na(np2)]<-0 #replace NAs as 0


#summary statistics
np2$DISTRICT<-as.factor(np2$DISTRICT)
summary(np2)
library(dplyr)
np2 <- np2 %>% mutate_at(c('COVID_Inc2',  
                'Sex_ratio21','Pop_den21',
                'DM_1718','BrA_1920',
                'COPD_1920','Health_fac_1718','avg_temp1822',
                'avg_precip1822','p_wo_safe_dr_w','HH_crowding',
                'p_unemployed_F','p_unemployed_M','Smoking_p_men','Smoking_rate_female',
                'p_wo_hw_fac', 
                'Acc_to_healthfac','Obesity_M','Obesity_F',
                'p_graduate','Gini_coeff','Median_age','HTN_1920'), as.integer)

subset.data <- subset(np2, select = c("COVID_Inc2",  
                "Sex_ratio21","Pop_den21",
                "DM_1718","BrA_1920",
                "COPD_1920","Health_fac_1718","avg_temp1822",
                "avg_precip1822","p_wo_safe_dr_w", "HH_crowding",
                "p_unemployed_F","p_unemployed_M","Smoking_p_men","Smoking_rate_female",
                "p_wo_hw_fac"
                ,"Acc_to_healthfac","Obesity_M","Obesity_F","p_graduate","Gini_coeff","Median_age","HTN_1920"))
str(subset.data)

#Standard deviation calculation
sapply(subset.data, sd)
colMeans(subset.data, na.rm=TRUE)
summary(subset.data)

sd(np2$COVID_cases2)
sd(np2$Tot_pop21)


## 2.Check for Overdispersion
require(qcc)
qcc.overdispersion.test(np2[!is.na(np2$COVID_Inc2),]$COVID_Inc2, type="poisson")
hist(np2$COVID_Inc2)

## 3.Negative Binomial Regression: 
### Univariable linear regression: Demographic variable

require(MASS)
names(np2)

#Standardization of the variables
np2[,c(5:26)] <- sapply(np2[,c(5:26)],scale)

covar1<-np2[,c(7:24,26)]
glm_result1=NULL
lmrs5= NULL
lm_coeff5= NULL
lm_pval5= NULL
pred1=NULL
mydata1=NULL

for(n in seq_along(covar1)){
  lmrs5[[n]] = summary(glm.nb(COVID_cases2~ offset(log(Tot_pop21)) +  covar1[[n]] , data=np2)) 
  lm_coeff5[[n]]<-lmrs5[[n]]$coefficients[[2]]
  lm_pval5[[n]]<-lmrs5[[n]]$coefficients[[2,4]]
}


lm_result7<-cbind(unlist(lm_coeff5), unlist(lm_pval5))
rownames(lm_result7)<-names(covar1)
colnames(lm_result7)<- c("coefficient","pval")

lm_result7
print(lm_result7[lm_pval5 < 0.25,])


### Multiple negative binomial (Model adjusted for covariates population density, sex ratio, climate variables and mean age)


require(MASS)
require(foreign)
require(ggplot2)

summary(v1<- glm.nb(COVID_cases2 ~ HH_crowding + p_graduate + p_unemployed_F + Gini_coeff +
                      p_unemployed_M  + avg_temp1822 + avg_precip1822 +
                       Acc_to_healthfac  + Health_fac_1718 + p_wo_hw_fac +  
                      BrA_1920 + COPD_1920 + DM_1718 + HTN_1920 + Obesity_M + Obesity_F +
                       Sex_ratio21 +  Median_age +
                      Pop_den21 + offset(log(Tot_pop21)),data=np2))


## VIF for NB

require(car)
 car::vif(v1)

require(MASS)
require(foreign)
require(ggplot2)


summary(v3<- glm.nb(COVID_cases2 ~ HH_crowding + p_unemployed_F + Gini_coeff +
                      p_unemployed_M  + avg_temp1822 + avg_precip1822 +
                       Acc_to_healthfac  + Health_fac_1718 + p_wo_hw_fac +  
                      BrA_1920 + COPD_1920 + DM_1718 + HTN_1920 + Obesity_M + Obesity_F +
                       Sex_ratio21 +  Median_age +
                      Pop_den21 + offset(log(Tot_pop21)),data=np2))

## VIF for NB
require(car)
 car::vif(v3)

cbind(coefest= coef(v3), confint(v3))


# Model selection

#backward elimination based on AIC.
v2=step(v3)
summary (v2)


#Model diagnostics using the DHARMa package.
The function plot() in DHARMa plots residuals against a predictor (by default against the fitted value, extracted from the DHARMa object, or any other predictor).

Outliers are highlighted in red (for information on definition and interpretation of outliers, see testOutliers).

To provide a vinp2al aid in detecting deviations from uniformity in y-direction, the plot function calculates an (optional) quantile regression of the residuals, by default for the 0.25, 0.5 and 0.75 quantiles. As the residuals should be uniformly distributed for a correctly specified model, the theoretical expectations for these regressions are straight lines at 0.25, 0.5 and 0.75, which are displayed as dashed black lines on the plot. Some deviations from these expectations are to be expected by chance, however, even for a perfect model, especially if the sample size is small. The function therefore tests if deviation of the fitted quantile regression from the expectation is significant, using testQuantiles. If so, the significant quantile regression will be highlighted as red, and a warning will be displayed in the plot.

#DHARMa is a very good package for assessing diagnostic plots for GLMs and GLMMs.
#There is a package called DHARMa, which implements simulations-based model-checking. The basic steps of the DHARMa workflow are:

#1.Simulate new data from the fitted model
#2.For each observation, calculate the cumulative density function (CDF) of the simulated data
#3.Calculate the residuals as the value of the CDF at the value of the observed data
#Assuming the model is correctly specified, all values of the CDF (ranging from 0 to 1) should appear with equal probability. The residuals should thus be uniformly distributed, which we can check by comparing the simulated residuals versus a uniform distribution using QQ-plots.


library(DHARMa)
negbin_pnu_simulation <-simulateResiduals (fittedModel= v2 )
plot(negbin_pnu_simulation)


# IN a DHARMa QQ plot, the residuals are expected to follow a uniform distribution instead of a normal distribution, and are standardized to values between 0 and 1. The QQ plot above shows that there is no over or under dispersion or zero inflation in the model. 

##4. Bayesian spatial regression analysis of COVID-19 incidence
### Neighborhood structure and Spatial Weight Matrix

np3<-rgdal::readOGR("hermes_npL_new_wgs_2.shp", encoding = 'utf-8')
np2.sf<- st_as_sf(np3)

# reprojection using european petroleum np2rvey group

np3<-spTransform(np3,"+init=epsg:4326")
require(spdep)
np3@data$utmX<-coordinates(np3)[,1]
np3@data$utmy<-coordinates(np3)[,2]
np3@data$utmX

# base map plotting 
# Add fill layer to UZB shape
library(tmap)
tm_shape(np3)+tmap_options(check.and.fix = TRUE) +  tm_fill() 
# Add border layer to UZB shape
tm_shape(np3) + tm_borders() 
# Add fill and border layers to UZB shape
tm_shape(np3) + tm_fill() + tm_borders() 

### Merging shapefile and data
library(dplyr)

# Key variable: DISTRICT
np3@data$DISTRICT <- as.factor(np3@data$DISTRICT)
np2$DISTRICT<- as.factor(np2$DISTRICT)
np3@data<-dplyr::left_join(np3@data, #left side data
                           np2, #the data you want to add
                           by=c("DISTRICT"="DISTRICT")) # key factor
head(np3@data) # check whether data is well merged
np3@data

# define weights matrix
library(maptools)
library(spdep)
library(sf)
library(spatialreg)
install.packages(spdep)
library(spdep)
library(spdep)
pdf("Neighborhood weights.pdf")
coord<-coordinates(np3)
pp<-poly2nb(np3, queen= FALSE)

# define weights matrix wt
wt <- nb2listw(pp, style = "W", zero.policy = T)
summary(wt)
xaa=coordinates (np3)[,1]
yaa=coordinates (np3)[,2]
coords<-cbind(xaa,yaa)

plot(np3, border="white",col="cadetblue")
plot(wt,coord, col="blue",cex=0.4,add=TRUE)
title(main="Neighborhood spatial weights")
box(col="black")
dev.off()

### Moran's I statistics for test of spatial autocorrelation

cor_sp1 = moran.mc(np3@data$COVID_Inc2 ,listw = wt, nsim = 999, na.action=na.exclude)
cor_sp1$statistic 
cor_sp1$p.value 


### Exploratory spatial data analysis

# COVID Incidence
library(maptools)

map_np3<- tm_shape(np3) + 
  tm_polygons(col = "COVID_Inc2", id="DISTRICT",n=6,
              palette = "Purples", legend.title= "Total COVID incidence", title="Incidence")+
  tm_legend(position = c("left", "bottom"), frame = TRUE)+
  tm_layout(frame = TRUE,main.title="Spatial distribution of COVID Incidence", main.title.position = "center", main.title.size=1)

map_np3+ tm_compass(type = "8star", position = c("right", "top")) + tm_scale_bar(breaks = c(0, 100, 200),  text.size = 0.3, position="center")



local1 <- localmoran(x = np3@data$COVID_Inc2, listw = nb2listw(pp, style = "W"))
local1

# Plot Local Moran

# binds renp2lts to our polygon shapefile

moran.map <- cbind(np3, local1)

tm_shape(moran.map, ) +
  tm_fill(col = "Ii",
          style = "quantile",
          title = "local moran statistic") +
  tm_layout(frame = TRUE,main.title="Local moran's I plot for COVID Incidence", main.title.position = "center", main.title.size=1)

### Getis ord Gi* clustering: 
#Install and load packages

packages = c('rgdal', 'spdep',  'tmap', 'tidyverse')
for (p in packages){
if(!require(p, character.only = T)){
install.packages(p)
}
library(p,character.only = T)
}

#Local G function

local_Gi <- spdep::localG(x=np3@data$COVID_Inc2, listw=wt)
local_Gi

#Calculate p-value for the local G statistics

pvals= pnorm(q=2*abs(local_Gi), lower.tail= FALSE) # for absolute value of z-statistics
pvals

#Plot significance results

options("install.lock"=FALSE)
install.packages("wesanderson", dependencies = TRUE, INSTALL_opts = "--no-lock")
pdf("Getis ord Gi clustering.pdf")
np3@data$local_Gi <- local_Gi
np3@data$pvals <- pvals
np3_sig <- np3[np3@data$pvals<=0.01,]
tmap=tm_shape(shp=np3) + tm_borders(alpha=1) + tm_shape(shp=np3_sig) + tm_borders(alpha=1) + tm_fill(col="local_Gi",palette = "-RdYlBu", title="local Gi", breaks=c(-2.58,-1.96,-1.65,1.65,1.96,2.58,8)) + tm_borders(alpha=1)
tmap
dev.off()


### 90% significant: Gi* z score > 1.645; 95% significant: Gi* z score > 1.960;
##99% significant: Gi* z score > 2.576; 99.9% significant: Gi* z score > 3.291


### Hierachial spatial analysis: 

install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
install.packages(c("sp", "spdep", "raster", "rgdal", "rgeos", "ggplot2", "leaflet", "DT", "dplyr", "SpatialEpi", "geoR"))
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(spdep)

names(inla.models()$likelihood)
inla.doc("familyname")
names(inla.models()$prior)
inla.doc("priorname")

library(maptools)
library(spdep)

#Convert Neighbor map to an INLA integrated map that INLA can use to fit a model.
np2_nb <- poly2nb(np3,queen=FALSE)
nb2INLA("np2.graph", np2_nb)
np2.adj <- paste(getwd(),"/np2.graph",sep="")
H <- inla.read.graph(filename="np2.graph")
image(inla.graph2matrix(H), xlab="", ylab="")
np2.adj

#Above is the adjency matrix for the Nepal COVID incidence. Rows and columns idenitfy areas; squares identify neighbors.
plot(np3, border="white",col="cadetblue")
plot(np2_nb,coord, col="blue",cex=0.4,add=TRUE)
title(main="Neighborhood spatial weights")

## Setting the prior
In R_INLA the default prior for the precision of the random effects (ui) is 1/σ^2 ~ gamma(1,0.00005)
We can change this prior by setting a Penalized complexity (PC) prior on the standard deviation σ.
Penalized complexity prior helps to avoid over-fitting by shrinking the marginal variance to zero.
We can specify the probability of σ being greater than 0.5 is small, equal to 0.01: P(σ>0.5) =0.01 (This is the default value for the PC prior in RINLA)


inla.doc("pc.prec")
#Hyperparameter with penalized complexity prior
hyperbym<-list(prec.unstruct=list(prior="pc.prec", param=c(0.5,0.01)), prec.spatial=list(prior="pc.prec", param=c(0.5,0.01)))

##We develop a id area to identify each district.
np3@data$idarea <- 1:nrow(np3@data)
np3@data

## Write formula.

#GLM
f0<-COVID_cases2 ~ HH_crowding + Gini_coeff + p_unemployed_M + 
    avg_temp1822 + Acc_to_healthfac + DM_1718 + HTN_1920 + Obesity_F + 
     Sex_ratio21 + Median_age + Pop_den21 + offset(log(Tot_pop21))
     
#BESAG ICAR
f<-COVID_cases2 ~ HH_crowding + Gini_coeff + p_unemployed_M + 
    avg_temp1822 + Acc_to_healthfac + HTN_1920 + Obesity_F + 
    DM_1718 + Sex_ratio21 + Median_age + Pop_den21 + offset(log(Tot_pop21)) + f(idarea, model="besag",
                          graph=np2.adj, scale.model= TRUE)
                          
#BYM with penalized complexity prior
f1<-COVID_cases2 ~ HH_crowding + Gini_coeff + p_unemployed_M + 
    avg_temp1822 + Acc_to_healthfac + HTN_1920 + Obesity_F + 
    DM_1718 + Sex_ratio21 + Median_age + Pop_den21 + offset(log(Tot_pop21)) + f(idarea, model="bym",
                          graph=np2.adj,scale.model= TRUE,hyper=hyperbym)
                          

np2data<- np3@data
np2data<-np3@data[,-32]

np2data1<- np2data
np2data1<-np2data[,-32]


# Fitting the INLA models and compute posteriors of the prediction. 
#GLM
res0<-inla(f0, data=np2data1, family="nbinomial", control.compute = list(dic=TRUE, waic=TRUE), control.predictor = list(compute=TRUE))
summary(res0)

#Besag
res<-inla(f, data=np2data1,family="nbinomial",control.compute = list(dic=TRUE, waic=TRUE), control.predictor = list(compute=TRUE))
summary(res)

#BYM with PC prior
res1<-inla(f1, data=np2data1,family="nbinomial", control.compute = list(dic=TRUE, waic=TRUE), control.predictor = list(compute=TRUE))
summary(res1)



#Comparing DIC and WAIC of the models
dic  <- c(res0$dic$dic, res$dic$dic, res1$dic$dic)   
waic <- c(res0$waic$waic, res$waic$waic, res1$waic$waic )
Z.out     <- cbind(dic, waic)
rownames(Z.out) <- c("GLM", 
                     "besag ICAR",
                     "bym" )
Z.out


# Calculating variance of the random effect using tmarginal(transforming the marginals)

marg.variance<-inla.tmarginal(function(x) 1/x,
                             res1$marginals.hyperpar$"Precision for idarea (iid component)")

#obtain summary statistics 
inla.zmarginal(marg.variance)

np3$RR<-res1$summary.fitted.values[,"mean"]
fitted.values<- res1$summary.fitted.values[,1]
observed.values<-np3@data$COVID_cases2
rmse <- sqrt(mean((observed.values-fitted.values)^2))
np3$LL<-res1$summary.fitted.values[,"0.025quant"]
np3$UL<-res1$summary.fitted.values[,"0.975quant"]

summary(np3@data[,c("RR","LL","UL")])


###Model validation

#Since the BYM model with pc prior has the least value for DIC and WAIC, INLA model res1 is chosen.
mu2      <- res1$summary.fitted.values[,"mean"]
E2       <- (np3$COVID_cases2 - mu2) / sqrt(mu2)

# plot validation

par(mfrow = c(2,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = mu2, 
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)



rmse2 <- sqrt(((observed.values-fitted.values)^2))
np3$RMSE <- rmse2

##Observed Vs. Predicted by INLA BYM
require(RColorBrewer)
spplot(np3, c("COVID_cases2", "RR"), col.regions=brewer.pal(9, "RdPu"), cuts=8,
       names.attr=c("Observed", "Predicted by INLA-BYM"))

#Predicted Vs. standard error in the predicted.
spplot(np3, c("RR", "RMSE"), col.regions=brewer.pal(9, "RdPu"), cuts=8,
       names.attr=c("Predicted by INLA-BYM", "SE plot"))


#Calculating relative risk

#BYM
# Calculating Relative risk by transforming the marginals.

tmarg <- inla.tmarginal(function(x)exp(x),res1$marginals.fixed$Acc_to_healthfac)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)

# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res1$marginals.fixed$Obesity_F)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)


# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res1$marginals.fixed$DM_1718)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)

# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res1$marginals.fixed$HH_crowding)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)

# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res1$marginals.fixed$Gini_coeff)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)


# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res1$marginals.fixed$p_unemployed_M)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)


# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res1$marginals.fixed$avg_temp1822)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)

# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res1$marginals.fixed$HTN_1920)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)

#Besag ICAR
# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res$marginals.fixed$Acc_to_healthfac)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)

# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res$marginals.fixed$Obesity_F)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)


# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res$marginals.fixed$DM_1718)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)

# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res$marginals.fixed$HH_crowding)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)

# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res$marginals.fixed$Gini_coeff)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)

# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res$marginals.fixed$p_unemployed_M)
inla.qmarginal(c(0.025,0.5,0.975),tmarg) 
inla.emarginal(function(x)x,tmarg)

# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res$marginals.fixed$avg_temp1822)
inla.qmarginal(c(0.025,0.5,0.975),tmarg) 
inla.emarginal(function(x)x,tmarg)

# Calculating Relative risk
tmarg <- inla.tmarginal(function(x)exp(x),res$marginals.fixed$HTN_1920)
inla.qmarginal(c(0.025,0.5,0.975),tmarg)
inla.emarginal(function(x)x,tmarg)
```
