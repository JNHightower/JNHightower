##Ensemble Modeling GWWA Project, 1/3/2011


##set working directory
setwd ('C:/Program Files/R/R2120/library/')
dir()


##Install libraries:

library(BIOMOD)
library(rpart)
library(MASS)
library(gbm)
library(gam)
library(nnet)
library(mda)
library(randomForest)
library(Design)
library (Hmisc)
library(reshape)
library(plyr)
library(boot)
library(maptools)
library(sp)
library(rgdal)


##----------------------------------------------------------------------------------------------------
## link BIOMOD Functions
setwd('C:/Program Files/R/R2120/library/BIOMOD/R/')
source("Models.R")
source("LoadModels.R")
source("pseudo.abs.R")
source ("SampleMat2.R")
source("Biomod.Models.R")
source("scope.R")
source("scopeExpSyst.R")
source("sre.R")
source("somers2.R")
source("multiple.plot.r")
source("Initial.State.R")
source("TSS.Stat.R")
source("KappaStat.R")
source("KappaSRE.R")
source("KappaRepet.R")
source("functionkeep.R")
source("response.plot.R")
source("CurrentPred.R")
source("CVnnet.R")
source("Ensemble.Forecasting.R")
source("Ensemble.Forecasting.raster.R")
source("Dispersal.limit.R")
source("Rescaler4.R")
source("Projection.raster.R")
source("ProjectionBestModel.R")
source("PredictionBestModel.R")
source("Projection.R")
source("ProbDensFunc.R")
source("printcp.R")
source("CutOff.Optimised.R")
source("testnull.R")
source("BIOMOD_internal.R")
source("response.plot.R")
source("LoadModels.R")
source("LoadProj.R")
source("KappaSRE.R")
source("level.plot.R")
source("BinaryTransformation.R")
source("FilteringTransformation.R")
##################################################################################

####################1km Data   ########################################
 

setwd ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/')
dir()

### All Data (1935-2009):
bird1<-read.csv("CoorXY_GW_BW_Current.csv", header=TRUE)
CoorXY<-bird1[,2:3]
head(CoorXY)
nrow(CoorXY)
plot(CoorXY[,1],CoorXY[,2])


bird<-read.csv("GW_BW_Current.csv", header=TRUE)
length(bird)
nrow(bird)
head(bird)
names(bird)
attach(bird)
plot(bird[,2],bird[,3])

## Only Modern Data (1998-2010):
setwd ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/Modern/')
dir()
bird1<-read.csv("CoorXY_GW_BW_Current_1998_2009.csv", header=TRUE)
CoorXY<-bird1[,2:3]
head(CoorXY)
nrow(CoorXY)
plot(CoorXY[,1],CoorXY[,2])


bird<-read.csv("GW_BW_Current_1998_2009.csv", header=TRUE)
length(bird)
nrow(bird)
head(bird)
names(bird)
attach(bird)
plot(bird[,2],bird[,3])

##   Only Historical Data (935-1997):
setwd ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/Historical/')
dir()
bird<-read.csv("GW_BW_Current.csv", header=TRUE)
length(bird)
nrow(bird)
head(bird)
names(bird)
attach(bird)
plot(bird[,2],bird[,3])

## Only Modern Data (1998-2010):
bird1<-read.csv("CoorXY_GW_BW_Current.csv", header=TRUE)
CoorXY<-bird1[,2:3]
head(CoorXY)
nrow(CoorXY)
plot(CoorXY[,1],CoorXY[,2])


###################################################################################################

level.plot(bird[,12], CoorXY[,1:2], show.scale=T,title='GWWA All Data')
level.plot(bird[,13], CoorXY[,1:2], show.scale=T,title='BWWA All Data')

multiple.plot(bird[,18:23], CoorXY[,1:2],cex=0.5)##BIOCLIM VARS.
multiple.plot(bird[,24:26], CoorXY[,1:2],cex=0.5)##PPT VARS.
multiple.plot(bird[,27:29], CoorXY[,1:2],cex=0.5)##TMAX VARS.
multiple.plot(bird[,30:32], CoorXY[,1:2],cex=0.5)##TMIN VARS.
multiple.plot(bird[,33:35], CoorXY[,1:2],cex=0.5)##ELEV, SLOPE, ASPECT

#############  Histograms and data 1km ###################################

Initial.State(sp.name="GW_BW",Response = bird[,12:13], Explanatory = bird[,18:35],
IndependentResponse = bird[,12:13], IndependentExplanatory = bird[,18:35])



ls()
head(DataBIOMOD)
Biomod.material
#####################################################################

Models(GLM = F, TypeGLM = "poly", Test = "AIC",GBM = T, No.trees = 1000, GAM = T,
Spline = 3, CTA =T, CV.tree = 50, ANN = T, CV.ann = 2, SRE = T, quant=0.025, FDA = F,
MARS = T, RF = T, NbRunEval = 3, DataSplit = 80, Yweights=NULL, Roc = T, Optimized.Threshold.Roc = T,
Kappa = T, TSS=T, KeepPredIndependent = F, VarImport=T, 
NbRepPA=3, strategy="random", nb.absences=1000)



CurrentPred(GLM=F, GBM=T, GAM=T, CTA=T, ANN=T, SRE=T, FDA=F, MARS=T,
RF=T,BinRoc=T, BinKappa=T, BinTSS=T, FiltKappa=T)

library(raster)
library(maptools)
library)gdal)
bio1= raster ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/Modern/bio1.asc')
bio5= raster ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/Modern/bio5.asc')
bio10= raster ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/Modern/bio10.asc')
bio15= raster ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/Modern/bio15.asc')

CurrentStack=stack(bio1, bio5, bio10, bio15)
library(raster)




Projection.raster(RasterProj = CurrentStack, Proj.name='CurrentRaster', GLM = F, 
GBM = T, GAM = T, CTA = F, ANN = F, SRE = F, FDA=F, MARS = F, RF = T, BinRoc = F, BinKappa = T, 
BinTSS = T, FiltRoc = T, FiltKappa = T, FiltTSS = T, repetition.models=T)

###############################################################################

##Calcuating the number of pseudo-absences (nb.absence) to run in model:
##The number of presences for each species in model:
length(Biomod.PA.data$GWWA)
length(Biomod.PA.data$BWWA)
##the number of presences for GWWA and BWWA:
sum(bird[,"GWWA"])
sum(bird[,"BWWA"])
##Hence, the number of absences available for calibration:
length(Biomod.PA.data$GWWA) - sum(bird[,"GWWA"])
length(Biomod.PA.data$BWWA) - sum(bird[,"BWWA"])

###################################################################################
## Evaluating model output:

Evaluation.results.Kappa
Evaluation.results.TSS
Evaluation.results.Roc

VarImportance
###################################################################################
##PA data generated
##Biomod.PA.data contains the amount of data available after the inner run of the pseudo-absence
##function. Biomod.PA.sample contains the rows to take from DataBIOMOD to get the data that
##has been used for the calibration of each species for each PA run. The last object is a result of
##the pesudo-absence function inner run and is of no importance here (but see the "Pseudo-absences"
##section of the Presentation manual for explanations).

##For example, let's see what data has been used for the calibration of the run PA1 :

our.lines <- Biomod.PA.sample$GWWA$PA1
par(mfrow=c(1,2))
level.plot(DataBIOMOD[, "GWWA"], CoorXY, title='original data', show.scale=F)
level.plot(DataBIOMOD[our.lines, "GWWA"], CoorXY[our.lines,], title='PA1', show.scale=F)

#####################################################################################
##Each algorithm (excepted SRE) generates an object storing the different parameterisation, the 
## importance of each variable for the model and other statistics. This output is essential as it allows
##generating predictions.
##These objects, the models themselves, are now stored out of the R workspace directly on the comput-
##ers' hard disk. They are named after the algorithm used and the species' names, i.e. Sp164 FDA for
##example. There is also extensions of the names concerning the repetitions and the pseudo-absences
##runs, so that one of our models will be Sp164 FDA PA1 rep2.
##Back loading the models and having them directly usable is very straightforward : simply use
##the load() function to have the model restored in the R workspace, with the same name plus the
##directory root. This is also the case with the other outputs stored outside of R (predictions and
##projections). The syntax is not always handy but easy to pick up :

#Example of the GBM
load("models/GWWA_GBM_PA1")
ls()
summary(GWWA_GBM_PA1)
load("models/BWWA_GBM_PA1")
summary(BWWA_GBM_PA1)


par(mfrow=c(2,2))
plot(BWWA_GBM_PA1)

##The gbm library also provides an experimental diagnostic tool that plots the ##fitted values versus the actual average values. Uses gam to estimate E(yjp). 
##Well-calibrated predictions imply that
##E(yjp) = p. The plot also includes a pointwise 95 band.
##This method can be applied to all models to visualise the relative goodness of ##fit of the model. The function requires the observed presence-absence of the 
##selected species and the predictions. Hence, you will need top load the predictions for this.

library(gbm)
load("pred/pred_BWWA")
##let's store the data that was used for calibration of the
##first PA run for Sp277 to simplify the code

data.used <- DataBIOMOD[Biomod.PA.sample$BWWA$PA1,"BWWA"]
calibrate.plot(data.used, BWWA[,"GLM",1,1]/1000)
dim(Biomod.PA.sample$BWWA$PA1)

load("models/GWWA_GBM_PA1")
plot(GWWA_GBM_PA1, i.var=2)

#####################################################################################################
##There are several useful outputs in CTA models. A critical one is frame which gives the 
##details of the node, the explained deviance by each node (dev) and the probability of 
##occurrences (yval).


##CTA

load("models/GWWA_CTA_PA1")
names(GWWA_CTA_PA1)
GWWA_CTA_PA1$frame

plot(GWWA_CTA_PA1, margin=0.05)
text(GWWA_CTA_PA1, use.n=T)
########################################################################################################
## RF
## The importance of each variable, as produced by random Forest, can be extracted.
load("models/GWWA_RF_PA1")
GWWA_RF_PA1$importance


##Here are the definitions of the variables' importance measures.
##- Mean Decrease Accuracy: For each tree, the prediction accuracy on the 
out-of-bag portion of the
##data is recorded. Then the same is done after permuting each predictor variable. 
##The difference between the two accuracies are then averaged across all trees, 
##and normalized by the standard error.
##- Mean Decrease Gini: The second measure is the total decrease in node
## impurities from splitting on
##the variable, averaged over all trees. For classication, 
##the node impurity is measured by the Gini index.

varImpPlot(GWWA_RF_PA1)

##########################################################################################
## Response Plots


##this one has already been loaded in a prior call
response.plot(GWWA_GBM_PA1, bird[18:35])


load("pred/Pred_GWWA")
names(GWWA_GAM_PA1)
library(gbm)
calibrate.plot(data.used,Pred_GWWA_indpdt[,"GAM",1,1]/1000)

level.plot(Pred_GWWA[,"RF",1,1], bird1[Biomod.PA.sample[[1]] [[1]],] , title='GWWA RF')

response.plot(GWWA_GAM_full,bird[24:35])

load("models/GWWA_GAM_PA1")
response.plot(GWWA_GAM_PA1, bird[24:35])

load("models/GWWA_RF_PA1")
response.plot(GWWA_RF_PA1, bird[24:35])

load("models/BWWA_GAM_PA1")
response.plot(BWWA_GAM_PA1,bird[24:35])

load("models/BWWA_RF_PA1")
response.plot(BWWA_RF_PA1,bird[24:35])

load("pred/Pred_GWWA_indpdt")
Pred_GWWA_indpdt[1500:2500,"GBM",,]
level.plot(Pred_GWWA_indpdt[,"GAM",1,1], CoorXY, title="GWWA")

load("pred/Pred_BWWA_indpdt")
Pred_BWWA_indpdt[1:500,"GBM",,]
level.plot(Pred_BWWA_indpdt[,"GAM",1,1], CoorXY, title="BWWA")

map.plot(Sp='GWWA',model='all',method='Kappa',format.type='probs',wanted='prediction')

################################################################################
## Assessing the best model for modern/ historical prediction:

CurrentPred(GLM=F, GBM=T, GAM=T, CTA=T, ANN=T, SRE=T, FDA=F, MARS=T,
RF=T,BinRoc=T, BinKappa=T, BinTSS=T, FiltRoc=T, FiltKappa=T, FiltTSS=T)

PredictionBestModel(GLM=F,GBM=T, GAM=T, CTA=T, ANN=T, SRE=T, FDA=F, MARS=T, RF=T,
method='all', Bin.trans = F, Filt.trans = T)


load("pred/BestModelByRoc")
BestModelByRoc
load("pred/BestModelByTSS")
BestModelByTSS
load("pred/BestModelByKappa")
BestModelByKappa


load("pred/Pred_GWWA_FiltKappa")

par(mfrow=c(1,3))
level.plot(Pred_GWWA_BinKappa[,"GAM",2,1], CoorXY, show.scale=F, title='probabilities')
level.plot(Pred_GWWA_BinKappa[,"GAM",2,1], CoorXY, show.scale=F, title='binary data')
level.plot(Pred_GWWA_FiltKappa[,"GAM",2,1], CoorXY, show.scale=F, title='filtered data')

load("pred/Pred_GWWA")

GWWA = cbind(CoorXY[Biomod.PA.sample[[1]][[1]],], 
Pred_GWWA[,"GAM",1,1])
write.csv(GWWA, "GWWA Prediction.csv") 

GWWA = cbind(CoorXY[Biomod.PA.sample[[1]][[2]],], 
Pred_GWWA[,"RF",1,1])
write.csv(GWWA, "GWWA Prediction_RF.csv") 

GWWA = cbind(CoorXY[Biomod.PA.sample[[1]][[2]],], 
Pred_GWWA[,"RF",1,1])
write.csv(GWWA, "GWWA Prediction_RF.csv") 

load("pred/Pred_BWWA")
load("pred/Pred_BWWA_FiltKappa")

dimnames(Pred_BWWA)
Biomod.PA.sample$BWWA$PA1
Biomod.PA.sample$GWWA$PA1

BWWA = cbind(CoorXY[Biomod.PA.sample[[2]][[1]],], 
Pred_BWWA[,"RF",1,1])
write.csv(BWWA, "BWWA Predictions_RF.csv") 

BWWA = cbind(CoorXY[Biomod.PA.sample[[2]][[1]],], 
Pred_BWWA[,"MARS",1,1])
write.csv(BWWA, "BWWA Predictions_MARS.csv") 

BWWA = cbind(CoorXY[Biomod.PA.sample[[2]][[1]],], 
Pred_BWWA[,"GBM",1,1])
write.csv(BWWA, "BWWA Predictions_GBM.csv") 

load("pred/PredBestModelBYRoc$BWWA")

#################################################################################

Projection(Proj = bird[,18:35], Proj.name='Current',
GLM=F, GBM=T, GAM=T, CTA=T, ANN=T, SRE=T, FDA=F, MARS=T,
RF=T, BinRoc=T, BinKappa=T, BinTSS=T, FiltRoc=T, 
FiltKappa=T, FiltTSS=T, repetition.models=T)


Ensemble.Forecasting(Proj.name= "Current", weight.method='Roc', PCA.median=T,
binary=T, bin.method='Roc', Test=F, decay=1.6, repetition.models=T)

load("proj.Current/Total_consensus_Current")

dim(Total_consensus_Current)
dim(Total_consensus_Current)
dimnames(Total_consensus_Current[,1,2])

GWWA = cbind(CoorXY[Biomod.PA.sample[[1]][[0]],], Total_consensus_Current[,1,1])
write.csv(GWWA, "GWWA Prediction.csv") 
GWWA = cbind(CoorXY,Pred_GWWA[,"GAM",2,1])
write.csv(GWWA, "GWWA Prediction.csv") 

BWWA = cbind(CoorXY[Biomod.PA.sample[[1]][[0]],], consensus_BWWA_Current[,"GAM",2,1]), 
write.csv(BWWA, "GWWA Prediction.csv") 
BWWA = cbind(CoorXY,Pred_GWWA[,"GAM",2,1])
write.csv(BWWA, "GWWA Prediction.csv") 



## Load future climate scenrio (A2a, A1BG, B2a--variable names have to be the same in all 
## input files##
setwd ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/')
A2a<-read.csv("GW_BW_A2a.csv", header=TRUE)
names(A2a)
attach(A2a)
setwd ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/models/')
Projection(Proj = A2a[,24:35], Proj.name='A2a', GLM = FALSE, GBM = TRUE, GAM = TRUE,
CTA = TRUE, RF = TRUE,BinRoc = TRUE, BinKappa = TRUE, BinTSS = TRUE, FiltRoc = TRUE, FiltKappa = TRUE, FiltTSS = TRUE,
repetition.models=TRUE)

## Check future projections made by Models:

load("proj.Future1/Proj_Future1_GWWA")
Proj_Future1_GWWA[740:760,,1,1]

#################################################################################

PredictionBestModel(GLM=FALSE, CTA=TRUE,GAM=TRUE,GBM=TRUE,RF=TRUE,method='all',
Bin.trans=TRUE,Filt.trans=TRUE)

PredictionBestModel(GLM=FALSE, CTA=TRUE,GAM=TRUE,GBM=TRUE,RF=TRUE,method='Roc',
Bin.trans=TRUE,Filt.trans=TRUE)

load("pred/BestModelByRoc")
BestModelByRoc

load("pred/BestModelByKappa")
BestModelByKappa

load("pred/BestModelByTSS")
BestModelByTSS

#################################################################################

ProjectionBestModel(Proj.name='Future1',Bin.trans=TRUE, Filt.trans=TRUE, method='all')
load(("proj.Future1/Proj_Future1_BestModelByRoc")
dim(Proj_Future1_BestModelByRoc)
dimnames(Proj_Future1_BestModelByRoc)[-1]
Proj_Future1_BestModelByRoc[100:500,,"GWWA"]


ProjectionBestModel(Proj.name='bird',Bin.trans=TRUE, Filt.trans=TRUE, method='all')
load(("proj.Future1/Proj_Future1_BestModelByRoc")
dim(Proj_Future1_BestModelByRoc)
dimnames(Proj_Future1_BestModelByRoc)[-1]
Proj_Future1_BestModelByRoc[100:500,,"GWWA"]

################################################################################

Ensemble.Forecasting(Proj.name= "Future1", weight.method='Roc', PCA.median=T,
binary=T, bin.method='Roc', Test=F, decay=1.6, repetition.models=T)

Ensemble.Forecasting(Proj.name= "bird", weight.method='Roc', PCA.median=T,
binary=T, bin.method='Roc', Test=F, decay=1.6, repetition.models=T)








##OUTPUTS
##This function will be run for all the species at once. It will produce an object per species. These
##objects are arrays of three dimensions :

##The second dimension is the repetition runs and the third dimension is the consensus methods.
##There is also an object called "Total consensus Future1" that makes a single output out of all the
##repetitions.

load("proj.Future1/consensus_GWWA_Future1")
dim(consensus_GWWA_Future1)

dimnames(consensus_GWWA_Future1)[-1]

data <- consensus_GWWA_Future1
par(mfrow=c(2,10))
par(mar=c(0.6,0.6,2,0.6))
level.plot(DataBIOMOD[,8], CoorXY, show.scale=F, title='GWWA', cex=0.5)
level.plot(data[,1,1], CoorXY, show.scale=F, title='GWWA_mean', cex=0.5)

data <- consensus_GWWA_bird
par(mfrow=c(2,10))
par(mar=c(0.6,0.6,2,0.6))
level.plot(DataBIOMOD[,8], CoorXY, show.scale=F, title='GWWA', cex=0.5)
level.plot(data[,1,1], CoorXY, show.scale=F, title='GWWA_mean', cex=0.5)


load("proj.Future1/Total_consensus_Future1")
dim(Total_consensus_Future1)

dimnames(Total_consensus_Future1)[-1]

##Now the second dimension is the species. Let's see and plot some of these :

Total_consensus_Future1[1:20,,1]



data <- Total_consensus_Future1
par(mfrow=c(2,10))
par(mar=c(0.6,0.6,2,0.6))
level.plot(DataBIOMOD[,8], CoorXY, show.scale=F, title='GWWA', cex=0.5)
level.plot(data[,1,1], CoorXY, show.scale=F, title='GWWA_mean', cex=0.5)
level.plot(data[,1,2], CoorXY, show.scale=F, title='GWWA_weighted_mean', cex=0.5)
level.plot(data[,1,3], CoorXY, show.scale=F, title='GWWA_median', cex=0.5)
level.plot(data[,1,6], CoorXY, show.scale=F, title='GWWA_TSS_mean', cex=0.5)
level.plot(DataBIOMOD[,9], CoorXY, show.scale=F, title='BWWA', cex=0.5)
level.plot(data[,2,1], CoorXY, show.scale=F, title='BWWA_mean', cex=0.5)
level.plot(data[,2,2], CoorXY, show.scale=F, title='BWWA_weighted_mean', cex=0.5)
level.plot(data[,2,3], CoorXY, show.scale=F, title='BWWA_median', cex=0.5)
level.plot(data[,2,4], CoorXY, show.scale=F, title='BWWA_TSS_mean', cex=0.5)

#################################################################################
#################################################################################


##BIOMOD Using rasters instead of data






bio1= raster ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/bio1.asc')

bio5= raster ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/bio5.asc')


bio10= raster ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/bio10.asc')

bio15= raster ('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/bio15.asc')


GWWA_stack=stack(conifer, elev, stand_age, tmax5)

###  Another approach ###############################################################################
## Project BIOMOD results into a raster file in order to export to ArcGIS
library (raster)
library (BIOMOD)
library (rgdal)
library(maptools)


setwd('C:/Program Files/R/R2120/library/BIOMOD/R/Model_Runs/July14/Sept26/Modern/')

bird1<-read.csv("CoorXY_GW_BW_Current_1998_2009.csv", header=TRUE)

Sp.Env<-bird[c(1:2, 6:11)]
CoorXY<-a[1:2]
attach (CoorXY)


Models(GLM = F, TypeGLM = "poly", Test = "AIC",GBM = T, No.trees = 1000, GAM = T,
Spline = 3, CTA = T, CV.tree = 50, ANN = F, CV.ann = 2, SRE = F, quant=0.025, FDA = F,
MARS = F, RF = T, NbRunEval = 3, DataSplit = 80, Yweights=NULL, Roc = T, Optimized.Threshold.Roc = T,
Kappa = T, TSS=T, KeepPredIndependent = T, VarImport=5, 
NbRepPA=2, strategy="circles", coor=CoorXY, distance=2, nb.absences=200)

Projection.raster(RasterProj= GWWA_stack, Proj.name="GWWA Model", GLM=TRUE, BinKappa=TRUE, FiltKappa=TRUE)



