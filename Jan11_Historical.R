##Ensemble Modeling GWWA Project, 1/3/2011


##set working directory
setwd('C:/Program Files/R/R2121/R2121/library/')
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

##----------------------------------------------------------------------------------------------------
## link BIOMOD Functions
setwd('C:/Program Files/R/R2121/R2121/library/BIOMOD/BIOMOD/R/')
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
source("BIOMOD-internal.R")
source("response.plot.R")
source("LoadModels.R")
source("LoadProj.R")
source("KappaSRE.R")
source("level.plot.R")
source("BinaryTransformation.R")
source("FilteringTransformation.R")
##################################################################################

####################1km Data   ########################################
 
setwd('C:/Program Files/R/R2121/R2121/library/BIOMOD/BIOMOD/R/Jan11_Historical/')
dir()

hist<-read.csv("CoorXY_hist.csv", header=TRUE)
CoorXY<-hist[,2:3]
head(CoorXY)
nrow(CoorXY)
plot(CoorXY[,1],CoorXY[,2])


goldhist<-read.csv("gw_hist.csv", header=TRUE)
length(goldhist)
nrow(goldhist)
head(goldhist)
names(goldhist)
attach(goldhist)
plot(goldhist[,2],goldhist[,3])

## Blue-winged warbler:

hist<-read.csv("CoorXY_hist.csv", header=TRUE)
CoorXY<-hist[,2:3]
head(CoorXY)
nrow(CoorXY)
plot(CoorXY[,1],CoorXY[,2])


bluehist<-read.csv("bw_hist.csv", header=TRUE)
length(bluehist)
nrow(bluehist)
head(bluehist)
names(bluehist)
attach(bluehist)
plot(bluehist[,2],bluehist[,3])

######################################################################################

level.plot(goldhist[,5], CoorXY[,1:2], show.scale=T,title='GWWA All Data')
level.plot(goldhist[,9], CoorXY[,1:2], show.scale=T,title='Elevation (m)')
points(goldhist[,5],XY=goldhist[,2:3],pch=2,cex=0.3,color='grey')
level.plot(goldhist[,6], CoorXY[,1:2], show.scale=T,title='Maximum May temperature (C)')
level.plot(goldhist[,11], CoorXY[,1:2], show.scale=T,title='Mean annual temperature (1/10C)')
level.plot(goldhist[,12], CoorXY[,1:2], show.scale=T,title='Minimum May temperaure (C)')
level.plot(goldhist[,13], CoorXY[,1:2], show.scale=T,title='Minimum June temperaure (C)')
level.plot(goldhist[,14], CoorXY[,1:2], show.scale=T, title='Minimum July temperature (C)')
multiple.plot(goldhist[,15:24], CoorXY[,1:2], cex=0.5)

#########################################################################################
Initial.State(sp.name="GWWA",Response = goldhist[,5], Explanatory = goldhist[,6:24],
IndependentResponse = goldhist[,5], IndependentExplanatory = goldhist[,6:24])
##########################################################################################
ls()
head(DataBIOMOD)
Biomod.material
##########################################################################################
Models(GLM = T, TypeGLM = "poly", Test = "AIC",GBM = T, No.trees = 2000, GAM = T,
Spline = 3, CTA = T, CV.tree = 50, ANN = T, CV.ann = 2, SRE = F, quant=0.025, FDA = T,
MARS = T, RF = T, NbRunEval = 3, DataSplit = 70, Yweights=NULL, Roc = T, 
Optimized.Threshold.Roc = T,Kappa = T, TSS=T, KeepPredIndependent = T, 
VarImport=5, NbRepPA=2, strategy="circles",coor=CoorXY,distance=0.5,nb.absences=500)

##############################################################################################
##############################################################################################


## Another version of Historical Data:
setwd('C:/Program Files/R/R2121/R2121/library/BIOMOD/BIOMOD/R/Jan9/Historical/')
dir()

bird1<-read.csv("CoorXY_GW_Historical.csv", header=TRUE)
CoorXY<-bird1[,2:3]
head(CoorXY)
nrow(CoorXY)
plot(CoorXY[,1],CoorXY[,2])


bird<-read.csv("GW_Historical.csv", header=TRUE)
length(bird)
nrow(bird)
head(bird)
names(bird)
attach(bird)
plot(bird[,2],bird[,3])

###################################################################################################

level.plot(bird[,5], CoorXY[,1:2], show.scale=T,title='GWWA All Data')
level.plot(bird[,9], CoorXY[,1:2], show.scale=T,title='Elevation (m)')
level.plot(bird[,6], CoorXY[,1:2], show.scale=T,title='Maximum May temperature (C)')
level.plot(bird[,11], CoorXY[,1:2], show.scale=T,title='Mean annual temperature (1/10C)')
level.plot(bird[,12], CoorXY[,1:2], show.scale=T,title='Minimum May temperaure (C)')
level.plot(bird[,13], CoorXY[,1:2], show.scale=T,title='Minimum June temperaure (C)')
level.plot(bird[,14], CoorXY[,1:2], show.scale=T,title='Minimum July temperaure (C)')

multiple.plot(bird[,7:14], CoorXY[,1:2],cex=0.5)##VARS.

#############  Histograms and data 1km ###################################
Initial.State(sp.name="GW",Response = bird[,5], Explanatory = bird[,6:14],
IndependentResponse = bird[,5], IndependentExplanatory = bird[,6:14])


Initial.State(sp.name="GW",Response = bird[,5], Explanatory = bird[,6:14],
IndependentResponse = bird[,5], IndependentExplanatory = bird[,6:14])


pseudo.abs(coor=bird[,2:3], status=bird[,6], strategy='random', env=data[,9:16], nb.points=1000,distance=5000, plot=T,
species.name= 'GW', acol='grey80', pcol='red', add.pres=T)

pseudo.abs(coor=bird[,2:3], status=bird[,6], strategy='sre', env=bird[,9:16], distance=5000, plot=T,
species.name= 'GW', nb.points=1000,acol='grey80', pcol='red', add.pres=T)

ls()
head(DataBIOMOD)
Biomod.material
#####################################################################

Models(GLM = T, TypeGLM = "poly", Test = "AIC",GBM = T, No.trees = 2000, GAM = T,
Spline = 3, CTA = T, CV.tree = 50, ANN = T, CV.ann = 2, SRE = F, quant=0.025, FDA = T,
MARS = T, RF = T, NbRunEval = 3, DataSplit = 70, Yweights=NULL, Roc = T, 
Optimized.Threshold.Roc = T,Kappa = T, TSS=T, KeepPredIndependent = T, 
VarImport=5, NbRepPA=2, strategy="circles",coor=CoorXY,distance=0.5,nb.absences=50)


###############################################################################
Evaluation.results.Roc
Evaluation.results.TSS
Evaluation.results.Kappa

VarImportance
###############################################################################
## Calculating the number of pseudo-absences (nb.absences) to run in model:
length(Biomod.PA.data$GW)
length(Biomod.PA.data$BW)
##the number of presence points for each dataset:
sum(bird[,6])
sum(bird[,7])
##Hence, the number of absences available for calibration:
length(Biomod.PA.data$GW)-sum(bird[,6])
length(Biomod.PA.data$BW)-sum(bird[,7])
###############################################################################

PredictionBestModel(GLM=TRUE, CTA=TRUE,GAM=TRUE,GBM=TRUE,RF=TRUE,ANN=TRUE, MARS=TRUE, SRE=F, FDA=TRUE,
method='all',Bin.trans=FALSE,Filt.trans=FALSE)

load("pred/BestModelByRoc")
BestModelByRoc

load("pred/BestModelByKappa")
BestModelByKappa

load("pred/BestModelByTSS")
BestModelByTSS
##################################################################################
load("pred/Pred_GW_indpdt")
Pred_GW_indpdt[100:500,"GAM",,]
level.plot(Pred_GW_indpdt[,"GAM",1,1], CoorXY, title="GWWA independent GAM")

load("pred/Pred_GW_indpdt")
Pred_GW_indpdt[100:500,"RF",,]
level.plot(Pred_GW_indpdt[,"RF",1,1], CoorXY, title="GWWA indepedent RF")

load("pred/Pred_GW_indpdt")
Pred_GW_indpdt[100:500,"GBM",,]
level.plot(Pred_GW_indpdt[,"GBM",1,1], CoorXY, title="GWWA")

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
load("pred/Pred_GW")
GW_GAM=cbind(CoorXY[Biomod.PA.sample[[1]][[2]],],Pred_GW[,"GAM",1,1])
write.csv(GW_GAM,"GW GAM.csv")

load("pred/Pred_GW")
GW_GBM=cbind(CoorXY[Biomod.PA.sample[[1]][[2]],],Pred_GW[,"GBM",1,1])
write.csv(GW_GBM,"GW GBM.csv")

load("pred/Pred_GW_indpdt")
GW_RF=cbind(Pred_GW_indpdt,CoorXY)
write.csv(GW_RF,"GW indepedent RF.csv")

load("pred/Pred_GW")
GW_RF=cbind(CoorXY[Biomod.PA.sample[[1]][[2]],],Pred_GW[,"RF",1,1])
write.csv(GW_RF,"GW RF.csv")

load("pred/Pred_GW")
GW_RF=cbind(CoorXY[Biomod.PA.sample[[1]][[2]],],Pred_GW[,"RF",1,1])
write.csv(GW_RF,"GW RF.csv")

load("pred/Pred_GW")
GW_RF=cbind(CoorXY[Biomod.PA.sample[[1]][[2]],],Pred_GW[,"RF",1,1])
write.csv(GW_RF,"Predicted Historical GW RF.csv")

load("pred/Pred_GW")
GW_GBM=cbind(CoorXY[Biomod.PA.sample[[1]][[1]],],Pred_GW[,"GBM",1,1])
write.csv(GW_GBM,"Predicted Historical GW GBM 2.csv")

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##GWWA Historial, Jan 12, 2012:

load("pred/Pred_GWWA")
GW_RF_Hist=cbind(CoorXY[Biomod.PA.sample[[1]][[2]],],Pred_GWWA[,"RF",1,1])
write.csv(GW_RF_Hist,"GW_RF_Historical.csv")

load("pred/Pred_GWWA")
GWWA_Hist_GAM=cbind(CoorXY[Biomod.PA.sample[[1]][[2]],],Pred_GWWA[,"GAM",1,1])
write.csv(GWWA_Hist_GAM,"GW_HIst_GAM.csv")

load("pred/Pred_GWWA")
GWWA_Hist_GBM=cbind(CoorXY[Biomod.PA.sample[[1]][[2]],],Pred_GWWA[,"GBM",1,1])
write.csv(GWWA_Hist_GBM,"GW_Hist_GBM.csv")

load("models/GWWA_RF_PA1")
GWWA_RF_PA1$importance
varImpPlot(GWWA_RF_PA1)

load("pred/Pred_GWWA")
data.used<-DataBIOMOD[Biomod.PA.sample$GWWA$PA1,"GW"]
calibrate.plot(data.used,Pred_GWWA[,"GBM",1,1]/1000)

load("pred/Pred_GWWA_indpdt")
level.plot(Pred_GWWA_indpdt[,"GAM",1,1],CoorXY,title="GWWA Historical GAM independent")
Pred_GWWA_indpdt[1:10,,,]

load("pred/Pred_GWWA_indpdt")
level.plot(Pred_GWWA_indpdt[,"RF",1,1],CoorXY,title="GWWA Historical RF independent")
Pred_GWWA_indpdt[1:10,,,]

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##Let's see what data has been used for calibration of the run PA1:

our.lines<-Biomod.PA.sample$GW$PA1
par(mfrow=c(1,2))
level.plot(DataBIOMOD[,"GW"],CoorXY,title='original data',show.scale=F)
level.plot(DataBIOMOD[our.lines,"GW"],CoorXY[our.lines,],title='PA1',show.scale=F)

###################################################################################

load("models/GW_GLM_PA1")
GW_GLM_PA1
summary(GW_GLM_PA1)
GW_GLM_PA1$anova
par(mfrow=c(2,2))
plot(GW_GLM_PA1)

load("pred/Pred_GW")
data.used<-DataBIOMOD[Biomod.PA.sample$GW$PA1,"GW"]
calibrate.plot(data.used,Pred_GW[,"GLM",1,1]/1000)

load("pred/Pred_GW")
data.used<-DataBIOMOD[Biomod.PA.sample$GW$PA1,"GW"]
calibrate.plot(data.used,Pred_GW[,"GBM",1,1]/1000)

load("models/GW_GBM_PA1")
plot(GW_GBM_PA1,i.var=1)
plot(GW_GBM_PA1,i.var=2)

load("models/GW_CTA_PA1")
names(GW_CTA_PA1)
plot(GW_CTA_PA1,margin=0.05)
text(GW_CTA_PA1,use.n=T)

load("models/GW_RF_PA1")
GW_RF_PA1$importance
varImpPlot(GW_RF_PA1)


load("pred/Pred_GW")
dim(Pred_GW)
Pred_GW[1:20,,1,1]

load("pred/Pred_GW_indpdt")
level.plot(Pred_GW_indpdt[,"CTA",1,1],CoorXY,title="GWWA Modern GAM independent")
Pred_GW_indpdt[1:10,,,]
#######################################################################################

calibrate.plot(data.used,Pred_GWWA_indpdt[,"GAM",1,1]/1000)

response.plot(GWWA_GAM_PA1,bird[15:20])

load("models/GWWA_GAM_PA1")
response.plot(GWWA_GAM_PA1, bird[15:20])

load("models/GWWA_RF_PA1")
response.plot(GWWA_RF_PA1, bird[15:20])

load("models/BWWA_GAM_PA1")
response.plot(BWWA_GAM_PA1,bird[15:20])

load("models/BWWA_RF_PA1")
response.plot(BWWA_RF_PA1,bird[15:20])

load("pred/Pred_GW_indpdt")
Pred_GWp_indpdt[1000:5000,"GAM",,]
level.plot(Pred_GW_indpdt[,"GAM",1,1], CoorXY, title="GWWA GAM")


load("pred/Pred_GW_indpdt")
Pred_GWp_indpdt[1000:5000,"GAM",,]
level.plot(Pred_GW_indpdt[,"GAM",1,1], CoorXY, title="GWWA GAM")


load("pred/Pred_GWWA_indpdt")
Pred_GWWA_indpdt[1000:5000,"RF",,]
level.plot(Pred_GWWA_indpdt[,"RF",1,1], CoorXY, title="GWWA RF")

load("pred/Pred_GWWA")
Pred_GWWA[100:5000,"RF",,]
level.plot(Pred_GWWA[,"RF",1,1], CoorXY, title="GWWA RF")

################################################################################
## Extract presence and absence predictions:

CurrentPred(GLM=TRUE, CTA=TRUE,GAM=TRUE,GBM=TRUE,RF=TRUE,ANN=TRUE,FDA=TRUE,MARS=TRUE,BinRoc=TRUE,
BinKappa=TRUE, BinTSS=TRUE)


load("pred/Pred_GWWA_BinRoc")
Pred_GWWA_BinRoc
GWWA_Hist_binary=cbind(Pred_GWWA_BinRoc, CoorXY)
write.csv(GWWA_Hist_GBM,"GW_Hist_BinRoc.csv")
######################################################################################
PredictionBestModel(GLM=T, GAM=T, GBM=T, ANN=F, FDA=F, MARS=F, RF=T, SRE=F, method='all', 
Bin.trans=T, Filt.trans=T)

load("pred/BestModelByRoc")
BestModelByRoc

load("pred/BestModelByKappa")
BestModelByKappa

load("pred/BestModelByTSS")
BestModelByTSS

#################################################################################
## Load future climate scenrio (A2a, B2a--variable names have to be the same in all 
## input files##

setwd('C:/Program Files/R/R2121/R2121/library/BIOMOD/BIOMOD/R/Jan9/')
modern<-read.csv("gw_mod_Jan9b.csv", header=TRUE)
names(modern)
attach(modern)

bird1<-read.csv("CoorXY_gw_modernb.csv", header=TRUE)

Projection(Proj = modern[,6:14], Proj.name='modern', GLM = T, GBM = T, GAM = T,
CTA = T, RF = T,ANN = T, MARS =T, FDA=T, BinRoc = T, BinKappa = T, BinTSS = T, FiltRoc = T, 
FiltKappa = T, FiltTSS = T, repetition.models=T)

## Check future projections made by Models:

load("proj.modern/Proj_modern_GW")
Proj_modern_GW[740:1000,,1,1]
Projected_modern_GW=cbind(Proj_modern_GW,CoorXY)

write.csv(Proj_modern_GW[,,3,1],CoorXY, file="GW projected modern.csv")
write.csv(Projected_modern_GW,file="GW projected modern v2.csv")
write.table(Proj_modern_GW[,,3,1],file="GW projected modern.txt")

interpolate(Projected_modern_GW)
#################################################################################

PredictionBestModel(GLM=TRUE, CTA=TRUE,GAM=TRUE,GBM=TRUE,RF=TRUE,method='all',
Bin.trans=FALSE,Filt.trans=FALSE)

load("pred/BestModelByRoc")
BestModelByRoc

load("pred/BestModelByKappa")
BestModelByKappa

load("pred/BestModelByTSS")
BestModelByTSS

#################################################################################

ProjectionBestModel(Proj.name='Modern',Bin.trans=TRUE, Filt.trans=TRUE, method='all')
load(("proj.Modern/Proj_Modern_BestModelByRoc")
dim(Proj_Modern_BestModelByRoc)
dimnames(Proj_Modern_BestModelByRoc)[-1]
Proj_Modern_BestModelByRoc[800:3100,,"GWWA"]

################################################################################

Ensemble.Forecasting(Proj.name= "modern", weight.method='Roc', PCA.median=T,
binary=T, bin.method='Roc', Test=F, decay=1.6, repetition.models=T)

##OUTPUTS
##This function will be run for all the species at once. It will produce an object per species. These
##objects are arrays of three dimensions :

##The second dimension is the repetition runs and the third dimension is the consensus methods.
##There is also an object called "Total consensus Future1" that makes a single output out of all the
##repetitions.

load("proj.modern/consensus_GW_modern")
dim(consensus_GW_modern)

dimnames(consensus_GWWA_p_Modern)[-1]

data <- consensus_GW_modern
par(mfrow=c(2,10))
par(mar=c(0.6,0.6,2,0.6))
level.plot(DataBIOMOD[,8], CoorXY, show.scale=T, title='GWWA', cex=0.5)
level.plot(Modern[,1,1], CoorXY, show.scale=F, title='GWWA_mean', cex=0.5)

dim(consensus_GW_modern)

write.table(consensus_GW_modern[,3,1],file="GW projected consensus.txt")

load("proj.Modern/Total_consensus_Modern")
dim(Total_consensus_Modern)

dimnames(Total_consensus_Modern)[-1]

##Now the second dimension is the species. Let's see and plot some of these :

Total_consensus_Modern[1:20,,1]



data <- Total_consensus_Modern
par(mfrow=c(2,10))
par(mar=c(0.6,0.6,2,0.6))
level.plot(DataBIOMOD[,12], CoorXY, show.scale=F, title='GWWA', cex=0.5)
level.plot(data[,1,1], CoorXY, show.scale=F, title='GWWA_mean', cex=0.5)
level.plot(data[,1,2], CoorXY, show.scale=F, title='GWWA_weighted_mean', cex=0.5)
level.plot(data[,1,3], CoorXY, show.scale=F, title='GWWA_median', cex=0.5)
level.plot(data[,1,6], CoorXY, show.scale=F, title='GWWA_TSS_mean', cex=0.5)
level.plot(DataBIOMOD[,13], CoorXY, show.scale=F, title='BWWA', cex=0.5)
level.plot(data[,2,1], CoorXY, show.scale=F, title='BWWA_mean', cex=0.5)
level.plot(data[,2,2], CoorXY, show.scale=F, title='BWWA_weighted_mean', cex=0.5)
level.plot(data[,2,3], CoorXY, show.scale=F, title='BWWA_median', cex=0.5)
level.plot(data[,2,4], CoorXY, show.scale=F, title='BWWA_TSS_mean', cex=0.5)

#################################################################################
#################################################################################


names(bird)

summary(bird$elev)
summary(bird$elev~bird[,6]=="2")
x<-bird[,21]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$elev[bird$Group=="1"],breaks= c(386:1132),main="GWWA on elev",xlab="elev",col="yellow")
hist(bird$elev[bird$Group=="2"],breaks= c(386:1132),main="GWWA A on elev",xlab="elev",col="blue")

summary(bird$ORNL)
summary(bird$ORNL~bird[,6]=="2")
x<-bird[,18]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$ORNL[bird$Group=="2"],breaks= c(0:802),main="GWWA absence on ornl",xlab="ornl",col="blue")
hist(bird$ORNL[bird$Group=="1"],breaks= c(0:802),main="GWWA on ornl",xlab="ornl",col="yellow")

summary(bird$USGS)
summary(bird$USGS~bird[,6]=="2")
x<-bird[,20]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$USGS[bird$Group=="1"],breaks= c(1:24),main="GWWA on usgs",xlab="usgs",col="yellow")
hist(bird$USGS[bird$Group=="2"],breaks= c(1:24),main="GWWA absence on usgs",xlab="usgs",col="blue")

summary(bird$Prcnt_decid)
summary(bird$Prcnt_decid~bird[,6]=="2")
x<-bird[,16]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$Prcnt_decid[bird$Group=="2"],breaks= c(0:80),main="GWWA absence on % dec",xlab="ornl",col="blue")
hist(bird$Prcnt_decid[bird$Group=="1"],breaks= c(0:80),main="GWWA on % dec",xlab="ornl",col="yellow")

summary(bird$slope)
summary(bird$slope~bird[,6]=="2")
x<-bird[,21]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$slope[bird$Group=="1"],breaks= c(0:17),main="GWWA on % dec",xlab="ornl",col="yellow")

summary(bird$Prcnt_evrgreen)
summary(bird$Prcnt_evrgreen~bird[,6]=="2")
x<-bird[,16]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$Prcnt_evrgreen[bird$Group=="2"],breaks= c(0:57),main="GWWA absence on % ev",xlab="ornl",col="blue")
hist(bird$Prcnt_evrgreen[bird$Group=="1"],breaks= c(0:57),main="GWWA on % ev",xlab="ornl",col="yellow")

summary(bird$Prcnt_dec_evergrn)
summary(bird$Prcnt_dec_evergrn~bird[,6]=="2")
x<-bird[,18]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$Prcnt_dec_evergrn[bird$Group=="2"],breaks= c(0:80),main="GWWA absence on % dec/ev",xlab="ornl",col="blue")
hist(bird$Prcnt_dec_evergrn[bird$Group=="1"],breaks= c(0:80),main="GWWA on % dec/ev",xlab="ornl",col="yellow")

summary(bird$Fraggy)
summary(bird$Fraggy~bird[,6]=="2")
x<-bird[,23]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$Fraggy[bird$Group=="2"],breaks= c(0:9),main="GWWA absence on Frag",xlab="Frag",col="blue")
hist(bird$Fraggy[bird$Group=="1"],breaks= c(0:9),main="GWWA on Frag ",xlab="Frag",col="yellow")

summary(bird$Distance_perf2)
summary(bird$Distance_perf2~bird[,6]=="2")
x<-bird[,23]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$Distance_perf2[bird$Group=="1"],breaks= c(0:5),main="GWWA on Frag ",xlab="Frag",col="yellow")

summary(bird$Prcnt_Int)
summary(bird$Prcnt_Int~bird[,6]=="2")
x<-bird[,27]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)

summary(bird$Prcnt_trans)
summary(bird$Prcnt_trans~bird[,6]=="2")
x<-bird[,28]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)


summary(bird$Dist_edge2)
summary(bird$Dist_edge2~bird[,6]=="2")
x<-bird[,37]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)

summary(bird$Prcnt_patch)
summary(bird$Prcnt_patch~bird[,6]=="2")
x<-bird[,38]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)


summary(bird$Dist_Int2)
summary(bird$Dist_Int2~bird[,6]=="2")
x<-bird[,30]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)

summary(bird$Prcnt_edge)
summary(bird$Prcnt_edge~bird[,6]=="2")
x<-bird[,35]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)

summary(bird$Dist_patch2)
summary(bird$Dist_patch2~bird[,6]=="2")
x<-bird[,40]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)

summary(bird$EVH)
summary(bird$EVH~bird[,6]=="2")
x<-bird[,45]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$EVH[bird$Group=="2"],breaks= c(11:111),main="GWWA absence on veg height",xlab="Veg height class",col="blue")
hist(bird$EVH[bird$Group=="1"],breaks= c(21:111),main="GWWA on veg height",xlab="Veg height class",col="yellow")


summary(bird$EVC)
summary(bird$EVC~bird[,6]=="2")
x<-bird[,46]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$EVC[bird$Group=="2"],breaks= c(11:127),main="BWWA on veg cover ",xlab="Veg cover class",col="blue")
hist(bird$EVC[bird$Group=="1"],breaks= c(21:128),main="GWWA on veg cover",xlab="Veg cover class",col="yellow")

summary(bird$Max_Temp)
summary(bird$Max_Temp~bird[,6]=="2")
x<-bird[,47]
y<-bird[,6]
myline.fit<-lm(y~x)
summary(myline.fit)
hist(bird$Max_Temp[bird$Group=="2"],breaks= c(23:30),main="BWWA on max temp ",xlab="Max temp",col="blue")
hist(bird$Max_Temp[bird$Group=="1"],breaks= c(23:30),main="GWWA on max temp",xlab="Max temp",col="yellow")

##############################################################################################

l<-vector(length=2)
l[1]="Present"
l[2]="Absent"
bird$Group<-factor(bird$Group,labels=l)

## Set defaults
## 2-to-1 cost matrix
clmat<-matrix(nrow=2,ncol=2)
clmat[1,1]<-0; clmat[2,1]<-1; clmat[1,2]<-2; clmat[2,2]<-0

## bucket and node size
##maxdepth controls how many layers are in final tree
my.control<-c(minsplit=10,minbucket=5,cp=0.001,maxcompete=2,maxsurrogate=1,
usesurrogate=2,xval=10,surrogatestyle=0,maxdepth=5)

#####################################################################################

## 1km Format ##############################

names(bird)

tree<-rpart(Group~land+Prcnt_tree+Lnd_diversity+Prcnt_decid+
Prcnt_evrgreen+Prcnt_dec_evergrn+ORNL+USGS+slope+elev+aspect+Fraggy+
Prcnt_canopy+Prcnt_Int+Prcnt_trans+Dist_Int2+Dist_perf2
+Dist_trans2+Dist_edge2+Prcnt_patch+Dist_patch2+
Dist_stream2+Dist_any_forest2+Dist_river2+EVH+EVC,
data=bird,method="class",control=my.control,
parms=list(loss=clmat))

tree<-rpart(Group~land+Prcnt_tree+Lnd_diversity2+Prcnt_decid+
Prcnt_evrgreen+Prcnt_dec_evergrn+ORNL+USGS+slope+elev+aspect+Fraggy+
Prcnt_canopy2+Prcnt_Int+Prcnt_trans+Dist_Int2+Dist_perf2
+Dist_trans2+Dist_edge2+Prcnt_patch+Dist_patch2+
Dist_stream2+Dist_any_forest2+Dist_river2+EVH+EVC,
data=bird,method="class",control=my.control,
parms=list(loss=clmat))


####################################################################################

#####################################################################################################


# Low-resolution plot
dev.set(which=1); plot(tree); text(tree)

# High-resolution plot
dev.set(which=2); post(tree,filename="")

# Cross-validation plots
dev.set(which=3); par(mfcol=c(2,1)); rsq.rpart(tree); par(mfcol=c(1,1))

# Tree pruning
pruned<-prune(tree,cp=0.0351)
dev.set(which=1); plot(pruned); text(pruned)



