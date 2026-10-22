##################################################################################
###### R Script: Using network analysis to study behavioural phenotypes: ######### 
######           an example using domestic dogs (Canis lupus familiaris) #########
##################################################################################

# -----------------------------------------------------------------
##### Install relevant packages and load libraries
# -----------------------------------------------------------------

library(devtools)
#install.packages("Amelia")
library(Amelia)
# install_github("SachaEpskamp/qgraph")
library(qgraph)
# install.packages("ggplot2")
library(ggplot2)
# install_github("SachaEpskamp/bootnet")
library(bootnet)
# install.packages("effsize")
library(effsize)
# install.packages("ordinal")
library(ordinal)
# install.packages("lme4")
library(lme4)
# install.packages("plyr")
library(plyr)


#--------------------------------------------------------------------
## a function for creating multiple plots from ggplot figures
## found here: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# -----------------------------------------------------------------

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#--------------------------------------------------------------------
###### set working directory & load data
# -----------------------------------------------------------------

setwd(" ") # choose the place where the data files are saved

myData = read.csv2("rawData.csv")

# -----------------------------------------------------------------
## remove questions with more than 10% of zero responses
## -- missing values/zeros may arise because the question was not as relevant/clear
# -----------------------------------------------------------------

missmap(myData[ , c(4:46)])
checkZero = myData[,c(4:46)]
checkZero[ is.na(checkZero) ] <- 6
apply(checkZero, 2, function(x) sum(x==0)) / nrow(checkZero)
# COU, COH , AQS, BPH, DST
myData_rmZero = myData[ , -c(4,31,36, 41 , 44)] 

# check for individuals who had many zero values myData_short 
checkZero2 = myData_rmZero[,4:41]
checkZero2[ is.na(checkZero2) ] <- 6
apply(checkZero2, 1, function(x) sum(x==0)) / ncol(checkZero2)
# most individuals had a small number of missing answers (< 3 questions)
# --  remove #170 who had 18% of 0 values; their other responses may not be as reliable
myData_rmZeroP = myData_rmZero[-c(170), ]

# -----------------------------------------------------------------
## check if missingness is based on work type or handlers
# -----------------------------------------------------------------

# create a new DF containing details about missingness

missingDF = data.frame(fullMissing = apply(myData_rmZeroP[ , -c(1:3, 42)], 1, function(x) (sum(is.na(x))))
								+ apply(myData_rmZeroP[ , -c(1:3, 42)], 1, function(x) (sum(na.omit(x==0))))  ,
                              fullNotMissing = 38 - (apply(myData_rmZeroP[ , -c(1:3, 42)], 1, function(x) (sum(is.na(x))))
								+ apply(myData_rmZeroP[ , -c(1:3, 42)], 1, function(x) (sum(na.omit(x==0))))),
                              workType = myData_rmZeroP$workType ,
                              handlerType = myData_rmZeroP$handlerType ,
                              id = myData_rmZeroP$newID )

# is missingness dependent on dog type?
summary(fitFullMissP <- glm( cbind(fullMissing, fullNotMissing) ~ factor(workType) , 
                    data = missingDF  , family=binomial("logit")))

# mean and SD for percentage missingness per dog type
detectMiss = as.matrix( apply(myData_rmZeroP[myData_rmZeroP$workType==0 , -c(1:3, 42)], 1, function(x) (sum(is.na(x))))
                    +apply(myData_rmZeroP[ myData_rmZeroP$workType==0, -c(1:3, 42)], 1, function(x) (sum(na.omit(x==0)))))

mean(detectMiss/38*100) ; sd(detectMiss/38*100)

patrolMiss = as.matrix( apply(myData_rmZeroP[ myData_rmZeroP$workType==1 , -c(1:3, 42)], 1, function(x) (sum(is.na(x))))
                        +apply(myData_rmZeroP[ myData_rmZeroP$workType==1 , -c(1:3, 42)], 1, function(x) (sum(na.omit(x==0)))))

mean(patrolMiss/38*100) ; sd(patrolMiss/38*100)

# is missingness dependent on handlers completing repeated surveys?
missingDF_sub = subset(missingDF, missingDF$handlerType == "patrol" | missingDF$handlerType == "detection") 

summary(fitMissP_2 <- glmer( cbind(fullMissing, fullNotMissing) ~ factor(handlerType) 
                    + (1 | id) ,
                    data = missingDF_sub  , family=binomial("logit")))

ICC = 0.619 / (0.619 + pi^2/3 ) # Intraclass Correlation  

# -----------------------------------------------------------------
# Remaining 0s are converted to NAs to treat as missing
# -----------------------------------------------------------------

pScores = myData_rmZeroP[, c(4:41) ] 
pScores[ pScores == 0 ] <- NA

# -----------------------------------------------------------------
# Remove questions that have many values to impute (> 5%)
# -----------------------------------------------------------------

pScores_rmvNA = pScores
missmap(pScores_rmvNA)
apply(pScores_rmvNA, 2, function(x) sum(is.na(x))) / 170
# Greater than 5% : ACO, ODR, FOV, TT
pScores_rmv5pc = pScores_rmvNA[, -c(2,29:30,37)]

# add work & handlerID in again
pScores_rmv5pc2 = cbind( myData_rmZeroP$workType, myData_rmZeroP$newID , pScores_rmv5pc )


# -----------------------------------------------------------------
## check if any correlations are very high & theoretically similar indicating redundancy
# -----------------------------------------------------------------

pScores_cor = cor_auto(pScores_rmv5pc2[,-c(1:2)] , forcePD = TRUE ) # polychoric correlations
pScores_corAbs = abs(pScores_cor)
pScores_corAbs[ pScores_corAbs < 0.8 ] <- NA

## SLP-BAL, SA-GWL, QK-ACT, CON-ADP, OBD-REC have correlations above 0.8 
## and are conceptually similar. 
## KEEP THE MOST SPECIFIC QUESTIONS: SLP,GWL,ACT,ADP,REC
imputeReady = pScores_rmv5pc2[ , -c(6,15,21,23,17)]

# -----------------------------------------------------------------
## How much missingness is there to impute?
# -----------------------------------------------------------------

amountMiss = as.matrix( apply(imputeReady[,-c(1:2)], 2, function(x) (sum(is.na(x))))/nrow(imputeReady))

mean(amountMiss)*100 ; sd(amountMiss)*100

colnames(imputeReady) <- c("workType","newID",colnames(imputeReady[,-c(1,2)]))

# -----------------------------------------------------------------
## impute missing values
# -----------------------------------------------------------------

set.seed(123) # set seed to be reproducible

imputeData = amelia(imputeReady , m=15 , parallel = 'multicore', ords=colnames(imputeReady[,-c(1:2)]) , cs = "newID" )

imputeList = list(as.data.frame(imputeData$imputations[1]), as.data.frame(imputeData$imputations[2]) , as.data.frame(imputeData$imputations[3]) , 
                  as.data.frame(imputeData$imputations[4]) , as.data.frame(imputeData$imputations[5]) , as.data.frame(imputeData$imputations[6]) ,
                  as.data.frame(imputeData$imputations[7]) , as.data.frame(imputeData$imputations[8]) , 
                  as.data.frame(imputeData$imputations[9]), as.data.frame(imputeData$imputations[10]),as.data.frame(imputeData$imputations[11]),
                  as.data.frame(imputeData$imputations[12]), as.data.frame(imputeData$imputations[13]), as.data.frame(imputeData$imputations[14]),
                  as.data.frame(imputeData$imputations[15]))

# add colnames 
for( i in seq(imputeList) ) { 
  colnames(imputeList[[i]]) <-c("workType", "id",
                           "CB","FHO","SLP","PS","PSV","FL","CUR","PLA",
                           "FIT","SOC","DA", "GWL","REC","OCD","OCN",
                           "ACT","ADP","ANX","STR","QU","FSH","GUS",
                           "TC","CC","FDA","TOY","TUG","FoH","WIL") 
		}

# save imputed data
write.csv(imputeList[[1]] , "MultipleImputationdata_one")
write.csv(imputeList[[2]] , "MultipleImputationdata_two")
write.csv(imputeList[[3]] , "MultipleImputationdata_three")
write.csv(imputeList[[4]] , "MultipleImputationdata_four")
write.csv(imputeList[[5]] , "MultipleImputationdata_five")
write.csv(imputeList[[6]] , "MultipleImputationdata_six")
write.csv(imputeList[[7]] , "MultipleImputationdata_seven")
write.csv(imputeList[[8]] , "MultipleImputationdata_eight")
write.csv(imputeList[[9]] , "MultipleImputationdata_nine")
write.csv(imputeList[[10]] , "MultipleImputationdata_ten")
write.csv(imputeList[[11]] , "MultipleImputationdata2_11")
write.csv(imputeList[[12]] , "MultipleImputationdata3_12")
write.csv(imputeList[[13]] , "MultipleImputationdata2_13")
write.csv(imputeList[[14]] , "MultipleImputationdata2_14")
write.csv(imputeList[[15]] , "MultipleImputationdata2_15")

# compute average of all & then round and save
imputeList_average = Reduce("+" , imputeList) / length(imputeList)
write.csv(round(imputeList_average) , "finalMIdata_averaged_rounded.csv")

# -----------------------------------------------------------------

#######################################################################
###################### Network analyses ###############################
#######################################################################

# -----------------------------------------------------------------

# Load MI data set
netData = read.csv("finalMIdata_averaged_rounded.csv")
netData$handlerType = myData_rmZeroP$handlerType # myData_rmZeroP must still be loaded in the R workspace
netData = as.data.frame(netData)

# -----------------------------------------------------------------
## check independence of responses for handlers filling out more than one survey
# -----------------------------------------------------------------

netData_Rep = subset(netData, netData$handlerType=="patrol" | netData$handlerType=="detection") 

# compare random effect to a model without using LRT and delta AIC
# for each descriptor
  fit1 <- clmm(factor(x) ~ factor(workType) + ( 1 | id ) , data = netData_Rep , link = "logit" )
  fit2 <- clm(factor(x) ~ factor(workType)  , data = netData_Rep , link = "logit" )
  anova(fit1, fit2)

# remove: CB (small variance) & FHO, OCD, OCN, ANX, QU, TC , CC, TUG 
netData_ind = netData[ , -c(4, 16 , 17, 20, 22, 25 , 26, 29)]
# -----------------------------------------------------------------
# functions to create for mode and variation ratio
# -----------------------------------------------------------------

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

variationRatio = function( x ) { 
  mode = Mode(x)
  freq = count(x == mode)$freq[2]
  N = length(x)
  varRatio = 1 - (freq/N)
  varRatio
}

# -----------------------------------------------------------------
## create descriptive statistics
# -----------------------------------------------------------------

descStat = data.frame(descriptor = colnames(netData_ind[-c(1:2,24)]), 
                      mode = apply(netData_ind[, -c(1:2,24)], 2, FUN=Mode),
                      varRatio = apply(netData_ind[, -c(1:2,24)],2,FUN=variationRatio),
                      modePatrol = apply(netData_ind[netData_ind$workType==1, -c(1:2,24)], 2, FUN=Mode),
                      varRatioPatrol = apply(netData_ind[netData_ind$workType==1,-c(1:2,24)], 2 , FUN = variationRatio),
                      modeDetect = apply(netData_ind[netData_ind$workType==0, -c(1:2,24)], 2, FUN=Mode),
                      varRatioDetect = apply(netData_ind[netData_ind$workType==0,-c(1:2,24)], 2 , FUN = variationRatio)
                      )

descStat[,c(-1)] = round(descStat[,c(-1)], 3)

# remove CB for low variation relative to other descriptors
finalData = netData_ind[, -3] 

# -----------------------------------------------------------------
## create detection and patrol dog specific polychoric correlation matrices
# -----------------------------------------------------------------

detectDogs = finalData[finalData$workType==0, -c(1:2,23)]
detectData = cor_auto(detectDogs)
patrolDogs = finalData[finalData$workType==1, -c(1:2,23)]
patrolData = cor_auto(patrolDogs)

# -----------------------------------------------------------------
## create detection and patrol dog specific GGM networks
# -----------------------------------------------------------------

EBICdetect = EBICglasso(detectData , n = 54 , gamma = 0.65 )
# write.csv(EBICdetect , "DetectionDog_AssociationMatrix.csv")
EBICpatrol = EBICglasso(patrolData , n = 116 , gamma = 0.65  )
# write.csv(EBICpatrol , "PatrolDog_AssociationMatrix.csv")

# patrol dog network; Figure 1a
pNetwork = qgraph(EBICpatrol, layout="spring", 
                     labels=colnames(patrolData) , 
                     posCol="darkblue" , negCol = "darkgoldenrod",
                     font = 2, 
                     vsize= 6,
                    borders = TRUE
		)

# patrol dog centrality measures
ncFullp = centralityTable(pNetwork , standardized = FALSE)

# detection dog network; Figure 1b
dNetwork = qgraph(EBICdetect, layout="spring", 
                     labels=colnames(detectData) , 
                     posCol="darkblue" , negCol = "darkgoldenrod",
                     font = 2, 
                     vsize= 6 , borders = TRUE
		)

# detection dog centrality measures
ncFulld = centralityTable(dNetwork , standardized = FALSE)

# -----------------------------------------------------------------
## make and save Figure 1 
# -----------------------------------------------------------------

tiff(filename = "NetworkGraphs.tiff", width = 6000 , height = 4000 , res = 300 ,
     compression = "lzw")
Layout <- averageLayout(pNetwork, dNetwork)
layout(t(1:2))
qgraph(pNetwork, layout = Layout)
text(-1, 1.1, "A", cex = 4 , font = 2)
qgraph(dNetwork, layout = Layout)
text(-1, 1.1, "B", cex = 4 , font = 2)
dev.off()

# -----------------------------------------------------------------
## Patrol dog non-parametric, node-wise and subject-wise bootstrapping
# -----------------------------------------------------------------

# NB: these bootstraps can take a long time to compute. 
## Set nCores = number of your computer cores for faster performance

nCores = 2

patrolBoot_Nonparametric <- bootnet(patrolDogs, 
                               nBoot = 2000, 
                               nCores = nCores, 
                               estFun = EBICglasso  , 
                               estArgs = list(n = nrow(patrolDogs) , gamma = 0.65) , 
                               prepFun = cor_auto , 
                               prepArgs = list() , 
                               labels = colnames(patrolDogs) , 
                               type = "nonparametric" , 
                               computeCentrality = TRUE , 
                               propBoot = 1 ,  
                               replacement = TRUE 
                               )

#save results
write.csv(patrolBoot_Nonparametric$bootTable , "bootTablePatrolBoot_Nonparametric.csv")
write.csv(summary(patrolBoot_Nonparametric) , "centralitySummaryPatrolBoot_Nonparametric.csv")

patrolBoot_Nodewise <- bootnet(patrolDogs, 
                   nBoot = 2000, 
                   nCores = nCores , 
                   estFun = EBICglasso  , 
                   estArgs = list(n = nrow(patrolDogs) , gamma = 0.65) , 
                   prepFun = cor_auto , 
                   prepArgs = list() , 
                   labels = colnames(patrolDogs) , 
                   type="node",
                   subNodes = 3:(ncol(patrolDogs)-1) , 
                   computeCentrality = TRUE , 
                   propBoot = 1 ,  
                   replacement = TRUE
                   )

# save results
write.csv(patrolBoot_Nodewise$bootTable , "bootTablePatrolBoot_Nodewise.csv")
write.csv(summary(patrolBoot_Nodewise) , "centralitySummaryPatrolBoot_Nodewise.csv")

patrolBoot_Subjectwise <- bootnet(patrolDogs, 
                              nBoot = 2000, 
                              nCores = 2 , 
                              estFun = EBICglasso  , 
                              estArgs = list(n = nrow(patrolDogs) , gamma = 0.65) , 
                              prepFun = cor_auto , 
                              prepArgs = list() , 
                              labels = colnames(patrolDogs) , 
                              type="person",
                              subPersons = round(seq(0.25,0.95,length=10)*nrow(patrolDogs)) , 
                              computeCentrality = TRUE , 
                              propBoot = 1 ,  
                              replacement = TRUE
)
 
write.csv(patrolBoot_Dogwise$bootTable , "bootTablePatrolBoot_Subjectwise.csv")
write.csv(summary(patrolBoot_Dogwise) , "centralitySummaryPatrolBoot_Subjectwise.csv")
 
# -----------------------------------------------------------------
## Detection dog non-parametric, node-wise and subject-wise bootstrapping
# -----------------------------------------------------------------
# 

nCores = 2 

detectBoot_Nonparametric = bootnet( detectDogs,
                                    nBoot = 2000 ,
                                    nCores = 2 ,
                                    estFun = EBICglasso,
                                    estArgs = list(n=nrow(detectDogs) , gamma = 0.65) ,
                                    prepFun = cor_auto ,
                                    prepArgs = list() ,
                                    labels = colnames(detectDogs) ,
                                    type = "nonparametric" ,
                                    computeCentrality = TRUE ,
                                    propBoot = 1 ,
                                    replacement = TRUE
)


write.csv(detectBoot_Nonparametric$bootTable , "bootTableDetectBoot_Nonparametric.csv")
write.csv(summary(detectBoot_Nonparametric) , "centralitySummaryDetectBoot_Nonparametric.csv")


detectBoot_Nodewise = bootnet( detectDogs,  
                               nBoot = 2000 , 
                               nCores = nCores , 
                               estFun = EBICglasso, 
                               estArgs = list(n=nrow(detectDogs) , gamma = 0.65) ,
                               prepFun = cor_auto , 
                               prepArgs = list() ,
                               labels = colnames(detectDogs) , 
                               type = "node" , 
                               subNodes = 3:(ncol(detectDogs)-1) ,
                               computeCentrality = TRUE , 
                               propBoot = 1 ,  
                               replacement = TRUE
)  

write.csv(detectBoot_Nodewise$bootTable , "bootTableDetectBoot_Nodewise.csv")
write.csv(summary(detectBoot_Nodewise) , "centralitySummaryDetectBoot_Nodewise.csv")

detectBoot_Subjectwise = bootnet( detectDogs,  
                              nBoot = 2000 , nCores = 2 , 
                              estFun = EBICglasso, 
                              estArgs = list(n=nrow(detectDogs) , gamma = 0.65) ,
                              prepFun = cor_auto , 
                              prepArgs = list() ,
                              labels = colnames(detectDogs) , 
                              type = "person" , 
                              subPersons = round(seq(0.25,0.95,length=10)*nrow(detectDogs))  ,
                              computeCentrality = TRUE , 
                              propBoot = 1 ,  
                              replacement = TRUE
)  

write.csv(detectBoot_Subjectwise$bootTable , "bootTableDetectBoot_Subjectwise.csv")
write.csv(summary(detectBoot_Subjectwise) , "centralitySummaryDetectBoot_Subjectise.csv")

#-------------------------------------------------------------
## Create Figure 2

# "Confidence intervals" are also calculated here but not reported
# because they have been shown to be biased.
# For more details, see Epskamp et al. (2016): https://arxiv.org/pdf/1604.08462v1.pdf
# (Thanks to Sacha Epskamp for pointing this out in an initial review) 
#-------------------------------------------------------------

centNP = read.csv( "centralitySummaryPatrolBoot_Nonparametric.csv")
centND = read.csv( "centralitySummaryDetectBoot_Nonparametric.csv")

btwnValP = subset( centNP , centNP$type == "betweenness" )
btwnValP$nodeName =c("ACT","ADP","CUR","DA","FDA","FIT","FL","FoH","FSH","GUS","GWL",
                     "PLA","PS","PSV","REC","SLP","SOC","STR","TOY","WIL") 
btwnValP$CIlower2 = btwnValP$sample - (2 * btwnValP$sd)
btwnValP$CIupper2 = btwnValP$sample + (2 * btwnValP$sd)
strengthValP = subset( centNP , centNP$type == "strength" )
strengthValP$nodeName=c("ACT","ADP","CUR","DA","FDA","FIT","FL","FoH","FSH","GUS","GWL",
                        "PLA","PS","PSV","REC","SLP","SOC","STR","TOY","WIL") 
strengthValP$CIlower2 = strengthValP$sample - (2 * strengthValP$sd)
strengthValP$CIupper2 = strengthValP$sample + (2 * strengthValP$sd)

btwnValD = subset( centND , centND$type == "betweenness" )
btwnValD$nodeName =c("ACT","ADP","CUR","DA","FDA","FIT","FL","FoH","FSH","GUS","GWL",
                     "PLA","PS","PSV","REC","SLP","SOC","STR","TOY","WIL") 
btwnValD$CIlower2 = btwnValD$sample - (2 * btwnValD$sd)
btwnValD$CIupper2 = btwnValD$sample + (2 * btwnValD$sd)
strengthValD = subset( centND , centND$type == "strength" )
strengthValD$nodeName =c("ACT","ADP","CUR","DA","FDA","FIT","FL","FoH","FSH","GUS","GWL",
                         "PLA","PS","PSV","REC","SLP","SOC","STR","TOY","WIL") 
strengthValD$CIlower2 = strengthValD$sample - (2 * strengthValD$sd)
strengthValD$CIupper2 = strengthValD$sample + (2 * strengthValD$sd)

btwn = data.frame(nodeName = rep(btwnValD$nodeName,2),
                  sample = c(btwnValP$sample , btwnValD$sample),
                  CIlowers = c(btwnValP$CIlower2, btwnValD$CIlower2),
                  CIuppers = c(btwnValP$CIupper2 , btwnValD$CIupper2),
                  group = rep(c("Patrol","Detection"),each = 20) )
btwn$CIlowers2 = btwn$CIlowers
btwn$CIlowers2[btwn$CIlowers2 < 0 ] <- 0

write.csv(btwn , "betweennessCentrality_PatrolDetection.csv")

btwn$group = as.character(btwn$group)
btwn$group= factor( btwn$group, levels=unique(btwn$group))

btwnBarPlot = ggplot( btwn , aes(x = nodeName , y = sample ,fill=group)) +
  geom_bar(position=position_dodge(), stat="identity", 
           colour='black', width=0.7) +
  scale_fill_manual("" , values=c("#00BFC4","#F8766D") ) +
  scale_y_continuous(breaks=seq(0,160,20)) + 
  xlab("") +
  ylab("Betweenness\n") + 
  theme(axis.title.x=element_text(size=30),
        axis.title.y=element_text(size=26), 
        axis.text.x = element_text(size=18,face="bold"), 
        axis.text.y = element_text(size=20) , 
        legend.key.height=unit(2.5,"line"), 
        legend.text=element_text(size=24))

btwnBarPlot

strength = data.frame(nodeName = rep(strengthValD$nodeName,2),
                      sample = c(strengthValP$sample , strengthValD$sample),
                      CIlowers = c(strengthValP$CIlower2, strengthValD$CIlower2),
                      CIuppers = c(strengthValP$CIupper2 , strengthValD$CIupper2),
                      group = rep(c("Patrol","Detection"),each = 20) )
strength$CIlowers2 = strength$CIlowers
strength$CIlowers2[strength$CIlowers2 < 0 ] <- 0

write.csv(strength , "strengthCentrality_PatrolDetection.csv")

strength$group = as.character(strength$group)
strength$group= factor( strength$group, levels=unique(strength$group))

strengthBarPlot = ggplot( strength , aes(x = nodeName , y = sample ,fill=group)) +
  geom_bar(position=position_dodge(), stat="identity", 
           colour='black', width=0.7) +
  scale_fill_manual("" , values=c("#00BFC4","#F8766D") ) +
  scale_y_continuous(limits = c(0,2.5), breaks=seq(0,4,0.5)) + 
  xlab("\nDescriptor") +
  ylab("Strength\n") + 
  theme(axis.title.x=element_text(size=30),
        axis.title.y=element_text(size=26), 
        axis.text.x = element_text(size=18,face="bold"), 
        axis.text.y = element_text(size=20) , 
        legend.key.height=unit(2.5,"line"), 
        #legend.title=element_text(size=22) , 
        #panel.grid.major = element_line(colour="black"),
        legend.text=element_text(size=24))

strengthBarPlot

# save plot using the multiplot function
tiff("CentralityBarPlots.tiff" , width=5000 , height = 2500 , res=300 , 
     compression="lzw")
multiplot(btwnBarPlot, strengthBarPlot)
dev.off()

#-------------------------------------------------------------
##### Calculate Cliff's Delta effect sizes
#-------------------------------------------------------------

patrolBoots <- read.csv("bootTablePatrolBoot_Nonparametric.csv")
detectBoots <- read.csv("bootTableDetectBoot_Nonparametric.csv")

bootPatrolBtwn = patrolBoots[ patrolBoots$type == "betweenness" , c("value","name","id")]
bootDetectBtwn = detectBoots[ detectBoots$type == "betweenness" , c("value","name","id")]
bootPatrolStrength = patrolBoots[ patrolBoots$type == "strength" , c("value","name","id")]
bootDetectStrength = detectBoots[ detectBoots$type == "strength" , c("value","name","id")]

patrolBootCent = data.frame(node = bootPatrolBtwn$id , 
                            boot = bootPatrolBtwn$name , 
                            betweenness = bootPatrolBtwn$value , 
                            strength = bootPatrolStrength$value , 
                            meanCent = rowMeans(cbind(bootPatrolBtwn$value, bootPatrolStrength$value))
)

detectBootCent = data.frame(node = bootDetectBtwn$id , 
                            boot = bootDetectBtwn$name , 
                            betweenness = bootDetectBtwn$value , 
                            strength = bootDetectStrength$value , 
                            meanCent = rowMeans(cbind(bootDetectBtwn$value, bootDetectStrength$value))
)

nodes = c("ACT","ADP","CUR","DA","FDA","FIT","FL","FoH","FSH","GUS","GWL", 
          "PLA","PS","PSV","REC","SLP","SOC","STR","TOY","WIL")

# make data frames of effect sizes
CliffsDeltaBtwn <- data.frame( node = nodes , 
                               effBtwn = rbind(  
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: ACT", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: ACT", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: ADP", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: ADP", "betweenness" ]),
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: CUR", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: CUR", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: DA", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: DA", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: FDA", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: FDA", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: FIT", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: FIT", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: FL", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: FL", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: FoH", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: FoH", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: FSH", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: FSH", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: GUS", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: GUS", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: GWL", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: GWL", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: PLA", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: PLA", "betweenness" ]),
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: PS", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: PS", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: PSV", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: PSV", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: REC", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: REC", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: SLP", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: SLP", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: SOC", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: SOC", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: STR", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: STR", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: TOY", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: TOY", "betweenness" ]), 
                                             conf.level = 0.95, 
                                             use.normal = FALSE),
                                 cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: WIL", "betweenness"]), 
                                             as.numeric(detectBootCent[detectBootCent$node == "N: WIL", "betweenness" ]),
                                             conf.level = 0.95, 
                                             use.normal = FALSE)
                               ))

CliffsDeltaStrength <- data.frame( node = nodes , 
                                   effStrength = rbind(  
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: ACT", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: ACT", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: ADP", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: ADP", "strength" ]),
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: CUR", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: CUR", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: DA", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: DA", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: FDA", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: FDA", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: FIT", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: FIT", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: FL", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: FL", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: FoH", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: FoH", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: FSH", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: FSH", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: GUS", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: GUS", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: GWL", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: GWL", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: PLA", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: PLA", "strength" ]),
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: PS", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: PS", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: PSV", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: PSV", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: REC", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: REC", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: SLP", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: SLP", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: SOC", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: SOC", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: STR", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: STR", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: TOY", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: TOY", "strength" ]), 
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE),
                                     cliff.delta(as.numeric(patrolBootCent[patrolBootCent$node == "N: WIL", "strength"]), 
                                                 as.numeric(detectBootCent[detectBootCent$node == "N: WIL", "strength" ]),
                                                 conf.level = 0.95, 
                                                 use.normal = FALSE)
                                   ))


upperCI_indexes = seq(2,40,2)
lowerCI_indexes = seq(1,39,2)
CIbtwn = matrix(unlist(CliffsDeltaBtwn$effBtwn.conf.int))
CIstrength = matrix(unlist(CliffsDeltaStrength$effStrength.conf.int))

# final matrix of individual and mean eff sizes & CIs 
CliffsDelta_matrix = data.frame( node = nodes , 
                                 effSizeBtwn = unlist(CliffsDeltaBtwn$effBtwn.estimate),
                                 effSizeStrength = unlist(CliffsDeltaStrength$effStrength.estimate) , 
                                 meanEff = rowMeans(cbind(unlist(CliffsDeltaBtwn$effBtwn.estimate),
                                                          unlist(CliffsDeltaStrength$effStrength.estimate))) , 
                                 lowerCIBtwn = CIbtwn[lowerCI_indexes,] , 
                                 lowerCIStrength = CIstrength[lowerCI_indexes,] , 
                                 upperCIBtwn = CIbtwn[upperCI_indexes,] , 
                                 upperCIStrength = CIstrength[upperCI_indexes,] ,
                                 meanLowerCI = rowMeans(cbind(CIbtwn[lowerCI_indexes,],CIstrength[lowerCI_indexes,])),
                                 meanUpperCI = rowMeans(cbind(CIbtwn[upperCI_indexes,],CIstrength[upperCI_indexes,]))
) 

write.csv(CliffsDelta_matrix , "CliffsDeltaMatrix.csv")

#-------------------------------------------------------------
## create Figure 3
#-------------------------------------------------------------

forestPlot <- ggplot(CliffsDelta_matrix , aes(x=meanEff , y = node )) + 
  geom_point(size=4) + 
  geom_errorbarh(aes(xmin=meanLowerCI,xmax=meanUpperCI), height=0.35) + 
  geom_vline(xintercept = c(-0.25,0.25) , linetype="dashed" ) +
  scale_x_continuous(limits = c(-1 , 1), breaks=seq(-1,1,.2)) + 
  xlab("\nCliff's Delta") +
  ylab("Descriptor\n") + 
  theme(axis.title.x=element_text(size=34,vjust=1),
        axis.title.y=element_text(size=34), 
        axis.text.y = element_text(size=28,face="bold"), 
        axis.text.x = element_text(size=29) )

forestPlot

ggsave("MeanEffectSizeDifferences.tiff" , last_plot() ,  scale = 1.5 , width=13 , height = 7)

# -----------------------------------------------------------------
## Create Figure 4
# -----------------------------------------------------------------

patrolNWcor <- plot(patrolBoot_Nodewise, statistics=c("betweenness","strength"), perNode=F, area = FALSE) +
                    theme_gray(base_size = 30) +
                    scale_color_brewer(guide=FALSE, palette = "Set1", direction = -1) +
                    scale_x_continuous(limits=c(3,19), breaks=seq(3,19,by=2)) +
                    ylab("Average correlation") +
                    ggtitle("Patrol dogs: node-wise") +
                    guides(fill=FALSE) +
                    annotate("text", x = 10 , y = 0.27 , label="Betweenness", size=10 , color="#377eb8") +
                    annotate("text", x = 8 , y = 0.90 , label="Strength", size=10, col="#e41a1c" ) + 
                    annotate("text", x = 3 , y = 0.95 , label = "A" , size = 14 , col="black", fontface="bold")

detectNWcor <- plot(detectBoot_Nodewise, statistics=c("betweenness","strength"), perNode=F, area = FALSE) +
                    theme_gray(base_size = 30) +
                    scale_color_brewer(guide=FALSE, palette = "Set1", direction = -1) +
                    scale_x_continuous(limits=c(3,19), breaks=seq(3,19,by=2)) +
                    ylab("Average correlation") +
                    ggtitle(("Detection dogs: node-wise")) +
                    guides(fill=FALSE) +
                    annotate("text", x = 12 , y = 0.15 , label="Betweenness", size=10 , color="#377eb8") +
                    annotate("text", x = 12 , y = 0.6 , label="Strength", size=10, col="#e41a1c" ) + 
                    annotate("text", x = 3 , y = 0.95 , label = "B" , size = 14 , col="black", fontface="bold")

PatrolAxis <- round(seq(0.25,0.95,length=10)*nrow(patrolDogs)) 

patrolSWcor <- plot(patrolBoot_Dogwise, statistics=c("betweenness","strength"), perNode=F, area = FALSE) +
                    theme_gray(base_size = 30) +
                    scale_color_brewer(guide=FALSE, palette = "Set1", direction = -1) +
                    scale_x_continuous(breaks=PatrolAxis , limits=c(29,110)) +
                    guides(fill=FALSE) +
                    ggtitle("Patrol dogs: subject-wise") +
                    annotate("text", x = 83 , y = 0.25 , label="Betweenness", size=10 , color="#377eb8") +
                    annotate("text", x = 70 , y = 0.85 , label="Strength", size=10, col="#e41a1c" ) +
                    xlab("Sampled dogs") +
                    ylab("Average correlation") +
                    annotate("text", x = 29 , y = 0.95 , label = "C" , size = 14 , col="black", fontface="bold")

Detectaxis <- round(seq(0.25,0.95,length=10)*nrow(detectDogs))

detectSWcor <- plot(detectBoot_Dogwise, statistics=c("betweenness","strength"), perNode=F, area = FALSE) + 
                    theme_gray(base_size = 30) + 
                    scale_color_brewer(guide=FALSE, palette = "Set1", direction = -1) +
                    scale_x_continuous(breaks=c(Detectaxis), limits=c(14,51)) +
                    scale_y_continuous(limits=c(-1,1), breaks=seq(-1,1,by=0.5))+
                    ylab("Average correlation") +
                    xlab("Sampled dogs") +
                    ggtitle("Detection dogs: subject-wise") +
                    guides(fill=FALSE) + 
                    annotate("text", x = 32, y = 0.15 , label="Betweenness", size=10 , color="#377eb8") +
                    annotate("text", x = 42 , y = -0.32 , label="Strength", size=10, col="#e41a1c" ) + 
                    annotate("text", x = 14 , y = 0.95 , label = "D" , size = 14 , col="black", fontface="bold")

tiff("MultiplotCentrality.tiff", height=4500 , width = 5000, 
     compression = "lzw", res=300 )
multiplot( plotlist = list(patrolNWcor, patrolSWcor,
                           detectNWcor, detectSWcor) , cols = 2)
dev.off()

# -----------------------------------------------------------------
## Figure S1
# -----------------------------------------------------------------

patrolNWrank_btwn <- plot(patrolBoot_Nodewise, statistics=c("betweenness"), perNode=T, legendNcol = 1) + 
                      theme_gray(base_size = 28) + 
                      scale_x_continuous(limits=c(3,19), breaks=seq(3,19,by=2)) +
                      scale_y_continuous(limits=c(0,105), breaks = seq(0,105,by=15)) +
                      xlab("Sampled nodes") +
                      ylab("Betweenness") +
                      theme(legend.title=element_blank(), 
                            strip.text = element_blank(),
                            legend.key.size = unit(1.5, 'lines'))

ggsave("PatrolBootNodewise_betweenness.tiff", last_plot(), height= 8, width = 10)

patrolNWrank_str <- plot(patrolBoot_Nodewise, statistics=c("strength"), perNode=T, legendNcol = 1) + 
                          theme_gray(base_size = 28) + 
                          scale_x_continuous(limits=c(3,19), breaks=seq(3,19,by=2)) +
                          scale_y_continuous(limits=c(0,2), breaks = seq(0,2,by=0.5)) +
                          xlab("Sampled nodes") +
                          ylab("Strength") +
                          theme(legend.title=element_blank(), 
                                strip.text = element_blank(),
                                legend.key.size = unit(1.5, 'lines'))

ggsave("PatrolBootNodewise_strength.tiff", last_plot(), height= 8, width = 10)

detectNWrank_btwn <- plot(detectBoot_Nodewise, statistics=c("betweenness"), perNode=T, legendNcol = 1) + 
                          theme_gray(base_size = 28) + 
                          scale_x_continuous(limits=c(3,19), breaks=seq(3,19,by=2)) +
                          scale_y_continuous(limits=c(0,90), breaks = seq(0,90,by=15)) +
                          xlab("Sampled nodes") +
                          ylab("Betweenness") +
                          theme(legend.title=element_blank(), 
                                strip.text = element_blank(),
                                legend.key.size = unit(1.5, 'lines'))

ggsave("DetectBootNodewise_betweenness.tiff", last_plot(), height= 8, width = 10)



detectNWrank_str <- plot(detectBoot_Nodewise, statistics=c("strength"), perNode=T, legendNcol = 1) + 
                          theme_gray(base_size = 28) + 
                          scale_x_continuous(limits=c(3,19), breaks=seq(3,19,by=2)) +
                          scale_y_continuous(limits=c(0,3), breaks = seq(0,3,by=0.5)) +
                          xlab("Sampled nodes") +
                          ylab("Strength") +
                          theme(legend.title=element_blank(), 
                                strip.text = element_blank(),
                                legend.key.size = unit(1.5, 'lines'))

ggsave("DetectBootNodewise_strength.tiff", last_plot(), height= 8, width = 10)

# -----------------------------------------------------------------
## Figure S2
# -----------------------------------------------------------------

PatrolAxis <- round(seq(0.25,0.95,length=10)*nrow(patrolDogs)) 

patrolSWrank_btwn <- plot(patrolBoot_Dogwise, statistics=c("betweenness"), perNode=T, legendNcol = 1) + 
                          theme_gray(base_size = 28) + 
                          scale_x_continuous(breaks=c(PatrolAxis), limits=c(29,110)) +
                          scale_y_continuous(limits=c(0,105), breaks = seq(0,105,by=15)) +
                          xlab("Sampled dogs") +
                          ylab("Betweenness") +
                          theme(legend.title=element_blank(), 
                                strip.text = element_blank(),
                                legend.key.size = unit(1.5, 'lines'))

ggsave("PatrolBootDogwise_betweenness.tiff", last_plot(), height= 8, width = 10)

patrolSWrank_str <- plot(patrolBoot_Dogwise, statistics=c("strength"), perNode=T, legendNcol = 1) + 
                         theme_gray(base_size = 28) + 
                         scale_x_continuous(breaks=c(PatrolAxis), limits=c(29,116)) +
                         scale_y_continuous(limits=c(0,4), breaks = seq(0,4,by=0.5)) +
                         xlab("Sampled dogs") +
                         ylab("Strength") +
                         theme(legend.title=element_blank(), 
                               strip.text = element_blank(),
                               legend.key.size = unit(1.5, 'lines'))

ggsave("PatrolBootDogwise_strength.tiff", last_plot(), height= 8, width = 10)

# 

Detectaxis <- round(seq(0.25,0.95,length=10)*nrow(detectDogs))

detectSWrank_btwn <- plot(detectBoot_Dogwise, statistics=c("betweenness"), perNode=T, legendNcol = 1) + 
                          theme_gray(base_size = 28) + 
                          scale_x_continuous(breaks=c(Detectaxis), limits=c(14,51)) +
                          scale_y_continuous(limits=c(0,75), breaks = seq(0,75,by=15)) +
                          xlab("Sampled dogs") +
                          ylab("Betweenness") +
                          theme(legend.title=element_blank(), 
                                strip.text = element_blank(),
                                legend.key.size = unit(1.5, 'lines'))

ggsave("DetectBootDogwise_betweenness.tiff", last_plot(), height= 8, width = 10)

detectSWrank_str <- plot(detectBoot_Dogwise, statistics=c("strength"), perNode=T, legendNcol = 1) + 
                          theme_gray(base_size = 28) + 
                          scale_x_continuous(breaks=c(Detectaxis), limits=c(14,51)) +
                          scale_y_continuous(limits=c(0,4), breaks = seq(0,4,by=0.5)) +
                          xlab("Sampled dogs") +
                          ylab("Strength") +
                          theme(legend.title=element_blank(), 
                                strip.text = element_blank(),
                                legend.key.size = unit(1.5, 'lines'))


ggsave("DetectBootDogwise_strength.tiff", last_plot(), height= 8, width = 10)

#-------------------------------------------------------------
##### Number of edges across nonparamteric bootstrap 
#-------------------------------------------------------------

patrolBoots = read.csv("bootTablePatrolBoot_Nonparametric.csv")
detectBoots = read.csv("bootTableDetectBoot_Nonparametric.csv")

dEdges = detectBoots[detectBoots$type=="edge",] 
dEdges[dEdges == 0 ] <- NA
dEdges2 = na.omit(dEdges)
splitD = split(dEdges2 , f = dEdges2$name)
nEdgesD = lapply(splitD , nrow)

pEdges = patrolBoots[patrolBoots$type=="edge",] 
pEdges[pEdges == 0 ] <- NA
pEdges2 = na.omit(pEdges)
splitP = split(pEdges2 , f = pEdges2$name)
nEdgesP = lapply(splitP , nrow)

sd(unlist(nEdgesD)) ; sd(unlist(nEdgesP))

