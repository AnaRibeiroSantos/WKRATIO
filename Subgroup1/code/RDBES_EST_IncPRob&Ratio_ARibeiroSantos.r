################################################################################################
# This code uses myExchangeFileH1 from ReadExchangeFileExample_WKRDBEST2 (Author:David Currie) ###
## This data is dummy data for H1 VS, FT, FO, SS, SA, FM and BV.
## This code does not present anything new. It was an opportunity for me to get to use RDBES format and play around with ratio and inclusion probabilities estimation.
## Ana Ribeiro Santos
########################################################################################

library(dplyr)

#load the data
source ("ReadExchangeFileExample_WKRDBEST2.r")
#lst(myExchangeFileH1)

# Because the data did not have sampled weight or total sample weight I had to estimated.
#add species name to FM table
myExchangeFileH1$FM$Species <- myExchangeFileH1$SA$SAspeciesCode[match(myExchangeFileH1$FM$SAid,myExchangeFileH1$SA$SAid)]
#head(myExchangeFileH1$FM)

#dim(myExchangeFileH1$FM)
#For this example we are goin to select just one species = 1015724 = Leuciscus aspius.
myExchangeFileH1$FM <- myExchangeFileH1$FM [myExchangeFileH1$FM$Species =="1015724",]

#calculate the weight per length. a and b paramters from Fisbase
a <- 0.01389
b <- 3

myExchangeFileH1$FM$Wt.at.len_g <- with (myExchangeFileH1$FM, ((as.numeric(FMclass)/10)^ b) * a) # this gives in grams, equation expected CM not MM hence divide by 10

#calculate the weight for each SAid (Sample ID)
sample_wght_SAid <- myExchangeFileH1$FM %>%
  group_by(SAid) %>%
  summarise (., SAsampleWeightLive_est = round (sum (Wt.at.len_g), 0)) %>% as.data.frame()

# merge it with SA table
myExchangeFileH1$SA <- merge (myExchangeFileH1$SA, sample_wght_SAid)
#head(myExchangeFileH1$SA)

#Calculate the totalweightlive based on the ratio_SA
myExchangeFileH1$SA$ratio_SA <- ifelse (is.na (myExchangeFileH1$SA$SAsampleWeightLive), with (myExchangeFileH1$SA, as.numeric (SAnumberTotal) / as.numeric (SAnumberSampled)),
                                        with (myExchangeFileH1$SA, as.numeric (SAtotalWeightLive) / as.numeric (SAsampleWeightLive)))

myExchangeFileH1$SA$SAtotalWeightLive_est <- with (myExchangeFileH1$SA, SAsampleWeightLive_est * ratio_SA)

##########################################################################################################################################################
# uses myExchnageFileH1 data and assumes simple random sampling within stages
# trips within domains (quarter, area and gear) and ages within lengths.
# In the code below, the hierarchy is as follows in order: VS, FT, "FO","SS","SA","FM","BV"


# set up the population totals
popDat <- as.data.frame(matrix(c(16,39,86,54, 32, 29, 1464760,3937180,6831707, 2876490, 3487500, 2345431),byrow=F,nrow=6,ncol=2))
names(popDat) <- c("tripNum","liveWt")
rownames(popDat) <-  c("Q2 27 GND", "Q3 27 DRB", "Q3 27 FPO", "Q4 27 FPO", "Q4 27 FYK", "Q4 27 GNS")

# calculate inclusion probabilities assuming SRS within strata
myExchangeFileH1$VS$VSincProb <- as.numeric (myExchangeFileH1$VS$VSnumberSampled)/ as.numeric (myExchangeFileH1$VS$VSnumberTotal)
myExchangeFileH1$BV$BVincProb <- as.numeric (myExchangeFileH1$BV$BVnumberSampled)/ as.numeric (myExchangeFileH1$BV$BVnumberTotal)
myExchangeFileH1$FM$FMincProb <- 1
myExchangeFileH1$SA$SAincProb <- as.numeric (myExchangeFileH1$SA$SAnumberSampled)/as.numeric (myExchangeFileH1$SA$SAnumberTotal) 
myExchangeFileH1$SS$SSincProb <- as.numeric (myExchangeFileH1$SS$SSnumberSampled)/as.numeric (myExchangeFileH1$SS$SSnumberTotal) 
myExchangeFileH1$FO$FOincProb <- as.numeric (myExchangeFileH1$FO$FOnumberSampled)/as.numeric (myExchangeFileH1$FO$FOnumberTotal) 
myExchangeFileH1$FT$FTincProb <- as.numeric (myExchangeFileH1$FT$FTnumberSampled)/as.numeric (myExchangeFileH1$FT$FTnumberTotal) 

# caclulate quarter from date, and then make a domain of quarter combined with area
# allocate domains to FM & BV as well
myExchangeFileH1$SA$FTid <- myExchangeFileH1$SS$FTid[match(myExchangeFileH1$SA$SSid,myExchangeFileH1$SS$SSid)]
myExchangeFileH1$SS$date <- myExchangeFileH1$FO$FOendDate[match(myExchangeFileH1$SS$FOid,myExchangeFileH1$FO$FOid)]
myExchangeFileH1$SA$date <- myExchangeFileH1$SS$date[match(myExchangeFileH1$SA$SSid,myExchangeFileH1$SS$SSid)]
myExchangeFileH1$SA$quarter <- paste("Q",(as.numeric(substr(myExchangeFileH1$SA$date,6,7))-1) %/% 3 + 1,sep="")
myExchangeFileH1$SA$domain <- paste(myExchangeFileH1$SA$quarter,myExchangeFileH1$SA$SAarea, myExchangeFileH1$SA$SAgear)
myExchangeFileH1$FM$domain <- myExchangeFileH1$SA$domain[match(myExchangeFileH1$FM$SAid,myExchangeFileH1$SA$SAid)]
myExchangeFileH1$BV$domain <- myExchangeFileH1$FM$domain[match(myExchangeFileH1$BV$FMid,myExchangeFileH1$FM$FMid)]
myExchangeFileH1$BV$numAtUnit <- 1

# add spp to FM
myExchangeFileH1$FM$spp <- myExchangeFileH1$SA$SAspeciesCode[match(myExchangeFileH1$FM$SAid,myExchangeFileH1$SA$SAid)]

# -----------------------------------------------------------------------
# a function to calculate inclusion probabilities for the units at the final stage
# of sampling given all the inclusion probabilities for the other stages
# -----------------------------------------------------------------------
getIncProb <- function(RDB,stages){
  nStages <- length(stages)
  if (any(stages %in% c("FM"))) {
    RDB[["FM"]][["FMincProb"]] <- 1
  }
  RDB[[stages[[1]]]][["incProb"]] <- RDB[[stages[[1]]]][[paste(stages[[1]],"incProb",sep="")]]
  for (i in 2:(nStages)) {
    indx <- RDB[[stages[[i]]]][[paste(stages[[i-1]],"id",sep="")]]
    indxPrev <- RDB[[stages[[i-1]]]][[paste(stages[[i-1]],"id",sep="")]]
    RDB[[stages[[i]]]][["incProbPrev"]] <- RDB[[stages[[i-1]]]][[paste("incProb",sep="")]][match(indx,indxPrev)]
    RDB[[stages[[i]]]][["incProb"]] <- RDB[[stages[[i]]]][["incProbPrev"]]*RDB[[stages[[i]]]][[paste(stages[[i]],"incProb",sep="")]]
  }
  return(RDB)
}


# -----------------------------------------------------------------------
# Number-at-Length
# -----------------------------------------------------------------------
# set up the stages in the sampling design used in the estimation
#stages <- list("FT","FO", "SS","SA","FM")
stages <- list("VS", "FT","FO", "SS","SA","FM")
#stages <- list("LE","SS","SA","FM")

# calculate the inclusion probabilities by length class given the other inclusion probabilities 
# in the higher stages
test <- getIncProb(myExchangeFileH1,stages)

# calculate a Horvitz Thompson estimate for total numbers at length by domain
# assuming srs within the domain - this is valid for the RDBshare data
estL <- tapply(as.numeric(test$FM$FMnumberAtUnit)/test$FM$incProb,list(test$FM$FMclass,test$FM$domain),sum)

#estL <- test$FM %>%
#  group_by(FMclass, domain, spp) %>%
#  summarise (., est.L = sum (as.numeric (FMnumberAtUnit)/incProb))

# calculate a HT estimate for total landed weight using the sampled landed weights
# assuming srs etc as before
estX <- tapply(test$SA$SAtotalWeightLive_est/test$SA$incProbPrev,list(test$SA$domain),sum)/1e3

# get the relevant population totals
popX <- popDat[match(names(estX),rownames(popDat)),"liveWt"]

# calculate the ratio estimates for numbers at length
estLR <- estL*matrix(rep(popX/estX,dim(estL)[1]),byrow=T,ncol=dim(estL)[2])


# compare the estLR with the ratio estimation using the ratio between the sum total sampled weight and the landings in the domain
test <- myExchangeFileH1
sumX <- tapply(myExchangeFileH1$SA$SAtotalWeightLive_est,list(myExchangeFileH1$SA$domain),sum) # sum te total sampled weight for each domain
testLR <- estL*matrix(rep(popX/sumX,dim(estL)[1]),byrow=T,ncol=dim(estL)[2]) # multiply the number at length by the ratio between landings and sampled weight

testLR
estLR
## They are super different!!! So im wondering if Im doing something wrong here!! The fact that Im working with made up data does not help at all!!

########################################################################################################################################
# JUST to test if using a different approach from Liz's code and function we would get the same number at length up to Vessel selection'
#############################################################################################################################################
VSinfo <- myExchangeFileH1$VS %>% transmute(VSid,SDid,
                                            VSincProb)%>%distinct()

FTinfo <-myExchangeFileH1$FT%>%transmute(FTid, VSid,  FTincProb )%>%distinct()


FOinfo<-myExchangeFileH1$FO%>%transmute(FOid,FTid, time=paste("Q",(as.numeric(substr(myExchangeFileH1$FO$FOendDate,6,7))-1) %/% 3 + 1,sep=""), 
                                        metier=FOgear, space=FOarea, FOincProb)%>%distinct()
FOinfo$domain <- paste (FOinfo$time, FOinfo$space, FOinfo$metier) 

SSinfo<-myExchangeFileH1$SS%>%transmute(SSid, FOid, SSincProb)%>%distinct()
#SAinfo<-myExchangeFileH1$SA%>%transmute(SAid,SSid,spp=SAspeciesCode,
#                              wsamp=SAsampleWeightMeasured,wsamptot=SAtotalWeightMeasured)%>%distinct()

SAinfo<- myExchangeFileH1$SA%>%transmute(SAid,SSid,spp=SAspeciesCode,SAincProb)%>%distinct()


FMinfo<-myExchangeFileH1$FM%>%transmute(FMid,SAid,len=FMclass,n=as.numeric(FMnumberAtUnit), FMincProb)

samp_est_stratum_incPro <-FMinfo%>% left_join (SAinfo)%>%left_join(SSinfo)%>%left_join(FOinfo)%>%left_join(FTinfo)%>% left_join(VSinfo) %>%
  group_by(domain, spp, len, n) %>%
  summarise (., final_prob = FMincProb * SAincProb * SSincProb * FOincProb * FTincProb * VSincProb) %>% 
  ungroup() %>%
  group_by(domain, spp, len) %>%
  summarise (., est.n = sum(n / final_prob))




