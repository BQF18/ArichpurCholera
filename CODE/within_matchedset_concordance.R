library(reshape)
##### load household level dataset
load("D:/cholerathesis/Data/MainData/Current/Qifang_manuscript_data/q1q4iccctrl.Rda")

##### list of variables
varlist <- c("defsharing_only","defecate_pitl","dist_more10","q11_supply",
                    "q11_tube","q4water2","soap","boil","pplrm_di","cluster")

##### subset household level dataset
clustertrend_cncdc <- q1q4iccctrl[varlist]  # subset variables

#### calculate spatial extent of each cluster  = "medianclustersize"
q1q4ctrlarealist <- data.frame() # empty dataframe that is going to store cluster + its cluster size
clusternames <- unique(q1q4iccctrl$cluster)  # all cluster IDs

for(i in 1:length(clusternames)) {
  q1.onecluster <- q1q4iccctrl[q1q4iccctrl$cluster==clusternames[i],]
  q1q4ctrlarealist[i,1]<-clusternames[i]  
  median.dist.onecluster <- median(dist(cbind(q1.onecluster$X,q1.onecluster$Y))) 
  q1q4ctrlarealist[i,2]<-median.dist.onecluster   
}

q1q4ctrlarealist <- rename(q1q4ctrlarealist,c(V1="cluster",V2="medianclustersize"))
summary(q1q4ctrlarealist$medianclustersize)


##### concordance within matched-sets

# all combinations of clusters 
clustertrend_cncdc <- q1q4iccctrl[c(varlist,"dataid")] 
within_hh_dat <- q1q4ctrlarealist

within_hh_function <- function(Var){
  for (clus in unique(clustertrend_cncdc$cluster)) {
    # For each cluster, count household pairs having the same exposure in the same matched-set = "cncdc_count"
    junksubset <- clustertrend_cncdc[which(clustertrend_cncdc$cluster==clus),c("dataid",Var)]
    # scenario 1: not a sinle pair of hh within a matched-set that both are non-NA
    if (  length(junksubset[,2][is.na(junksubset[,2])==FALSE]) <2) {
      cncdc_count = NA
    }
    else {
      # discard the hh in the matchedset that answered NA for this variable
      junksubset2 <- junksubset[!is.na(junksubset[,2]),]
      junkcombn <- combn(junksubset2[,1],2)
      # and calculate the total pairs of households that have the same exposure
      cncdc_count=0
      for (j in 1:ncol(junkcombn)) {
        junk_count <- ifelse(clustertrend_cncdc[which(clustertrend_cncdc$dataid==junkcombn[1,j]),Var]  == 
                               clustertrend_cncdc[which(clustertrend_cncdc$dataid==junkcombn[2,j]),Var] , 1, 0)
        cncdc_count = junk_count + cncdc_count
      }
    }
    # in terms of proportion
    within_hh_dat [which(within_hh_dat$cluster == clus),Var] <- cncdc_count/ncol(junkcombn)
  }
   
  # bootstrap
  set.seed(1234)
  boot.iter=1000
  bootlist.trend.within<- lapply(1:boot.iter, function(i)  sample(nrow(within_hh_dat),replace=TRUE))
  # bootstrap linear regression
  b1<- rep(NA,1000)
  b0<- rep(NA,1000)
  for (b in 1:boot.iter) {
    junk <- within_hh_dat[unlist(bootlist.trend.within[b]),]  
    b1[b] = lm(substitute(i ~ medianclustersize, list(i=as.name(Var))),data=junk)$coeff[2]
  }
  quantile(b1,probs=c(0.025,0.975))*100
}

within_hh_function(Var=varlist[1])
