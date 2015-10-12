library(reshape)
##### load household level dataset, households in matched-sets only
load("Data/q1q4_ctrl.Rda")

#### de-identified the household data and deleted the GPS coordinates of each household.
#### so the following 10 lines of code would NOT run.
#### load the following data instead
load("Data/q1q4spatialextent.Rda")
# #### calculate spatial extent of each cluster  = "medianclustersize"
# q1q4spatialextent <- data.frame() # empty dataframe that is going to store cluster + its cluster size
# clusternames <- unique(q1q4_ctrl$cluster)  # all cluster IDs
# 
# for(i in 1:length(clusternames)) {
#   q1.onecluster <- q1q4_ctrl[q1q4_ctrl$cluster==clusternames[i],]
#   q1q4spatialextent[i,1]<-clusternames[i]  
#   median.dist.onecluster <- median(dist(cbind(q1.onecluster$X,q1.onecluster$Y))) 
#   q1q4spatialextent[i,2]<-median.dist.onecluster   
# }
# q1q4spatialextent <- rename(q1q4spatialextent,c(V1="cluster",V2="medianclustersize"))


##### FUNCTION: concordance within matched-sets

within_hh_dat <- q1q4spatialextent

within_hh_function <- function(Var,
                               boot.iter=1000){
  for (clus in unique(q1q4_ctrl$cluster)) {
    # In each matched-set, count household pairs having the same exposure in the same matched-set = "cncdc_count"
    junksubset <- q1q4_ctrl[which(q1q4_ctrl$cluster==clus),c("dataid",Var)]
    # If a matched-set contains zero pair of hhs that both answered the question 
    if (  length(junksubset[,2][is.na(junksubset[,2])==FALSE]) <2) {
      cncdc_count = NA
    }
    else {
      # if at least one pair that both answered the question
      # discard the hh in the matched-set that answered NA for this variable
      junksubset2 <- junksubset[!is.na(junksubset[,2]),]
      junkcombn <- combn(junksubset2[,1],2)
      # and calculate the total pairs of households that have the same exposure
      cncdc_count=0
      for (j in 1:ncol(junkcombn)) {
        junk_count <- ifelse(q1q4_ctrl[which(q1q4_ctrl$dataid==junkcombn[1,j]),Var]  == 
                               q1q4_ctrl[which(q1q4_ctrl$dataid==junkcombn[2,j]),Var] , 1, 0)
        cncdc_count = junk_count + cncdc_count
      }
    }
    # get proportion
    within_hh_dat [which(within_hh_dat$cluster == clus),Var] <- cncdc_count/ncol(junkcombn)
  }
   
  # bootstrap
  set.seed(1234)
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




##### CALCULATE concordance within matched-sets for each variable
# list of variables
varlist <- c("defsharing_only","defecate_pitl","dist_more10","q11_supply",
             "q11_tube","q4water2","soap","boil","pplrm_di")

within_hh_function(Var=varlist["defsharing_only"])

