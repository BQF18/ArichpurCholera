library(reshape)
##### load household level dataset
load("D:/cholerathesis/Data/MainData/Current/Qifang_manuscript_data/q1q4iccctrl.Rda")

##### list of variables
namelist_cncdc <- c("defsharing_only","defecate_pitl","dist_more10","q11_supply",
                    "q11_tube","q4water2","soap","boil","pplrm_di","cluster")

##### subset household level dataset
clustertrend_cncdc <- q1q4iccctrl[namelist_cncdc]  # subset variables

#### calculate spatial extent of each cluster  = "medianclustersize"
q1q4ctrlarealist <- data.frame() # empty dataframe that is going to store cluster + its cluster size
clusternames <- unique(q1q4iccctrl$cluster)  # all cluster IDs

for(i in 1:length(clusternames)) {
  #q1.onecluster <- q1q4iccctrl[q1q4iccctrl$cluster==clusternames[i],]
  #max.dist.onecluster <- max(dist(cbind(q1.onecluster$X,q1.onecluster$Y)))
  #q1q4ctrlarealist[i,1]<-max.dist.onecluster
  q1q4ctrlarealist[i,1]<-clusternames[i]   ###!
  median.dist.onecluster <- median(dist(cbind(q1.onecluster$X,q1.onecluster$Y))) 
  q1q4ctrlarealist[i,2]<-median.dist.onecluster     ###!
}

#q1q4ctrlarealist <- rename(q1q4ctrlarealist,c(V1="clustersize",V2="cluster",V3="medianclustersize"))
q1q4ctrlarealist <- rename(q1q4ctrlarealist,c(V1="cluster",V2="medianclustersize"))
summary(q1q4ctrlarealist$medianclustersize)


##### concordance within matched-sets

# all combinations of clusters 
clustertrend_cncdc <- q1q4iccctrl[c(namelist_cncdc,"dataid")] 
clustertrend_cncdc_m2 <- q1q4ctrlarealist
namelist_cncdc2 <- namelist_cncdc[!namelist_cncdc %in% "cluster"] # exclude cluster id from variable list

for(Var in  namelist_cncdc2) {
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
    clustertrend_cncdc_m2 [which(clustertrend_cncdc_m2$cluster == clus),Var] <- cncdc_count/ncol(junkcombn)
  }
}

###### 

#for(Var in  c(namelist_cncdc2)) {
  # prep for bootstrap
  set.seed(1234)
  boot.iter=1000
  bootlist.trend.within<- lapply(1:boot.iter, function(i)  sample(nrow(clustertrend_cncdc_m2),replace=TRUE)) # 500columns
  # bootstrap linear regression
  b1<- rep(NA,1000)
  b0<- rep(NA,1000)
  for (b in 1:boot.iter) {
    junk <- clustertrend_cncdc_m2[unlist(bootlist.trend.within[b]),]  
    b1[b] = lm(substitute(i ~ medianclustersize, list(i=as.name(Var))),data=junk)$coeff[2]
    b0[b] = lm(substitute(i ~ medianclustersize, list(i=as.name(Var))),data=junk)$coeff[1] 
  }
  quantile(b1,probs=c(0.025,0.975))*100
  #### only if we need to plot...comment out the rest
  junk <- b1 %o% clustertrend_cncdc_m2$medianclustersize  +  b0
  totalhat <- t(junk)   # imagine dist is cbind to the dataset
     
  # calculate CI at each dist. 95%CI for each row....for linear association 
  junk95ci <- apply(totalhat,1,quantile,probs=c(0.025,0.975),na.rm=TRUE)
  within.linear.boot<- cbind(clustertrend_cncdc_m2[c("cluster","clustersize","medianclustersize",Var)], t(junk95ci))
  within.linear.boot<-as.data.frame(within.linear.boot)
  names(within.linear.boot)
  names(within.linear.boot) <- c("cluster","clustersize","medianclustersize",
                                Var,"low","high")# contain dist. and CI
#}



##############################################################




#### trend: by variable vs. distance between centroid ####
#### stored in clustertrend_centroid_m 
#####
# create a centroid variable for each cluster with controls only
# variables wanted to include
centroidvarlist = namelist_cncdc
namelist_centroid = c("X","Y","longitude","latitude",centroidvarlist)
clustertrend_centroid <- q1q4iccctrl[namelist_centroid]  # subset variables
clustertrend_centroid$cluster <- as.numeric(clustertrend_centroid$cluster)
# includes the cluster centroid GPS info for each cluster 
# calculates the centroid for each cluster
clustertrend_centroid_m<-ddply(clustertrend_centroid,.(cluster),colwise(mean, na.rm=TRUE)  ) # nan b/c all answers are NA
# reformat "cluster"
clustertrend_centroid_m$cluster<-formatC(clustertrend_centroid_m$cluster, width = 4, format = "d", flag = "0") 


# method two: concordance pair #####

# for each combination of clusters, subset all households within the two clusters 
# check concordance value, and output that number 

# calculated these above in method 1 section: all combinations of clusters 
clustertrend_centroid$cluster<-formatC(clustertrend_centroid$cluster, width = 4, format = "d", flag = "0") 
# all combinations of clusters 
ccombn <- combn(clustertrend_centroid_m$cluster,2)
# transform to dataframe
ccombn.df <- as.data.frame(ccombn,stringsAsFactors=FALSE)   # 2rows
ccombnt.df <- as.data.frame(t(ccombn),stringsAsFactors=FALSE)   #2columns
ccombnt.df$dist <- 0

#for (Var in centroidvarlist[!centroidvarlist %in% "cluster"]) { 
  ccombnt.df$pdiff <- 0
  for (i in 1:ncol(ccombn.df)) {  # types of cluster combinations
    junka = clustertrend_centroid_m[which(clustertrend_centroid_m$cluster==ccombn.df[1,i]),2:3] # XY of one cluster
    junkb = clustertrend_centroid_m[which(clustertrend_centroid_m$cluster==ccombn.df[2,i]),2:3] # XY of another cluster
    ccombnt.df$dist[i] <- dist( rbind (junka, junkb) ) # calculate the distance btw centroid...and store in a column in ccombnt.df
    #subset cluster1
    #subset cluster2
    junkc = clustertrend_centroid[which(clustertrend_centroid$cluster==ccombn.df[1,i]),Var]
    junkd = clustertrend_centroid[which(clustertrend_centroid$cluster==ccombn.df[2,i]),Var]
    # for the non-NA values
    junkc<-junkc[!is.na(junkc)]
    junkd<-junkd[!is.na(junkd)]
    # loop through, calculate concordance
    countz=0
    for (j in 1:length(junkc)) {
      z<-junkd==junkc[j]
      countz<-length(z[z==TRUE])+countz
    }
    zz <- countz/ (length(junkc)*length(junkd))
    zz <- ifelse(zz==Inf, NA,zz)
    ccombnt.df$pdiff[i] <- zz
  }
  ccombnt.df <- rename(ccombnt.df, c(pdiff = Var))
#}
summary(ccombnt.df$dist)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#5.833  257.000  444.700  442.900  609.300 1040.000 
sd(ccombnt.df$dist)
# [1] 220.5144



# method for all variables ####
# for every variable...get a value
# average all the value/variable

junkm <- ccombnt.df[,4:ncol(ccombnt.df)]
ccombnt.df$total <- apply(junkm, 1, mean,na.rm=TRUE)


#### trend: by variable vs. distance between centroid..linear...bootstrap...3/4distance ####
# subset to variables needed for bootstrap
ccombnt.bs <- ccombnt.df[c("dist",Var)]
# subset to 3/4 of the total dist?
ccombnt.bs34<-ccombnt.bs[which(ccombnt.bs$dist<780),]
#ccombnt.df34<-ccombnt.df[which(ccombnt.df$dist<780),]

set.seed(1234)
bootlist.trend.btw <- lapply(1:boot.iter, function(i)  sample(nrow(ccombnt.bs34),replace=TRUE))

# bootstrap linear association
b1<- rep(NA,1000)
b0<- rep(NA,1000)
for (b in 1:boot.iter) {
  junk <- ccombnt.bs34[unlist(bootlist.trend.btw[b]),]  # for each iteration based on bs34
  b1[b] = lm(substitute(i ~ dist, list(i=as.name(Var))),data=junk)$coeff[2]
  b0[b] = lm(substitute(i ~ dist, list(i=as.name(Var))),data=junk)$coeff[1]
}
quantile(b1,probs=c(0.025,0.975))*100
junk <- b1 %o% ccombnt.bs34$dist  +  b0
totalhat <- t(junk)   # imagine dist is cbind to the datasest 1054*500

# calculate CI at each dist. 95%CI for each row....for linear association 
junk95ci <- apply(totalhat,1,quantile,probs=c(0.025,0.975),na.rm=TRUE)
between.linear.boot<- cbind(ccombnt.bs34$dist, t(junk95ci))
between.linear.boot<-as.data.frame(between.linear.boot)
names(between.linear.boot)
names(between.linear.boot) <- c("dist","low","high")# contain dist. and CI
# combine estimates onto between.linear.boot...dist becomes ordered with this line of code..unordered before
between.linear.boot <- merge( between.linear.boot , ccombnt.bs34, by="dist")



##### calculate optimal know
opfun <- function(k){
  #spline.lm <- lm(substitute(i ~ bs(dist,degree=1, knots=c(k)),list(i=as.name(Var))),
  #                data=between.linear.boot)
  spline.lm <- 
    lm(q4water2 ~ bs(dist,degree=1, knots=c(k)),
                  data=between.linear.boot)
  #plot(between.linear.boot$dist,predict(spline.lm),type="l")
  AIC(spline.lm)
}
optmin<-optimize(opfun, lower=min(between.linear.boot$dist),
         upper=max(between.linear.boot$dist),
         maximum=F)$minimum



ccombnt.bs <- ccombnt.df[c("dist",Var)]
# subset to 3/4 of the total dist?
ccombnt.bsopt<-ccombnt.bs[which(ccombnt.bs$dist<optmin),]

set.seed(1234)
bootlist.trend.btw <- lapply(1:boot.iter, function(i)  sample(nrow(ccombnt.bsopt),replace=TRUE))

# bootstrap linear association
b1<- rep(NA,1000)
b0<- rep(NA,1000)
for (b in 1:boot.iter) {
  junk <- ccombnt.bsopt[unlist(bootlist.trend.btw[b]),]  # for each iteration based on bs34
  b1[b] = lm(substitute(i ~ dist, list(i=as.name(Var))),data=junk)$coeff[2]
  b0[b] = lm(substitute(i ~ dist, list(i=as.name(Var))),data=junk)$coeff[1]
}
quantile(b1,probs=c(0.025,0.975))*100
splm<- lm(substitute(i ~ dist, list(i=as.name(Var))),data=ccombnt.bsopt)#$coeff[2]
lines(ccombnt.bsopt$dist,predict(splm),type="l")

# junk <- b1 %o% ccombnt.bsopt$dist  +  b0
# totalhat <- t(junk)   # imagine dist is cbind to the datasest 1054*500
# 
# # calculate CI at each dist. 95%CI for each row....for linear association 
# junk95ci <- apply(totalhat,1,quantile,probs=c(0.025,0.975),na.rm=TRUE)


######################
###¡¡plot linear btw + within
#
between.linear.boot$species <- "between matched-set concordance"
within.linear.boot$species <- "within matched-set concordance"

twotrend <- ggplot() +
  geom_point(data=between.linear.boot, aes_string(x="dist", y=Var) ,shape=19,size=1) +   
  geom_smooth(data=between.linear.boot, aes_string(x="dist", y=Var, color="species") , method="loess", size=1) +
  ylab("Concordance of\n type of a latrine") +
  xlab("Distance between matched-sets (m)") +
  #geom_smooth(data=within.linear.boot, aes_string(x="medianclustersize", y=Var, color="species"),method=lm,size=1,fill=NA) +
  scale_color_manual(values=c("within matched-set concordance"="blue", 
                            "between matched-set concordance"="red")) +
  #ylim(c(0,1.2))+
  scale_y_continuous(labels=c(0,0.2,0.4,0.6,0.8,1.0,1.2),breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),
                     limits=c(0,1.1))+
  theme(legend.title=element_blank(),
        legend.position="bottom",#c(.65,-.92),
        legend.direction="vertical",
        legend.background = element_rect(fill=NA),
        legend.text = element_text(size = 20),
        axis.title.y = element_text(size = rel(1.8),vjust=0.8),
        axis.title.x = element_text(size = rel(1.8),vjust=-0.8),
        axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18))
print(twotrend)


between.linear.boot
lms<-lm(formula=q11_supply ~ bs(dist,degree=1,knot=c(400)),data=between.linear.boot)
summary(lms)
plot(between.linear.boot$dist,predict(lms),type="l")
