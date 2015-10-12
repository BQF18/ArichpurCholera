library(plyr)
library(splines)
#### load either household level dataset or individual level data
load("Data/q1q4_ctrl.Rda")
# or ..
load("Data/q2_ctrl.Rda")

###### visualize between matched-set concordance vs. ..
###### distance between matched-set with LOESS.  (Figure 3 in the manuscript)
###### Var=exposure
###### Varlabel = label of exposure in the plot
###### Data = either the household or individual level dataset
between.linear.boot<-between_function(Varlabel="wulala",Var="defsharing_only",Data=q1q4iccctrl)

##### calculate the "knot" for exposures showing clustering beyond matched-set

names(between.linear.boot)[4] <- "Fakename"  # rename it
optim <- optimize(opfun,Data=between.linear.boot,lower=min(between.linear.boot$dist),
                 upper=max(between.linear.boot$dist),
                 maximum=F)$minimum
names(between.linear.boot)[4] <- "defsharing_only"  # name it back

##### find the linear association between concordance and distance 
##### before concordance levels off...(Results in Figure 2E)
lmknot(optim=optim,Var="defsharing_only",boot.iter=1000,Data=between.linear.boot)




########  ALL FUNCTIONS USED ABOVE
#### visualize between matched-set concordance vs. distance between matched-set with LOESS
between_function <- function(Var,Varlabel,Data){
  # format data
  Data$cluster <- as.numeric(Data$cluster)
  #### distance between centroid
  clustertrend_centroid_m <- ddply(Data,.(cluster),colwise(mean, na.rm=TRUE)  ) # nan b/c all answers are NA
  # format "cluster"
  clustertrend_centroid_m$cluster<-formatC(clustertrend_centroid_m$cluster, width = 4, format = "d", flag = "0")
  Data$cluster<-formatC(Data$cluster, width = 4, format = "d", flag = "0") 
  # all combinations of matched-sets 
  ccombn <- combn(clustertrend_centroid_m$cluster,2)
  # transform to dataframe
  ccombnt.df <- as.data.frame(t(ccombn),stringsAsFactors=FALSE)   #2columns
  ccombnt.df$dist <- 0
  ccombnt.df$pdiff <- 0
  for (i in 1:nrow(ccombnt.df)) { 
    junka = clustertrend_centroid_m[which(clustertrend_centroid_m$cluster==ccombn[1,i]),c("X","Y")] # XY of one cluster
    junkb = clustertrend_centroid_m[which(clustertrend_centroid_m$cluster==ccombn[2,i]),c("X","Y")] # XY of another cluster
    ccombnt.df$dist[i] <- dist( rbind (junka, junkb) ) # calculate the distance btw centroid...
    #subset matched-set 1
    #subset matched-set 2
    junkc = Data[which(Data$cluster==ccombn[1,i]),Var]
    junkd = Data[which(Data$cluster==ccombn[2,i]),Var]
    # delete NAs
    junkc<-junkc[!is.na(junkc)]
    junkd<-junkd[!is.na(junkd)]
    # calculate concordance
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
  
  ######################
  # graph between matched-set concordance
  ccombnt.df$species <- "between matched-set concordance"
  
  plot1 <- ggplot() +
    geom_point(data=ccombnt.df, aes_string(x="dist", y=Var) ,shape=19,size=1) +   
    geom_smooth(data=ccombnt.df, aes_string(x="dist", y=Var, color="species") , method="loess", size=1) +
    ylab(paste("Concordance of\n",Varlabel,sep="")) +
    xlab("Distance between matched-sets (m)") +
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
  print(plot1)
  ccombnt.df
}
####  find the knot
opfun <- function(k,Data){
  spline.lm <- lm(Fakename ~ bs(dist,degree=1, knots=c(k)),
                  data=Data)
  AIC(spline.lm)
}
####  bootstrap
lmknot <- function(optim,Var,boot.iter,Data){
  ccombnt.bs <- Data[c("dist",Var)]
  ccombnt.bsopt<-ccombnt.bs[which(ccombnt.bs$dist<optim),]
  
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
  aa<-quantile(b1,probs=c(0.025,0.975))*100
  bb<-lm(substitute(i ~ dist, list(i=as.name(Var))),data=ccombnt.bsopt)$coeff[2]*100
  print(c(aa,bb))
}


