# Script tested on 8 Oct 2022 with this configuration

# > sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.4
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#   [1] ggrepel_0.9.1     semPlot_1.1.6     lavaan_0.6-12     corrplot_0.92    
# [5] stargazer_5.2.3   progress_1.2.2    xtable_1.8-4      dplyr_1.0.10     
# [9] vars_1.5-6        lmtest_0.9-40     urca_1.3-3        strucchange_1.5-3
# [13] sandwich_3.0-2    MASS_7.3-57       reshape2_1.4.4    ggplot2_3.3.6    
# [17] pheatmap_1.0.12   ranger_0.14.1     glmnet_4.1-4      gtrendsR_1.5.1   
# [21] yuima_1.15.15     mvtnorm_1.1-3     cubature_2.0.4.4  expm_0.999-6     
# [25] data.table_1.14.2 quantmod_0.4.20   TTR_0.24.3        xts_0.12.1       
# [29] zoo_1.8-10        Matrix_1.4-1     
# 
# loaded via a namespace (and not attached):
#   [1] minqa_1.2.4         colorspace_2.0-3    deldir_1.0-6       
# [4] ellipsis_0.3.2      htmlTable_2.4.1     corpcor_1.6.10     
# [7] base64enc_0.1-3     rstudioapi_0.14     farver_2.1.1       
# [10] fansi_1.0.3         codetools_0.2-18    splines_4.2.1      
# [13] mnormt_2.1.0        knitr_1.40          glasso_1.11        
# [16] Formula_1.2-4       nloptr_2.0.3        cluster_2.1.3      
# [19] png_0.1-7           compiler_4.2.1      backports_1.4.1    
# [22] fastmap_1.1.0       cli_3.4.0           htmltools_0.5.3    
# [25] prettyunits_1.1.1   tools_4.2.1         OpenMx_2.20.6      
# [28] igraph_1.3.4        coda_0.19-4         gtable_0.3.1       
# [31] glue_1.6.2          Rcpp_1.0.9          carData_3.0-5      
# [34] vctrs_0.4.1         nlme_3.1-157        lisrelToR_0.1.5    
# [37] iterators_1.0.14    calculus_0.3.3      psych_2.2.5        
# [40] xfun_0.33           stringr_1.4.1       openxlsx_4.2.5     
# [43] lme4_1.1-30         lifecycle_1.0.2     gtools_3.9.3       
# [46] XML_3.99-0.10       scales_1.2.1        hms_1.1.2          
# [49] kutils_1.70         parallel_4.2.1      RColorBrewer_1.1-3 
# [52] curl_4.3.2          pbapply_1.5-0       gridExtra_2.3      
# [55] rpart_4.1.16        latticeExtra_0.6-30 stringi_1.7.8      
# [58] foreach_1.5.2       sem_3.1-15          checkmate_2.1.0    
# [61] boot_1.3-28         zip_2.2.1           shape_1.4.6        
# [64] rlang_1.0.5         pkgconfig_2.0.3     arm_1.12-2         
# [67] lattice_0.20-45     purrr_0.3.4         labeling_0.4.2     
# [70] htmlwidgets_1.5.4   tidyselect_1.1.2    plyr_1.8.7         
# [73] magrittr_2.0.3      R6_2.5.1            generics_0.1.3     
# [76] Hmisc_4.7-1         pillar_1.8.1        foreign_0.8-82     
# [79] withr_2.5.0         rockchalk_1.8.157   survival_3.3-1     
# [82] abind_1.4-5         nnet_7.3-17         tibble_3.1.8       
# [85] crayon_1.5.1        fdrtool_1.2.17      interp_1.1-3       
# [88] utf8_1.2.2          jpeg_0.1-9          grid_4.2.1         
# [91] qgraph_1.9.2        pbivnorm_0.6.0      digest_0.6.29      
# [94] mi_1.1              RcppParallel_5.1.5  munsell_0.5.0    


rm(list=ls())

library(xts)
library(quantmod)
library(data.table)
library(yuima)
library(gtrendsR)
library(glmnet)
library(ranger)    
library(pheatmap)  
library(ggplot2)   
library(reshape2)  
library(vars)     
library(dplyr)
library(xtable)      
library(progress)    
library(stargazer)   
library(corrplot)    
library(lavaan)      
library(semPlot)     
library(ggrepel)    

start.date <- "2019-11-01"   
start.period <- "2019-11-01/" 
full.period <- "2020-01-01/2020-10-11"



# colors
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
cutoff.distance <- 0.5  
cols <- makeColorRampPalette(c("white", "steelblue",    # distances 0 to 3 colored from white to red
                               "steelblue", "red"), # distances 3 to max(distmat) colored from green to black
                             cutoff.distance  ,
                             100)

# Download the data locally or via this command from github
# load("iData.rda")
load(url("https://github.com/siacus/swbCovid/blob/main/rda/iData.rda?raw=true"))



var.ita <- c("swbi", "iDeaths", "iCases",
             "FTSEMIB", "iCoronaVirus", "iCoronaVirusNews",
             "iCovid", "iCovidNews",  "iRt", "iUnemployment",
             "iUnemploymentNews", "iEconomy", "iEconomyNews", 
             "iGDP", "iGDPNews", "iStress", "iDepression",
             "iHealth", "iSolitude", "iInsomnia", "iResidential",
             "iWorkplace",
             "ipm25", "itemperature", "iWuhan", "iAdultContent",
             "ilockdown",
             "iFB.CLI","iFB.ILI","iFB.MC","iFB.DC","iFB.HF")


# progressive
itmpL <- cbind(itmp, stats::lag(itmp[,"swbi"],1)) 
colnames(itmpL)[ncol(itmpL)] <- "swbiLag"
var.itaL <- c(var.ita,"swbiLag")

dati <- data.frame(na.approx(itmpL["/2020-10-11",var.itaL],rule=2))
dati <- data.frame( scale( dati ))

iperiods <- c("2019-11",
              "2019-12",
              "2020-01",
              "2020-02",
              "2020-03",
              "2020-04",
              "2020-05",
              "2020-06",
              "2020-07",
              "2020-08",
              "2020-09",
              "2020-01/2020-09"
)

incc <- length(var.itaL)-1
inpp <- length(iperiods)
iPV <- matrix(NA, inpp,incc)
iCC <- matrix(NA, inpp,incc)
ik <- 0
for(iper in iperiods){
  ik <- ik+1
  iSub <- data.frame(na.approx(itmpL[iper,var.itaL],rule=2))
  for(i in 1:ncol(iSub)){
    ux <- which(is.na(iSub[,i]))
    if(length(ux)>0){
      iSub[ux,i] <- 0
    }
  } 
  for(j in 1:incc){
    itmpc <- na.omit(cbind(iSub[,1],iSub[,j+1]))
    ict <- cor.test(itmpc[,1],itmpc[,2],method = "spearman")
    iCC[ik,j] <- ict$estimate
    iPV[ik,j] <- ict$p.value
  }
}
iCC[is.na(iCC)] <- 0
iPV[is.na(iPV)] <- 1
iCC[iPV>0.05] <- 0


rownames(iCC) <- rownames(iPV) <- 
  c("Nov 2019", "Dec 2019", "Jan 2020", "Feb 2020",
                  "Mar 2020", "Apr 2020", "Mag 2020", "Jun 2020",
                  "Jul 2020", "Aug 2020", "Sep 2020", "Year 2020")
colnames(iCC) <- colnames(iPV) <- colnames(iSub)[-1]
corcol <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
iCr <- iCC[,-grep("swb",colnames(iCC))]

icol.order <- try(hclust(dist(t(iCr)))$order, TRUE)

# Monthly correlation plots (IT)
corrplot(iCr[,icol.order],col=corcol(100),tl.col="black", tl.srt=45)
pdf("corItaly.pdf",width=9,height=6)
corrplot(iCr[,icol.order],col=corcol(100),tl.col="black", tl.srt=45)
dev.off()

# DynENET (IT)
idates <- rownames(dati)
lambda <- seq(1e-6,1,length=250)
imsem <- NULL
imse1 <- NULL
imsear <- NULL
intercepti <- FALSE
#for(ialpha in c(0,0.5,1)){
for(ialpha in c(0.5)){
  icf1 <- NULL
  icfm <- NULL
  ircf1 <- NULL
  ircfm <- NULL
  ip1 <- NULL
  ipm <- NULL
  il1 <- NULL
  ilm <- NULL
  iar <- NULL
  ni <- nrow(dati)
  ioffset <- 31
  pb <- progress_bar$new(total = ni-ioffset)
  for(i in ioffset:(ni-1)){
    pb$tick()
    rg <- (i-ioffset+1):i
    y <- dati$swbi[rg]
    tmp.iar <- try(as.numeric(predict(arima(y,order=c(1,0,1)),n.ahead = 1)$pred), TRUE)
    if(class(tmp.iar)[1]!="try-error"){
      iar <- c(iar,tmp.iar)
    } else {
      iar <- c(iar,as.numeric(NA))
    }
    x <- as.matrix(dati[rg,colnames(dati)!="swbi"])
    x.new <- as.matrix(dati[i+1,colnames(dati)!="swbi"],nrow=1)
    set.seed(123)
    cvi <- cv.glmnet(x=x, y=y, alpha=ialpha, intercept=intercepti, lambda=lambda, grouped=FALSE)
    set.seed(123)
    inet <- glmnet(x=x, y=y, alpha=ialpha,  intercept=intercepti, lambda=lambda)
    icf.1se <- coef(inet, s=cvi$lambda.1se)
    icf.min <- coef(inet, s=cvi$lambda.min)
    t1 <- data.frame(date=idates[i], t(as.matrix(icf.1se)))
    colnames(t1)[2] <- "Intercept"
    t2 <- data.frame(date=idates[i], t(as.matrix(icf.min)))
    colnames(t2)[2] <- "Intercept"
    icf1 <- rbind(icf1, t1)
    icfm <- rbind(icfm, t2)
    ip1 <- c(ip1, predict(inet, newx = x.new,s = cvi$lambda.1se))
    ipm <- c(ipm, predict(inet, newx = x.new,s = cvi$lambda.min))
    il1 <- c(il1, cvi$lambda.1se)
    ilm <- c(ilm, cvi$lambda.min)

    vv <- c("swbi",colnames(x)[which(as.numeric(icf.1se)!=0) -1])
    if(length(vv)>1){
      set.seed(123)
      rang <- ranger(swbi~., data=dati[rg,vv], importance="impurity",num.trees = 2500)
      imp <- sort(importance(rang), decreasing = TRUE)
      rank <- (length(imp):1)/length(imp)
      names(rank) <- names(imp)
      tmprank <- t1[,-2]
      tmprank[1,-1] <- 0
      tmprank[1, names(rank)] <- rank
      ircf1 <- rbind(ircf1, tmprank)
    }
  
    vv <- c("swbi",colnames(x)[which(as.numeric(icf.min)!=0) -1])
    if(length(vv)>1){
      set.seed(123)
      rang <- ranger(swbi~., data=dati[rg,vv], importance="impurity")
      imp <- sort(importance(rang), decreasing = TRUE)
      rank <- (length(imp):1)/length(imp)
      names(rank) <- names(imp)
      tmprank <- t2[,-2]
      tmprank[1,-1] <- 0
      tmprank[1, names(rank)] <- rank
      ircfm <- rbind(ircfm, tmprank)
    }
  }
  ipt <- dati$swbi[(ioffset+1):ni]
  id <- data.frame(true=ipt,p1=ip1,pm=ipm,ar=iar)
  imse1 <- c(imse1, mean((ipt-ip1)^2))
  imsem <- c(imsem, mean((ipt-ipm)^2))
  imsear <- c(imsear, mean((ipt-iar)^2,na.rm=TRUE))
  print(cor(id,use="pairwise"))
}
imse1
imsem
imsear
imse1/imsear

dati <- data.frame(na.approx(itmpL["/2020-10-11",var.itaL],rule=2))
dati <- data.frame( scale( dati ))
iyear <- year(as.Date(rownames(dati)))
imonth <- month(as.Date(rownames(dati)))
imod <- list()
k <- 0
for(iy in 2020){
  for(im in 1:9){
    idx <- which(iyear==iy & imonth==im)
    if(length(idx)>0){
      print(c(iy,im))
      k <- k+1
      isub <- dati[idx,]
      mod <- step(lm(swbi~ -1 + ., data=isub))
      imod[[k]] <- list(year=iy,mod=mod, month=im)
    }
  }
}

imodFull <- step(lm(swbi~ -1 + ., data=dati[iyear==2020 & imonth<=10,]))

mi1 <- imodFull 
mi2 <- imod[[1]]$mod
mi3 <- imod[[2]]$mod
mi4 <- imod[[3]]$mod
mi5 <- imod[[4]]$mod
mi6 <- imod[[5]]$mod
mi7 <- imod[[6]]$mod
mi8 <- imod[[7]]$mod
mi9 <- imod[[8]]$mod

IM <- stargazer::stargazer(
    mi1,mi2,mi3,mi4,mi5,
    mi6,mi7,mi8,mi9,
    no.space=TRUE, align=TRUE,
    omit.stat=c("LL","ser","f"),
    column.labels=c("Jan-Sep",
                    "Jan", "Feb", "Mar", "Apr", "May",
                    "Jun", "Jul", "Aug", "Sep"
    ),
    model.numbers=FALSE,digits = 2, dep.var.labels="SWB-I")
  

# Benchmarking DynENET with ARIMA(1,0,1) - (IT)
itoplot <- xts(id[,c("true","p1","ar")],order.by=as.Date(idates[(ioffset+1):ni]))
colnames(itoplot) <- c("SWB-I","Elastic Net","ARIMA")

pdf("iforecast.pdf")
autoplot(itoplot, facets = NULL) +xlab("")
dev.off()

ixcfm <- xts(icfm[,-1], order.by=as.Date(icfm[,1]))
ixcf1 <- xts(icf1[,-1], order.by=as.Date(icf1[,1]))


iA <- icf1[,var.ita[-1]]
iA <- apply(iA,2,function(u) u!=0)
iB <- ircf1[,var.ita[-1]]
inn <- ncol(iA) 
isumm <- NULL
for(i in 1:inn){
  iavg <- sum(iA[,i])
  iravg <- mean(iB[,i])
  iravgpos <- mean(iB[iB[,i]>0,i])
  isumm <- rbind(isumm, cbind(iavg,iravg,iravgpos))
}
colnames(isumm) <- c("weeks","rank","posrank")
rownames(isumm) <- colnames(iA)
isumm <- data.frame(isumm)
isumm$var <- colnames(iA)
isumm$country <- "Italy"
save(isumm,file="isumm.rda")

# IF Space (IT)
pdf(file="isummary.pdf",width=9,height=6)
ggplot(data=isumm,aes(x=weeks,y=rank,label=var)) + 
  geom_point() + 
  geom_text_repel(  point.padding = unit(0.35, "lines")) +
  xlab(sprintf("number of times variable is selected over %d analyses",nrow(iA)))+
  ylab("average relative rank")
dev.off()  

pdf(file="isummaryPos.pdf",width=9,height=6)
ggplot(data=isumm,aes(x=weeks,y=posrank,label=var)) + 
  geom_point() + 
  geom_text_repel(  point.padding = unit(0.35, "lines")) +
  xlab(sprintf("number of times variable is selected over %d analyses",nrow(iA)))+
  ylab("average relative rank")
dev.off()  



bigmat <- ircf1[,-1]   
rownames(bigmat) <- ircfm$date
bigmat <- bigmat[,-match("swbiLag",colnames(bigmat))]
bigmat[is.na(bigmat) | bigmat==0] <- 0
dp <- colSums(bigmat)
dp
dim(bigmat)
bigmat <- t(bigmat)
gtrend <- c("iCoronaVirus"  ,   "iCoronaVirusNews", 
               "iCovid", "iCovidNews",  "iRt",
            "iUnemployment", "iUnemploymentNews",
            "iEconomy", "iEconomyNews",
            "iGDP","iGDPNews",
            "iStress","iDepression","iHealth",
            "iSolitude","iInsomnia","iWuhan","iAdultContent")
pandemic <- c("iDeaths","iCases")
finance <- c("FTSEMIB")
mobility <- c("iResidential","iWorkplace")
airquality <- c("ipm25","itemperature")
policy <- ("ilockdown") 
lagSWBI <- ("swbiLag") 
facebook <- c("iFB.CLI","iFB.ILI","iFB.MC","iFB.DC","iFB.HF")
group <- rep("", nrow(bigmat))
group[match(pandemic, rownames(bigmat))] <- "Pandemic"
group[match(gtrend, rownames(bigmat))] <- "GTrend"
group[match(finance, rownames(bigmat))] <- "Finance"
group[match(mobility, rownames(bigmat))] <- "Mobility"
group[match(airquality, rownames(bigmat))] <- "Environment"
group[match(policy, rownames(bigmat))] <- "Policy"
group[match(lagSWBI, rownames(bigmat))] <- "AR effect"
group[match(facebook, rownames(bigmat))] <- "Facebook"

annotation_row = data.frame( DataClass = factor(group) )
rownames(annotation_row) = rownames(bigmat)
row.order <- try(hclust(dist(bigmat))$order, TRUE)
dat_new <- bigmat[row.order, ] # re-order matrix accoring to clustering

nba.m <- reshape2::melt(dat_new)
mycol <- colorRampPalette(c("white", "steelblue", "red"))(12)
nba.m$Var2 <- as.Date(nba.m$Var2)

bb <- ggplot(nba.m, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value), colour = "gray50") + 
  scale_fill_gradient2(midpoint=0.5, low = "white", mid = "steelblue",high = "red")+
  labs(x = "", y = "Selected variables", fill="Rank")+
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  #scale_x_date(date_minor_breaks = "1 month")
  scale_x_date(date_breaks = "1 week",  
               date_labels = "%d %b", 
               limits = as.Date(range(colnames(bigmat))))+
  ggtitle("SWB-Italy")
bb

pdf(sprintf("importance-%s.pdf","Italy"),width=10, height=5)
print(bb)
dev.off()

clabels <- colnames(bigmat)
ipos <- unique(c(seq(1,length(clabels),by=7),length(clabels)))

clabels[-ipos] <- ""
pheatmap(bigmat, cluster_cols = FALSE,
             #  clustering_distance_rows =  "correlation",
            #  color = c("white", "steelblue", "red"),
             color=cols,
             # legend = FALSE,
         cellwidth=3,
             labels_col = clabels,
             annotation_row = annotation_row,
             treeheight_row = 0, treeheight_col = 0,
             fontsize = 12,
             fontsize_col = 12,
             border_color = "gray",
             main = sprintf("%s : relative importance of predictors selected by Dynamic Elastic Net","Italy"),
             angle_col=90,filename = sprintf("importance2-%s.pdf","Italy"),
            width=17,height =10
         )

bigmat <- icf1[,-(1:2)]
rownames(bigmat) <- icfm$date

bigmat <- bigmat[,-match("swbiLag",colnames(bigmat))]

icols <- makeColorRampPalette(c("red", "white",    # distances 0 to 3 colored from white to red
                                "white", "green"), # distances 3 to max(distmat) colored from green to black
                              1-max(bigmat)/diff(range(bigmat))  ,
                              100)

bigmat[is.na(bigmat) | bigmat==0] <- 0
dp <- colSums(bigmat)
dp
dim(bigmat)
bigmat <- t(bigmat)
group <- rep("", nrow(bigmat))
group[match(pandemic, rownames(bigmat))] <- "Pandemic"
group[match(gtrend, rownames(bigmat))] <- "GTrend"
group[match(finance, rownames(bigmat))] <- "Finance"
group[match(mobility, rownames(bigmat))] <- "Mobility"
group[match(airquality, rownames(bigmat))] <- "Air Quality"
group[match(policy, rownames(bigmat))] <- "Policy"
group[match(lagSWBI, rownames(bigmat))] <- "AR effect"
group[match(facebook, rownames(bigmat))] <- "Facebook"

annotation_row = data.frame( DataClass = factor(group) )
rownames(annotation_row) = rownames(bigmat)
row.order <- try(hclust(dist(bigmat))$order, TRUE)
dat_new <- bigmat[row.order, ] # re-order matrix accoring to clustering

nba.m <- reshape2::melt(dat_new)

mycol <- colorRampPalette(c("red", "white", "blue"))(7)

Mx <- max(range(bigmat))*1.05
mx <- min(range(bigmat))*1.05
sq <- c(seq(mx,-1e-5, length= 4),  0,seq(1e-5,Mx,length=4))
nba.m$cat <- cut(nba.m$value,breaks=sq)
lv <- levels(nba.m$cat)
nba.m$cat <- as.character(nba.m$cat)
nba.m$cat[nba.m$value==0] <- "0"
nba.m$cat[nba.m$cat == "(-1e-05,0]"] <- lv[3]
nba.m$cat[nba.m$cat == "(0,1e-05]"] <- lv[6]
nba.m$cat <- sub("-1e-05","0",nba.m$cat)
nba.m$cat <- sub("1e-05","0",nba.m$cat)
lv2 <- lv
lv2 <- sub("-1e-05","0",lv2)
lv2 <- sub("1e-05","0",lv2)
lv2 <- c(lv2[1:3],"0",lv2[-(1:5)])

nba.m$cat2 <- factor(nba.m$cat, levels = lv2, ordered = TRUE )
table(nba.m$cat2)
nba.m$Var2 <- as.Date(nba.m$Var2)
hclust(dist(bigmat)) -> oo
nba.m$Var1 <- factor(nba.m$Var1, 
                     levels = rev(oo$labels[oo$order]),ordered = TRUE)

aa <- ggplot(nba.m, aes(Var2, Var1)) + 
  geom_tile( aes(fill = cat2),
             colour = "darkgray") + 
  scale_fill_manual(values=mycol, 
                    breaks=levels(nba.m$cat2) )+
  labs(x = "", y = "Variable", fill="Size of\nstandardized\ncoefficents",
       title="SWB-Italy")+
  scale_x_date(date_breaks = "1 week",  
               date_labels = "%d %b", 
               limits = as.Date(range(colnames(bigmat))))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
aa
pdf("Coefficients-Italy.pdf",width=17, height=10, pointsize = 12)
print(aa)
dev.off()




######### JAPAN  

# Download the data locally or via this command from github
# load("jData.rda")
load(url("https://github.com/siacus/swbCovid/blob/main/rda/jData.rda?raw=true"))

var.jap <- c("swbj", "jDeaths", "jCases", "NIKKEI",
             "jCoronaVirus", "jCoronaVirusNews",
             "jCorona", "jCoronaNews", "jCovid",
             "jCovidNews", "jUnemployment", "jUnemploymentNews",
             "jEconomy", "jEconomyNews", "jGDP", "jGDPNews",
             "jStress", "jDepression", "jHealth", "jSolitude",
             "jInsomnia", "jResidential", "jWorkplace",
             "jpm25", "jtemperature",
             "jWuhan", "jAdultContent", "jlockdown",
             "jFB.CLI", "jFB.ILI", "jFB.MC", "jFB.DC", "jFB.HF"
             ) 
# progressive
jtmpL <- cbind(jtmp, stats::lag(jtmp[,"swbj"],1))
colnames(jtmpL)[ncol(jtmpL)] <- "swbjLag"
var.japL <- c(var.jap,"swbjLag")



datj <- data.frame(na.approx(jtmpL["/2020-09-20",var.japL],rule=2))
datj <- data.frame( scale(datj) )

jperiods <- c("2019-11",
              "2019-12",
              "2020-01",
              "2020-02",
              "2020-03",
              "2020-04",
              "2020-05",
              "2020-06",
              "2020-07",
              "2020-08",
              "2020-09",
              "2020-01/2020-09"
)

jncc <- length(var.japL)-1
jnpp <- length(jperiods)
jPV <- matrix(NA, jnpp,jncc)
jCC <- matrix(NA, jnpp,jncc)
jk <- 0
for(jper in jperiods){
  jk <- jk+1
  jSub <- data.frame(na.approx(jtmpL[jper,var.japL],rule=2))
  for(i in 1:ncol(jSub)){
    ux <- which(is.na(jSub[,i]))
    if(length(ux)>0){
      jSub[ux,i] <- 0
    }
  } 
  for(j in 1:jncc){
    jtmpc <- na.omit(cbind(jSub[,1],jSub[,j+1]))
    jct <- cor.test(jtmpc[,1],jtmpc[,2],method = "spearman")
    jCC[jk,j] <- jct$estimate
    jPV[jk,j] <- jct$p.value
  }
}
jCC[is.na(jCC)] <- 0
jPV[is.na(jPV)] <- 1
jCC[jPV>0.05] <- 0


rownames(jCC) <- rownames(jPV) <- 
  c("Nov 2019", "Dec 2019", "Jan 2020", "Feb 2020",
    "Mar 2020", "Apr 2020", "Mag 2020", "Jun 2020",
    "Jul 2020", "Aug 2020", "Sep 2020", "Year 2020")
colnames(jCC) <- colnames(jPV) <- colnames(jSub)[-1]
corcol <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
jCr <- jCC[,-grep("swb",colnames(jCC))]

jcol.order <- try(hclust(dist(t(jCr)))$order, TRUE)

# Monthly correlation plot (JP)
corrplot(jCr[,jcol.order],col=corcol(100),tl.col="black", tl.srt=45)
pdf("corJapan.pdf",width=9,height=6)
corrplot(jCr[,jcol.order],col=corcol(100),tl.col="black", tl.srt=45)
dev.off()



# DynENET (JP)

jdates <- rownames(datj)
lambda <- seq(1e-6,1,length=250)

jmsem <- NULL
jmse1 <- NULL
jmsear <- NULL
interceptj <- FALSE
#for(jalpha in c(0,0.5,1)){
for(jalpha in c(0.5)){
  jcf1 <- NULL
  jcfm <- NULL
  jrcf1 <- NULL
  jrcfm <- NULL
  jp1 <- NULL
  jpm <- NULL
  jl1 <- NULL
  jlm <- NULL
  jar <- NULL
  nj <- nrow(datj)
  joffset <- 31
  pb <- progress_bar$new(total = nj-joffset)
  for(i in joffset:(nj-1)){
    pb$tick()
    rg <- (i-joffset+1):i
    y <- datj$swbj[rg]
    tmp.jar <- try(as.numeric(predict(arima(y,order=c(1,0,1)),n.ahead = 1)$pred), TRUE)
    if(class(tmp.jar)[1]!="try-error"){
      jar <- c(jar,tmp.jar)
    } else {
      jar <- c(jar, as.numeric(NA))
    }
    x <- as.matrix(datj[rg,colnames(datj)!="swbj"])
    x.new <- as.matrix(datj[i+1,colnames(datj)!="swbj"],nrow=1)
    set.seed(123)
    cvj <- cv.glmnet(x=x, y=y, alpha=jalpha,intercept=interceptj,nfolds=10,lambda=lambda,grouped=FALSE)
    set.seed(123)
    jnet <- glmnet(x=x, y=y, alpha=jalpha,intercept=interceptj,lambda=lambda)
    jcf.1se <- coef(jnet, s=cvj$lambda.1se)
    jcf.min <- coef(jnet, s=cvj$lambda.min)
    t1 <- data.frame(date=jdates[i], t(as.matrix(jcf.1se)))
    colnames(t1)[2] <- "Intercept"
    t2 <- data.frame(date=jdates[i], t(as.matrix(jcf.min)))
    colnames(t2)[2] <- "Intercept"
    jcf1 <- rbind(jcf1, t1)
    jcfm <- rbind(jcfm, t2)
    jp1 <- c(jp1, predict(jnet, newx = x.new,s = cvj$lambda.1se))
    jpm <- c(jpm, predict(jnet, newx = x.new,s = cvj$lambda.min))
    jl1 <- c(jl1, cvj$lambda.1se)
    jlm <- c(jlm, cvj$lambda.min)
  
    vv <- c("swbj",colnames(x)[which(as.numeric(jcf.1se)!=0) -1])
    if(length(vv)>1){
      set.seed(123)
      rang <- ranger(swbj~., data=datj[rg,vv], importance="impurity",num.trees = 2500)
      imp <- sort(importance(rang), decreasing = TRUE)
      rank <- (length(imp):1)/length(imp)
      names(rank) <- names(imp)
      tmprank <- t1[,-2]
      tmprank[1,-1] <- 0
      tmprank[1, names(rank)] <- rank
      jrcf1 <- rbind(jrcf1, tmprank)
    }
  
    vv <- c("swbj",colnames(x)[which(as.numeric(jcf.min)!=0) -1])
    if(length(vv)>1){
      set.seed(123)
      rang <- ranger(swbj~., data=datj[rg,vv], importance="impurity")
      imp <- sort(importance(rang), decreasing = TRUE)
      rank <- (length(imp):1)/length(imp)
      names(rank) <- names(imp)
      tmprank <- t2[,-2]
      tmprank[1,-1] <- 0
      tmprank[1, names(rank)] <- rank
      jrcfm <- rbind(jrcfm, tmprank)
    }
  }
  jpt <- datj$swbj[(joffset+1):nj]
  jd <- data.frame(true=jpt,p1=jp1,pm=jpm,ar=jar)
  jmse1 <- c(jmse1, mean((jpt-jp1)^2))
  jmsem <- c(jmsem, mean((jpt-jpm)^2))
  jmsear <- c(jmsear, mean((jpt-jar)^2,na.rm=TRUE))
  print(cor(jd,use="pairwise"))
}
jmse1
jmsem
jmsear
jmse1/jmsear



jA <- jcf1[,var.jap[-1]]
jA <- apply(jA,2,function(u) u!=0)
jB <- jrcf1[,var.jap[-1]]
jnn <- ncol(jA) 
jsumm <- NULL
for(j in 1:jnn){
  javg <- sum(jA[,j])
  jravg <- mean(jB[,j])
  jravgpos <- mean(jB[jB[,j]>0,j])
  jsumm <- rbind(jsumm, cbind(javg,jravg,jravgpos))
}
colnames(jsumm) <- c("weeks","rank","posrank")
rownames(jsumm) <- colnames(jA)
jsumm <- data.frame(jsumm)
jsumm$var <- colnames(jA)
jsumm$country <- "Japan"
save(jsumm,file="jsumm.rda")


# IF Space (JP)
pdf(file="jsummary.pdf",width=9,height=6)
ggplot(data=jsumm,aes(x=weeks,y=rank,label=var)) + 
  geom_point() + 
  geom_text_repel(  point.padding = unit(0.35, "lines")) +
  xlab(sprintf("number of times variable is selected over %d analyses",nrow(jA)))+
  ylab("average relative rank")
dev.off()  

pdf(file="jsummaryPos.pdf",width=9,height=6)
ggplot(data=jsumm,aes(x=weeks,y=posrank,label=var)) + 
  geom_point() + 
  geom_text_repel(  point.padding = unit(0.35, "lines"), show_guide=F) +
  xlab(sprintf("number of times variable is selected over %d analyses",nrow(jA)))+
  ylab("average relative rank")+
  labs(color="Country")
dev.off() 




rbind(isumm,jsumm) -> summ

pdf(file="summary.pdf",width=9,height=6)
ggplot(data=summ,aes(x=weeks,y=rank,label=var,color=country)) +
  geom_point() +
  geom_text_repel(  point.padding = unit(0.35, "lines")) +
  xlab("number of times variable is selected")+
  ylab("average relative rank")
dev.off()

pdf(file="summaryPos.pdf",width=9,height=6)
ggplot(data=summ,aes(x=weeks,y=posrank,label=var,color=country)) +
  geom_point() +
  geom_text_repel(  point.padding = unit(0.35, "lines")) +
  xlab("number of times variable is selected over")+
  ylab("average relative rank")
dev.off()


pdf(file="summary.pdf",width=10,height=7)
ggplot(data=summ,aes(x=weeks,y=rank,label=var,color=country)) + 
     geom_point() + 
     geom_text_repel(  point.padding = unit(0.75, "lines"), show_guide=F) +
     xlab("number of times variable is selected")+
     ylab("average relative rank")+
     labs(color="Country")

    
dev.off()

pdf(file="summaryPos.pdf",width=10,height=7)
ggplot(data=summ,aes(x=weeks,y=posrank,label=var,color=country)) + 
     geom_point() + 
     geom_text_repel(  point.padding = unit(0.35, "lines")) +
     xlab("number of times variable is selected over")+
     ylab("average relative rank")

dev.off()

# Benchmarking DynENET against ARIMA(1,0,1) - (JP)
datj <- data.frame(na.approx(jtmpL["/2020-09-20",var.japL],rule=2))
datj <- data.frame( scale(datj) )
jyear <- year(as.Date(rownames(datj)))
jmonth <- month(as.Date(rownames(datj)))
jmod <- list()
k <- 0
for(jy in 2020){
  for(jm in 1:10){
    jdx <- which(jyear==jy & jmonth==jm)
    if(length(jdx)>0){
      print(c(jy,jm))
      k <- k+1
      jsub <- datj[jdx,]
      mod <- step(lm(swbj~ -1 + ., data=jsub))
      jmod[[k]] <- list(year=jy,mod=mod, month=jm)
    }
  }
}

jmodFull <- step(lm(swbj~ -1 + ., data=datj[jyear==2020 & jmonth<=10,]))


mj1 <- jmodFull 
mj2 <- jmod[[1]]$mod
mj3 <- jmod[[2]]$mod
mj4 <- jmod[[3]]$mod
mj5 <- jmod[[4]]$mod
mj6 <- jmod[[5]]$mod
mj7 <- jmod[[6]]$mod
mj8 <- jmod[[7]]$mod
mj9 <- jmod[[8]]$mod

JM <- stargazer::stargazer(
  mj1,mj2,mj3,mj4,mj5,
  mj6,mj7,mj8,mj9,
  no.space=TRUE, align=TRUE,
  omit.stat=c("LL","ser","f"),
  column.labels=c("Jan-Sep",
                  "Jan", "Feb", "Mar", "Apr", "May",
                  "Jun", "Jul", "Aug", "Sep"
  ),
  model.numbers=FALSE,digits = 2, dep.var.labels="SWB-J")



jtoplot <- xts(jd[,c("true","p1","ar")],order.by=as.Date(jdates[(joffset+1):nj]))
colnames(jtoplot) <- c("SWB-J","Elastic Net","ARIMA")

pdf("jforecast.pdf")
autoplot(jtoplot, facets = NULL) +xlab("")
dev.off()

jxcfm <- xts(jcfm[,-1], order.by=as.Date(jcfm[,1]))
jxcf1 <- xts(jcf1[,-1], order.by=as.Date(jcf1[,1]))


bigmat <- jrcf1[,-1]  


rownames(bigmat) <- jrcfm$date
bigmat <- bigmat[,-match("swbjLag",colnames(bigmat))]


bigmat[is.na(bigmat) | bigmat==0] <- 0
dp <- colSums(bigmat)
dp
dim(bigmat)
bigmat <- t(bigmat)
gtrend <- c("jCoronaVirus"  ,   "jCoronaVirusNews", 
            "jCorona",          "jCoronaNews" ,     
            "jCovid", "jCovidNews",  
            "jUnemployment","jUnemploymentNews",
            "jEconomy","jEconomyNews",
            "jGDP","jGDPNews",
            "jStress","jDepression","jHealth",
            "jSolitude","jInsomnia","jWuhan",
            "jAdultContent")
pandemic <- c("jDeaths","jCases")
finance <- c("NIKKEI")
mobility <- c("jResidential","jWorkplace")
airquality <- c("jpm25","jtemperature")
policy <- ("jlockdown") 
lagSWBJ <- ("swbjLag") 
facebook <- c("jFB.CLI","jFB.ILI","jFB.MC","jFB.DC","jFB.HF")

group <- rep("", nrow(bigmat))
group[match(pandemic, rownames(bigmat))] <- "Pandemic"
group[match(gtrend, rownames(bigmat))] <- "GTrend"
group[match(finance, rownames(bigmat))] <- "Finance"
group[match(mobility, rownames(bigmat))] <- "Mobility"
group[match(airquality, rownames(bigmat))] <- "Environment"
group[match(policy, rownames(bigmat))] <- "Policy"
group[match(lagSWBJ, rownames(bigmat))] <- "AR effect"
group[match(facebook, rownames(bigmat))] <- "Facebook"

annotation_row = data.frame( DataClass = factor(group) )
rownames(annotation_row) = rownames(bigmat)
row.order <- try(hclust(dist(bigmat))$order, TRUE)
dat_new <- bigmat[row.order, ] # re-order matrix accoring to clustering

nba.m <- reshape2::melt(dat_new)
mycol <- colorRampPalette(c("white", "steelblue", "red"))(12)
nba.m$Var2 <- as.Date(nba.m$Var2)


bb <- ggplot(nba.m, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value), colour = "gray50") + 
  scale_fill_gradient2(midpoint=0.5, low = "white", mid = "steelblue",high = "red")+
  labs(x = "", y = "Selected variables", fill="Rank")+
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  #scale_x_date(date_minor_breaks = "1 month")
  scale_x_date(date_breaks = "1 week",  
               date_labels = "%d %b", 
               limits = as.Date(range(colnames(bigmat))))+
  ggtitle("SWB-Japan")
bb

pdf(sprintf("importance-%s.pdf","Japan"),width=9, height=6)
print(bb)
dev.off()

clabels <- colnames(bigmat)
jpos <- unique(c(seq(1,length(clabels),by=7),length(clabels)))
clabels[-jpos] <- ""

pheatmap(bigmat, cluster_cols = FALSE,
         #  clustering_distance_rows =  "correlation",
         #  color = c("white", "steelblue", "red"),
         color=mycol,
         # legend = FALSE,
         cellwidth=3,
         labels_col = clabels,
         annotation_row = annotation_row,
         treeheight_row = 0, treeheight_col = 0,
         fontsize = 12,
         fontsize_col = 12,
         border_color = "gray",
         main = sprintf("%s : relative importance of predictors selected by Dynamic Elastic Net","Japan"),
         angle_col="90",filename = sprintf("importance2-%s.pdf","Japan"),
         width=17,height =10
)
 


bigmat <- jcf1[,-(1:2)]
rownames(bigmat) <- jcfm$date
bigmat <- bigmat[,-match("swbjLag",colnames(bigmat))]

jratio <- max(bigmat)/diff(range(bigmat))
                           
makeColorRampPalette2 <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}

jcols <- makeColorRampPalette2(c("darkred", "red",    # distances 0 to 3 colored from white to red
                                "lightgreen", "green"), # distances 3 to max(distmat) colored from green to black
                              1-jratio,
                              20)

bigmat[is.na(bigmat) | bigmat==0] <- 0
dp <- colSums(bigmat)
dp
dim(bigmat)
bigmat <- t(bigmat)

group <- rep("", nrow(bigmat))
group[match(pandemic, rownames(bigmat))] <- "Pandemic"
group[match(gtrend, rownames(bigmat))] <- "GTrend"
group[match(finance, rownames(bigmat))] <- "Finance"
group[match(mobility, rownames(bigmat))] <- "Mobility"
group[match(airquality, rownames(bigmat))] <- "Air Quality"
group[match(policy, rownames(bigmat))] <- "Policy"
group[match(lagSWBJ, rownames(bigmat))] <- "AR effect"
group[match(facebook, rownames(bigmat))] <- "Facebook"

annotation_row = data.frame( DataClass = factor(group) )
rownames(annotation_row) = rownames(bigmat)
row.order <- try(hclust(dist(bigmat))$order, TRUE)
dat_new <- bigmat[row.order, ] # re-order matrix accoring to clustering


nba.m <- reshape2::melt(dat_new)
mycol <- colorRampPalette(c("red", "white", "blue"))(7)

Mx <- max(range(bigmat))*1.05
mx <- min(range(bigmat))*1.05
sq <- c(seq(mx,-1e-5, length= 4),  0,seq(1e-5,Mx,length=4))
nba.m$cat <- cut(nba.m$value,breaks=sq)
lv <- levels(nba.m$cat)
nba.m$cat <- as.character(nba.m$cat)
nba.m$cat[nba.m$value==0] <- "0"
nba.m$cat[nba.m$cat == "(-1e-05,0]"] <- lv[3]
nba.m$cat[nba.m$cat == "(0,1e-05]"] <- lv[6]
nba.m$cat <- sub("-1e-05","0",nba.m$cat)
nba.m$cat <- sub("1e-05","0",nba.m$cat)
lv2 <- lv
lv2 <- sub("-1e-05","0",lv2)
lv2 <- sub("1e-05","0",lv2)
lv2 <- c(lv2[1:3],"0",lv2[-(1:5)])

nba.m$cat2 <- factor(nba.m$cat, levels = lv2, ordered = TRUE )
table(nba.m$cat2)
nba.m$Var2 <- as.Date(nba.m$Var2)
hclust(dist(bigmat)) -> oo
nba.m$Var1 <- factor(nba.m$Var1, 
                     levels = rev(oo$labels[oo$order]),ordered = TRUE)

aa <- ggplot(nba.m, aes(Var2, Var1)) + 
  geom_tile( aes(fill = cat2),
             colour = "darkgray") + 
  scale_fill_manual(values=mycol, 
                    breaks=levels(nba.m$cat2) )+
  labs(x = "", y = "Variable", fill="Size of\nstandardized\ncoefficents",
       title="SWB-Japan")+
  scale_x_date(date_breaks = "1 week",  
               date_labels = "%d %b", 
               limits = as.Date(range(colnames(bigmat))))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
aa
pdf("Coefficients-Japan.pdf",width=17, height=10,pointsize = 12)
print(aa)
dev.off()


# SEM

# Download the data locally or via this command from github
# load("iData.rda")
load(url("https://github.com/siacus/swbCovid/blob/main/rda/iData.rda?raw=true"))

var.ita <- c("swbi", "iDeaths", "iCases",
             "FTSEMIB", "iCoronaVirus", "iCoronaVirusNews",
             "iCovid", "iCovidNews",  "iRt", "iUnemployment",
             "iUnemploymentNews", "iEconomy", "iEconomyNews", 
             "iGDP", "iGDPNews", "iStress", "iDepression",
             "iHealth", "iSolitude", "iInsomnia", "iResidential",
             "iWorkplace",
             "ipm25", "itemperature", "iWuhan", "iAdultContent",
             "ilockdown",
             "iFB.CLI","iFB.ILI","iFB.MC","iFB.DC","iFB.HF")



dati <- data.frame(na.approx(itmp["2020-01-01/2020-09-30",var.ita],rule=2))
dati <- data.frame( scale( dati ))

ita_formula ="
VirusSearch =~   iCoronaVirus + iCovid +iRt    
PsySearch =~  iStress    + iInsomnia    + iSolitude +iDepression
HealthStatus =~ iFB.CLI + iFB.ILI  
Mobility =~ iResidential + iWorkplace + ilockdown 
Finance =~ FTSEMIB + iFB.HF + iUnemployment 
SocDist =~  iFB.MC + iFB.DC 
WellBeing =~ swbi
WellBeing ~ VirusSearch + HealthStatus + Mobility + Finance + SocDist
Mobility ~~ PsySearch + SocDist +iCases
PsySearch  ~ WellBeing  
iAdultContent ~ WellBeing
"

sem.ita <- sem(ita_formula, dati, std.lv = TRUE)
pars.ita <- summary(sem.ita, fit.measures=FALSE)



pars.ita$PE$exo <- NULL 


mymod <- sem.ita

pdf("sem-ita-paper.pdf",width = 20,height=15)
semPaths(mymod, style ="lisrel", 
         whatLabels = "est", edge.label.cex = .5,  
         #nodeLabels = mylab,
         nCharNodes = 20,sizeMan2=2,
         nDigits=2, #levels=c(-20,-10,0,10,20),
         label.prop=0.7, edge.label.color = "black", rotation = 2, 
         equalizeManifests = TRUE, optimizeLatRes = TRUE, node.width = 1.5, 
         edge.width = 0.6, shapeMan = "rectangle", shapeLat =  "ellipse", 
         shapeInt = "triangle", sizeMan = 4, sizeInt = 4, sizeLat = 6, 
         curve=1.5, unCol = "#070b8c",layout = "tree",
         residuals=FALSE, exoCov = TRUE, covAtResiduals = FALSE,springLevels=TRUE,
         exoVar=FALSE, intercepts=FALSE)
dev.off()



# Download the data locally or via this command from github
# load("jData.rda")
load(url("https://github.com/siacus/swbCovid/blob/main/rda/jData.rda?raw=true"))


var.jap <- c("swbj", "jDeaths", "jCases", "NIKKEI",
             "jCoronaVirus", "jCoronaVirusNews",
             "jCorona", "jCoronaNews", "jCovid",
             "jCovidNews", "jUnemployment", "jUnemploymentNews",
             "jEconomy", "jEconomyNews", "jGDP", "jGDPNews",
             "jStress", "jDepression", "jHealth", "jSolitude",
             "jInsomnia", "jResidential", "jWorkplace",
             "jpm25", "jtemperature",
             "jWuhan", "jAdultContent", "jlockdown",
             "jFB.CLI", "jFB.ILI", "jFB.MC", "jFB.DC", "jFB.HF"
) 



datj <- data.frame(na.approx(jtmp["2020-01-01/2020-09-30",var.jap],rule=2))
datj <- data.frame( scale(datj) )


jpn_formula ="
VirusSearch =~  jCoronaVirus + jCovid       + jCorona  
PsySearch =~  jStress    + jInsomnia    + jSolitude +jDepression
HealthStatus =~ jFB.CLI + jFB.ILI  
Mobility =~ jResidential + jWorkplace + jlockdown 
Finance =~ NIKKEI + jFB.HF + jUnemployment 
SocDist =~  jFB.MC + jFB.DC 
WellBeing =~ swbj
WellBeing ~ VirusSearch + HealthStatus + Mobility + Finance + SocDist  
Mobility ~~ PsySearch + SocDist  +jCases
PsySearch  ~ WellBeing
jAdultContent ~ WellBeing
"

sem.jpn <- sem(jpn_formula, datj, std.lv = TRUE)
pars.jpn <- summary(sem.jpn, fit.measures=FALSE)



pars.jpn$PE$exo <- NULL 

mymod <- sem.jpn


pdf("sem-jpn-paper.pdf",width = 20,height=15)
semPaths(mymod, style ="lisrel", 
         whatLabels = "est", edge.label.cex = .5,  
         #nodeLabels = mylab,
         nCharNodes = 20,sizeMan2=2,
         nDigits=2, #levels=c(-20,-10,0,10,20),
         label.prop=0.7, edge.label.color = "black", rotation = 2, 
         equalizeManifests = TRUE, optimizeLatRes = TRUE, node.width = 1.5, 
         edge.width = 0.6, shapeMan = "rectangle", shapeLat =  "ellipse", 
         shapeInt = "triangle", sizeMan = 4, sizeInt = 4, sizeLat = 6, 
         curve=1.5, unCol = "#070b8c",layout = "tree",
         residuals=FALSE, exoCov = TRUE, covAtResiduals = FALSE,springLevels=TRUE,
         exoVar=FALSE, intercepts=FALSE)
dev.off()



