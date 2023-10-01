############################## -*- Mode: Ess-R -*- ############################
## 
## Filename: test.r
## 
## Title     :
## 
## Created   : <2022-09-23 (Fri) 15:24:46 by Hyo-Jun Kim>
## Time-stamp: <2022-09-23 (Fri) 16:27:53 by Hyo-Jun Kim>
## 
## URL: 
## Keywords: 
## 
######################################################################
## 
### Commentary: 1. Try to derive a good comparison target
##              2. Use a Classification Model to tell the similarity between DNN and DEA
## 
## 
######################################################################
## 
### Change Log:
## 
## 
######################################################################
## 
### Code:

library(nonparaeff)
library(dplyr)
library(readr)
library(tibble)
library(rlang)
library(geometry)


half.rnorm <- function(...){
  re <- abs(rnorm(...))
  return(re)
}

#### DGP function for traditional eff  ####

dgp <- function(ndmu = NA, noutput = NA, ninput = NA, orientation = 1,
                lower.bound = 0, upper.bound = 1, theta.dist = NA, ...)
{
  if( is.na(ndmu) | is.na(noutput) | is.na(ninput) | orientation > 2){
    stop("ndmu or noutput or ninput should not be NA, and orientation should be 1 or 2.")
  }
  
  ys <- runif(ndmu*noutput, min = lower.bound, max = upper.bound)
  xs <- runif(ndmu*ninput, min = lower.bound, max = upper.bound)
  
  true.eff <- exp(-theta.dist(ndmu,...))
  ori <- rep(orientation, ndmu)
  
  y.names <- paste("y", 1:noutput, sep= ".")
  x.names <- paste("x", 1:ninput, sep= ".")
  
  if(orientation == 1){
    aug.xs <- xs/true.eff
    my.dat <- as.data.frame(matrix(c(ys, xs, aug.xs, true.eff, ori), nrow = ndmu))
    aug.x.names <- paste("aug.x", 1:ninput, sep = ".")
    names(my.dat) <- c(y.names, x.names, aug.x.names, "true.eff", "ori")
    
  } else {
    aug.ys <- ys*true.eff
    my.dat <- as.data.frame(matrix(c(ys, xs, aug.ys, true.eff, ori), nrow = ndmu))
    aug.y.names <- paste("aug.y", 1:noutput, sep = ".")
    names(my.dat) <- c(y.names, x.names, aug.y.names, "true.eff", "ori")
  }
  return(my.dat)
}

## example
dgp(ndmu = 64, noutput = 5, ninput = 5, orientation = 1, theta.dist =  half.rnorm)


#########################################################################################################################
#########################################################################################################################


#### Use DNN with DGP ####
## In R we provide only DGP model you could check DNN model in python

## DGP Function for DNN
dgp <- function(ndmu = NA, noutput = NA, ninput = NA, orientation = 1,
                lower.bound = 0, upper.bound = 1, theta.dist = NA, ...)
{
  if( is.na(ndmu) | is.na(noutput) | is.na(ninput) | orientation > 2){
    stop("ndmu or noutput or ninput should not be NA, and orientation should be 1 or 2.")
  }
  
  ys <- runif(ndmu*noutput, min = lower.bound, max = upper.bound)
  xs <- runif(ndmu*ninput, min = lower.bound, max = upper.bound)
  
  true.eff <- exp(-theta.dist(ndmu,...))
  #ori <- rep(orientation, ndmu)
  
  y.names <- paste("y", 1:noutput, sep= ".")
  x.names <- paste("x", 1:ninput, sep= ".")
  
  if(orientation == 1){
    aug.xs <- xs/true.eff
    my.dat <- as.data.frame(matrix(c(ys, aug.xs, true.eff), nrow = ndmu))
    #aug.x.names <- paste("aug.x", 1:ninput, sep = ".")
    names(my.dat) <- c(y.names, x.names, "true.eff")
    
  } else {
    aug.ys <- ys*true.eff
    my.dat <- as.data.frame(matrix(c(aug.ys, xs, true.eff), nrow = ndmu))
    #aug.y.names <- paste("aug.y", 1:noutput, sep = ".")
    names(my.dat) <- c(y.names, x.names, "true.eff")
  }
  return(my.dat)
}
## i: ndmu
## j: noutput
## k: ninput
## l: orientation
## m: theta.dist
## n: sample.seed
getwd()

#for(i in c(100,500,1000,5000, 10000, 20000, 40000, 80000,160000,320000,640000,1280000)){
for(i in c(100)){  
  for(j in c(1,2,3,4,5)){
    for(k in c(1,2,3,4,5)){
      for(l in c(1,2)){
        
        sample.dat <- data.frame(matrix(ncol = j+k+1))
        names(sample.dat) <- c( paste("y", 1:j, sep= "."),  paste("x", 1:k, sep= "."), "true.eff")
        
        for(m in c(half.rnorm)){
          for(n in 1){
            
            set.seed(n)
            sample.dat <- rbind(sample.dat,dgp(ndmu = i, noutput = j,
                                               ninput = k,
                                               orientation = l,
                                               theta.dist = m))
            
          }
        }
        
        sample.dat <- sample.dat[-1,]
        
        #sample.dat <- add_column(sample.dat, halfnorm = rep(c(1,0), 25600))
        #sample.dat <- add_column(sample.dat, exp = rep(c(0,1), 25600))
        
        readr :: write_csv(sample.dat, path 
                           = paste0('C:/Users/ocean/Desktop/0902_Re/Data/dataset_before_categorized/',i,'_',j,k,'_',l,'_dat.csv'))
      }
    }
  }
}


#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

#### extract high efficiency DMU from validation Set to use at DEA####
## we tried to extract from the x ## we will add outlier instead
#########################################################################################################################


#for(i in c(100,500,1000,5000, 10000, 20000, 40000, 80000,160000,320000,640000,1280000)){
for(i in c(100)){
  for(j in c(1,2,3,4,5)){
    for(k in c(1,2,3,4,5)){
      for(l in c(1,2)){
        sample.dat <- read.csv(file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/dataset_before_categorized/',i,'_',j,k,'_',l,'_dat.csv'))
        extract.dat <- sample.dat[int(nrow(sample.dat)*0.8+1):int(nrow(sample.dat)),]
        
        if(l == 1){
          if(k == 1){
            extract.dat$average.x <- extract.dat[,int(j+1)]
            
          } else {
            extract.dat$average.x <- rowMeans(extract.dat[,int(j+1):int(j+k)])
          }
          attach(extract.dat)
          topbottom_5p <- quantile(average.x, probs = c(0.05, 0.95))
          detach(extract.dat)
          quantile5_x.dat <- extract.dat[extract.dat$average.x >= topbottom_5p[1] & extract.dat$average.x < topbottom_5p[2] ,]
          DNN.quantile5_x.dat <-quantile5_x.dat[,1:int(j+k+1)]
          write.csv(DNN.quantile5_x.dat,file=paste0('C:/Users/ocean/Desktop/0902_Re/Data/test_dataset/test_df/',i,'_',j,k,'_',l,'_test_dat.csv'))
          
          
        } else {
          if(j == 1){
            extract.dat$average.y <- extract.dat[,1]
          } else {
            extract.dat$average.y <- rowMeans(extract.dat[,1:int(j)])
          }
          attach(extract.dat)
          topbottom_5p <- quantile(average.y, probs = c(0.05, 0.95))
          detach(extract.dat)
          quantile5_y.dat <- extract.dat[extract.dat$average.y >= topbottom_5p[1] & extract.dat$average.y < topbottom_5p[2] ,]
          DNN.quantile5_y.dat <-quantile5_y.dat[,1:int(j+k+1)]
          write.csv(DNN.quantile5_y.dat,file=paste0('C:/Users/ocean/Desktop/0902_Re/Data/test_dataset/test_df/',i,'_',j,k,'_',l,'_test_dat.csv'))
          
        }
        
      }
    }
  }
}



#########################################################################################################################
#########################################################################################################################
#########################################################################################################################


#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

#### extract high efficiency DMU from validation Set to use at DEA####
## we tried to extract from the x 
#########################################################################################################################
## extract data which average of x is top 10%
## estimate DEA but DEA Model doesn't know the Real Efficiency

sample.dat <- read.csv(file = paste0('C:/Users/ocean/Desktop/0902_Re/Data_outlier/dataset_before_categorized/100_11_1_dat.csv'))
extract.dat <- sample.dat[int(nrow(sample.dat)*0.8+1):int(nrow(sample.dat)),]
sample.dat
extract.dat

sample.id <- sample(as.numeric((rownames(extract.dat))[1]):as.numeric((rownames(extract.dat))[nrow(extract.dat)]), int(nrow(extract.dat)/5), replace = FALSE)
typeof(sample.id)
sample.id
new.id <- as.vector(sample.id)
new.id
typeof(new.id)
dat(extract.dat[new.id,])



#for(i in c(160000,320000,640000,1280000)){
for(i in c(100,500,1000,5000)){
  for(j in c(1,2,3,4,5)){
    for(k in c(1,2,3,4,5)){
      for(l in c(1,2)){
        sample.dat <- read.csv(file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/dataset_before_categorized/',i,'_',j,k,'_',l,'_dat.csv'))
        extract.dat <- sample.dat[int(nrow(sample.dat)*0.8+1):int(nrow(sample.dat)),]
        
        if(l == 1){
          if(k == 1){
            extract.dat$average.x <- extract.dat[,int(j+1)]
            
          } else {
            extract.dat$average.x <- rowMeans(extract.dat[,int(j+1):int(j+k)])
          }
          attach(extract.dat)
          topbottom_5p <- quantile(average.x, probs = c(0.05, 0.95))
          detach(extract.dat)
          quantile5_x.dat <- extract.dat[extract.dat$average.x >= topbottom_5p[1] & extract.dat$average.x < topbottom_5p[2] ,]
          DNN.quantile5_x.dat <-quantile5_x.dat[,1:int(j+k+1)]
          write.csv(DNN.quantile5_x.dat,file=paste0('C:/Users/ocean/Desktop/0902_Re/Data/test_dataset/test_df/',i,'_',j,k,'_',l,'_test_dat.csv'))
          
          
        } else {
          if(j == 1){
            extract.dat$average.y <- extract.dat[,1]
          } else {
            extract.dat$average.y <- rowMeans(extract.dat[,1:int(j)])
          }
          attach(extract.dat)
          topbottom_5p <- quantile(average.y, probs = c(0.05, 0.95))
          detach(extract.dat)
          quantile5_y.dat <- extract.dat[extract.dat$average.y >= topbottom_5p[1] & extract.dat$average.y < topbottom_5p[2] ,]
          DNN.quantile5_y.dat <-quantile5_y.dat[,1:int(j+k+1)]
          write.csv(DNN.quantile5_y.dat,file=paste0('C:/Users/ocean/Desktop/0902_Re/Data/test_dataset/test_df/',i,'_',j,k,'_',l,'_test_dat.csv'))
          
        }
        
      }
    }
  }
}


#########################################################################################################################
#########################################################################################################################
## estimate DEA which extracted data( average of x is more than bottom 5% and less than top 5%)
#for(i in c(160000,320000,640000,1280000)){
#  for(j in c(1,2,3,4,5)){
#    for(k in c(1,2,3,4,5)){
#      for(l in c(1,2)){
#        if( j + k < 4){
#          next
#        }
# 32만 DMU 에 대해서 trad_34_2 까지 완료

#for(i in c(10000, 20000, 40000, 80000, 160000, 320000)){
#  for(j in c(1,2,3,4,5)){
#    for(k in c(1,2,3,4,5)){
#      for(l in c(1,2)){
#        if( j + k < 4){
#          next
#        }

## 초기 모델: first model
for(i in c(100,500,1000,5000)){
  for(j in c(1,2,3,4,5)){
    for(k in c(1,2,3,4,5)){
      for(l in c(1,2)){
        if( j + k < 4){
          next
        }
#        if( k == 4 & l == 1){
#          next
#        }
        
        sample.dat <- read.csv(file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/test_dataset/test_df/',i,'_',j,k,'_',l,'_test_dat.csv'))
        sample.dat <- sample.dat[,2:int(ncol(sample.dat)-1)]
        
        
        if (l == 1){
          re.dea.crs.quant <- dea(base = sample.dat, frontier = sample.dat, noutput = j, rts = 1, orientation = l,
                                  onlytheta = T)
          write.csv(re.dea.crs.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/',i,'_',j,k,'_',l,'_trad_crs.csv'))
          
          ##################
          re.dea.vrs.quant <- dea(base = sample.dat, frontier = sample.dat, noutput = j, rts = 2, orientation = l,
                                  onlytheta = T)
          write.csv(re.dea.vrs.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/',i,'_',j,k,'_',l,'_trad_vrs.csv'))
          
          ##################
          
        } else {
          re.dea.crs.quant <- 1/dea(base = sample.dat, frontier = sample.dat, noutput = j, rts = 1, orientation = 2,
                                    onlytheta = T)
          write.csv(re.dea.crs.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/',i,'_',j,k,'_',l,'_trad_crs.csv'))
          
          ##################
          re.dea.vrs.quant <- 1/dea(base = sample.dat, frontier = sample.dat, noutput = j, rts = 2, orientation = 2,
                                    onlytheta = T)
          write.csv(re.dea.vrs.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/',i,'_',j,k,'_',l,'_trad_vrs.csv'))
          
          
        }
        
        eff.unknown.quant <- data.frame(ext_crs_unknown = re.dea.crs.quant$eff,
                                        ext_vrs_unknown = re.dea.vrs.quant$eff)
        write.csv(eff.unknown.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/',i,'_',j,k,'_',l,'_trad_eff.csv'))
        
        
      }
    }
  }
}


# 128_432까지 됨[441-452]

## 최적화 모델(미완성: 지속적으로 메모리 문제 발생) ##

for(i in c(1280000)){
  for(j in c(5)){
    for(k in c(1,2,3,4,5)){
      for(l in c(1,2)){
        
        if( j + k < 4){
          next
        }
        
        sample.dat <- read.csv(file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/test_dataset/test_df/',i,'_',j,k,'_',l,'_test_dat.csv'))
        sample.dat <- sample.dat[,2:int(ncol(sample.dat)-1)]
        
        my.front.idx <- sort(unique(as.vector(convhulln(sample.dat))))
        my.convhull.dat <- sample.dat[my.front.idx, ]
        
        
        if (l == 1){
          
          tmp <- dea(my.convhull.dat, noutput = j, rts = 2, onlytheta = T)
          real.front <- my.convhull.dat[tmp$eff > 1 - 1e-4, ]
          e = 1
          sample.split <- sample.dat[100*(e-1)+1:100*e,]
          re.dea.crs.quant <- dea(base = sample.split, frontier = real.front, noutput = j, rts = 1, orientation = l,
                                  onlytheta = T)
          for(e in 2:int(i*0.18/100)){
            re.dea.crs.quant <- rbind(re.dea.crs.quant, 
                                      dea(base = sample.split, frontier = real.front, noutput = j, 
                                          rts = 1, orientation = l,onlytheta = T))
          }
          write.csv(re.dea.crs.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/',i,'_',j,k,'_',l,'_trad_crs.csv'))
          
          ##################
          
          tmp <- dea(my.convhull.dat, noutput = j, rts = 2, onlytheta = T)
          real.front <- my.convhull.dat[tmp$eff > 1 - 1e-4, ]
          e = 1
          re.dea.vrs.quant <- dea(base = sample.split, frontier = real.front, noutput = j, rts = 2, orientation = l,
                                  onlytheta = T)
          for(e in 2:int(i*0.18/100)){
            re.dea.vrs.quant <- rbind(re.dea.vrs.quant,
                                      dea(base = sample.split, frontier = real.front, noutput = j, rts = 2, orientation = l,
                                          onlytheta = T))
          }
          write.csv(re.dea.vrs.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/',i,'_',j,k,'_',l,'_trad_vrs.csv'))
          
          ##################
          
        } else {
          tmp <- dea(my.convhull.dat, noutput = j, rts = 2, orientation = 2, onlytheta = T)
          real.front <- my.convhull.dat[tmp$eff < 1 + 1e-3, ]
          e = 1
          sample.split <- sample.dat[100*(e-1)+1:100*e,]
          re.dea.crs.quant <- 1/dea(base = sample.split, frontier = real.front, noutput = j, rts = 1, orientation = 2,
                                    onlytheta = T)
          for(e in 2:int(i*0.18/100)){
            re.dea.crs.quant <- rbind(re.dea.crs.quant,
                                      1/dea(base = sample.split, frontier = real.front, noutput = j, rts = 1, orientation = 2,
                                            onlytheta = T))
          }
          write.csv(re.dea.crs.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/',i,'_',j,k,'_',l,'_trad_crs.csv'))
          
          
          tmp <- dea(my.convhull.dat, noutput = j, rts = 2, orientation = 2, onlytheta = T)
          real.front <- my.convhull.dat[tmp$eff < 1 + 1e-3, ]
          e = 1
          re.dea.vrs.quant <- 1/dea(base = sample.split, frontier = real.front, noutput = j, rts = 2, orientation = 2,
                                    onlytheta = T)
          for(e in 2:int(i*0.18/100)){
            re.dea.vrs.quant <- rbind(re.dea.vrs.quant,
                                      1/dea(base = sample.split, frontier = real.front, noutput = j, rts = 2, orientation = 2,
                                            onlytheta = T))
          }
          write.csv(re.dea.vrs.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/',i,'_',j,k,'_',l,'_trad_vrs.csv'))
          
          
        }
        
        eff.unknown.quant <- data.frame(ext_crs_unknown = re.dea.crs.quant$eff,
                                        ext_vrs_unknown = re.dea.vrs.quant$eff)
        
        write.csv(eff.unknown.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/',i,'_',j,k,'_',l,'_trad_eff.csv'))
        
        
      }
    }
  }
}


e = 1
re.dea.vrs.quant <- dea(base = sample.split, frontier = real.front, noutput = j, rts = 2, orientation = l,
                        onlytheta = T)
for(e in 2:int(i*0.18/100)){
  re.dea.vrs.quant <- rbind(re.dea.vrs.quant,
                            dea(base = sample.split, frontier = real.front, noutput = j, rts = 2, orientation = l,
                                onlytheta = T))
}
write.csv(re.dea.vrs.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/640000_44_1_trad_vrs.csv'))
eff.unknown.quant <- data.frame(ext_crs_unknown = re.dea.crs.quant$eff,
                                ext_vrs_unknown = re.dea.vrs.quant$eff)

write.csv(eff.unknown.quant, file = paste0('C:/Users/ocean/Desktop/0902_Re/Data/trad_eff/640000_44_1_trad_eff.csv'))
