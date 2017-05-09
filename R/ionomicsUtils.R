#################################################################
#This file contains commonly used R ionomics functions
#Author: gziegler
#################################################################

###########################################
#remove outliers from the provided vector and return vector with outliers changed to NA
#x <- vector of data points to remove outliers from
#mcut <- Number of MADs a data point need to be from the median to be considered an outlier, default is 6.2
# See: Davies, P.L. and Gather, U. (1993).
# "The identification of multiple outliers" (with discussion)
# J. Amer. Statist. Assoc., 88, 782-801.
`is.outlier` <- function (x,mcut=6.2) {
  y <- na.omit(x)
  lims <- median(y) + c(-1, 1) * mcut * mad(y, constant = 1)
  for(j in 1:length(x)){
    if(is.na(x[j]) | x[j] < lims[1] | x[j] > lims[2]){
      x[j] <- NA
    }
  }
  return(x)	
}

#########################################
#calls is.outlier for all of the data columns (cols) of the provided data.frame (x)
# and returns the data frame with NA in place of outliers
#
#x <- dataset
#mcut <- Number of MADs a data point need to be from the median to be considered an outlier, default is 6.2
#by  <- column name to group data by for outlier removal (e.g. by line, by run, etc.), if not provided
#       then by whole dataset
#cols <- vector of column numbers in x to remove outliers from
#Example usage: outlierRemoveDataset(expdata,15,"ICP.run",grep("_intensity",colnames(expdata)))
`outlierRemoveDataset` <- function (x,mcut=6.2,by=NA,cols){
  for(i in cols){
    if(is.na(by)){
      x[,i] <- is.outlier(x[,i],mcut)
    }else{
      for(j in unique(x[,by])){
        if(is.na(j)){
          x[is.na(x[,by]), i] <- is.outlier(x[is.na(x[,by]), i],mcut)
        }else{
          x[x[,by] == j & !(is.na(x[,by])), i] <- is.outlier(x[x[,by] == j & !(is.na(x[,by])), i],mcut)
        }
      }
    }
  }
  return(x)
}

#This will calculate the RSD breaks for data points

#add run breaks
#make axis start and end at end of points
`scatterPlot` <- function(scatter,rsd=NA,shape=NA,main,xlab="Sample No.",ylab="Concentration",runBreaks=NA){
  #suppressPackageStartupMessages(require(ggplot2))
  #suppressPackageStartupMessages(require(grid))
  data <- data.frame(scatter = scatter,rsd = rsd, shape = shape)  
  smooth <- stats::lowess(na.omit(data[,1]),f=0.1)$y
  if(any(is.na(data$scatter))){
    for(na.val in which(is.na(data$scatter))){
      smooth <- append(smooth,NA,after=na.val-1)
    }
  }
  data$smooth <- smooth
  if(!sum(is.na(data$rsd))==length(data$rsd)){
    data$colorCol <- cut(as.numeric(data$rsd),breaks=c(2*(0:5),Inf),labels=c("<2","2-4","4-6","6-8","8-10",">10"),include.lowest=TRUE)
  }
  data$x <- 1:nrow(data)
  p1 <- ggplot(data=data,aes(x=x,y=scatter))  
    if(!sum(is.na(data$rsd))==length(data$rsd) & !sum(is.na(data$shape))==length(data$shape)){ #both shape and rsd
      p1 <- p1 + geom_point(aes(colour=colorCol,shape=shape))
      p1 <- p1 + scale_colour_manual(values=colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(7),name="RSD")
      p1 <- p1 + scale_shape_manual(name = "Type", values = c(17,16))      
    }else if(!sum(is.na(data$rsd))==length(data$rsd)){ #only rsd
      p1 <- p1 + geom_point(aes(colour=colorCol))
      p1 <- p1 + scale_colour_manual(values=colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(7),name="RSD")      
    }else if(!sum(is.na(data$shape))==length(data$shape)){ #only shape
      p1 <- p1 + geom_point(aes(shape=shape))
      p1 <- p1 + scale_shape_manual(name = "Type", values = c(17,16))
    }else{ #neither
      p1 <- p1 + geom_point()
    }
    if(nrow(data) > 50){
      p1 <- p1 + geom_line(aes(y=smooth),size=1.3,colour="blue")
    }
    p1 <- p1 + labs(x = xlab,y = ylab,title = main)
    p1 <- p1 + theme(legend.background=element_rect(),legend.position="bottom",legend.box="horizontal",legend.direction="horizontal")
    p1 <- p1 + scale_x_continuous(limits=c(0,nrow(data)),expand=c(0,0))
    if(!is.na(runBreaks)){
      p1 <- p1 + geom_vline(xintercept = runBreaks)
    }
    p1 <- p1 + geom_hline(yintercept = mean(data[,1],na.rm=TRUE),colour="red",size=1)
    #print(p1)
    return(p1)
    #grid.text(0.2, unit(0.035,"npc"), label=paste("mean=",signif(mean(data[,1],na.rm=TRUE),2),sep=""),gp=gpar(col="red",fontsize=12))
}

# Multiple plot function
#
# I did not write this function, it is from:
# http://wiki.stdout.org/rcookbook/Graphs/Multiple%20graphs%20on%20one%20page%20%28ggplot2%29/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
# usage: multiplot(p1, p2, p3, p4, cols=2)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  #require(grid)
  
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

#pull out max and min sample and boxplot

#sorted boxplot of sample replicates

#xyplot of replicates

`weightEstimation` <- function(x,y,predictx,maxRSD=25,maxEl=10,df=500/9){
  #calculate normalized concentrations for provided data
  normX <- x
  for(i in colnames(x)){ #foreach element
    normX[,i] <- (x[,i]/y)*(df)
  }  
  elRSD <- sapply(normX,sd,na.rm=TRUE)/sapply(normX,mean,na.rm=TRUE)*100
  elMean <- sapply(normX,mean,na.rm=TRUE)
  meetCutoff <- elRSD[which(elRSD < maxRSD)]
  if(length(meetCutoff) < maxEl){
    numEl <- length(meetCutoff)
  }else{
    numEl <- maxEl
  }
  bestElements <- names(sort(meetCutoff)[1:numEl])
  
  weightEstimates <- predictx
  for(i in colnames(predictx)){
    weightEstimates[,i] <- (predictx[,i]/elMean[i])*(df)
  }
  
  #Changed these steps to use medians to be a little more robust to a single bad concnetration measurement
  #so each pass removes values further than 1 mad from the median
#  firstPass.avgWeight <- rowMeans(weightEstimates[bestElements], na.rm = TRUE) #calculate the means
  #firstPass.avgWeight <- apply(weightEstimates[bestElements],1,median,na.rm=T) #calculate the medians
  #remove outliers at 3 MAD    
  weightEstimates <- as.data.frame(t(apply(weightEstimates[bestElements],1,is.outlier,mcut=3)))
#   for(i in colnames(weightEstimates)){
#     weightEstimates[i][weightEstimates[i] > 1.5*firstPass.avgWeight | weightEstimates[i] < .5*firstPass.avgWeight] = NA
#   }
#   
#   secondPass.avgWeight <- rowMeans(weightEstimates[bestElements], na.rm = TRUE)
#   for(i in colnames(weightEstimates)){
#     weightEstimates[i][weightEstimates[i] > 1.5*secondPass.avgWeight | weightEstimates[i] < .5*secondPass.avgWeight] = NA
#   }
  
  finalPrediction <- rowMeans(weightEstimates[bestElements], na.rm = TRUE)
  return(list(weights=as.vector(finalPrediction),els=bestElements))
}

#return the linear model equation of x~y formatted for plotting using lm_eqn(df))
#Example usage: p + annotate("text", label=lm_eqn(df), parse=TRUE, x=Inf, y=Inf, hjust=1.1, vjust=1.5)
lm_eqn = function(df){
  m = lm(df[,1] ~ df[,2])
  paste("italic(y)==",format(coef(m)[1], digits = 2),"+",format(coef(m)[2], digits = 2),"%.%italic(x)*\",\"~~italic(r)^2==",format(summary(m)$r.squared, digits = 3),sep="") 
}


#label outliers from a linear model in a plot (e.g. those with a std. resid. > 2) (note that plot.lm defaults to labelling 3 ids, regardless of extremity)
#df <- data.frame(x = or.data$B11_normConc[1:100],y = or.data$B11_corrConc[1:100],label=or.data$sample[1:100])
`labelOutliers` <- function(df){
  #require(ggplot2)
  df <- na.omit(df)
  m = lm(df[,"x"] ~ df[,"y"])
  df.fortified = fortify(m)
  
  #select extreme values
  df$extreme = ifelse(abs(df.fortified$`.stdresid`) > 2, 1, 0)
  p <- ggplot(df,aes(x=x,y=y))+geom_point()+geom_text(data=df[df$extreme==1,],aes(x=x,y=y,label=label),size=3,hjust=-.3)+annotate("text", label=lm_eqn(df), parse=TRUE, x=Inf, y=Inf, hjust=1.1, vjust=1.5)
  return(p)
}

#Visualize where replicates occur in an ionomics run
#takes a vector of when each replicate occurred (e.g. if two samples A and B were run twice in this order A,B,A,B the vector would be 1,1,2,2
# if the replicates were run A,A,B,B the vector would be 1,2,1,2)
# #example usage
# expdata <- read.table("../BL12-015 - 06PR/outfiles/120717_Data.csv", sep = ",", header = TRUE, na.strings="NA",stringsAsFactors = FALSE)
# numReps <- as.data.frame(table(expdata$sample))
# for(i in numReps$Var1){
#   expdata$repNum[expdata$sample == i & !is.na(expdata$sample)] <- c(1:numReps$Freq[numReps$Var1 == i])
# }  
# visReps(expdata$repNums)
`visReps` <- function(vec){
  #require(ggplot2)
  vec <- na.omit(vec)
  df <- data.frame(sampleNum=factor(1:length(vec),ordered=TRUE),repNum=as.factor(vec),one=1)
  p <- ggplot(df,aes(x=sampleNum,y=one,fill=repNum,colour=repNum,width=.5))+geom_bar(stat="identity",position="dodge")+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                      axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                      #axis.title.x=element_blank(),
                                                                      axis.title.y=element_blank(),
                                                                      panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                                      panel.grid.minor=element_blank(),plot.background=element_blank())  
  print(p)
}
#histogram of NormConcs, colored by binned RSDs
# data <- expdata[expdata$Co59_normConc < 0.015 & expdata$Co59_normConc > 0,c("Co59_normConc","Co59_RSD")]
# colnames(data) <- c("scatter","rsd")
# data$colorCol <- cut(as.numeric(data$rsd),breaks=c(2*(0:5),Inf),labels=c("<2","2-4","4-6","6-8","8-10",">10"),include.lowest=TRUE)
# hist_rsd <- ggplot(data,aes(x=scatter,fill=colorCol))
# hist_rsd + geom_bar(binwidth=(range(data$scatter,na.rm=T)[2]-range(data$scatter,na.rm=T)[1])/60)
# 
# #or look at the distributions broken down by rsd
# hist_rsd + geom_density(alpha=0.2)

###Take list of elemental symbols or symbol isotope pairs and return a list of the full element name###
##Symbols not in the list will be returned unchanged with a warning
`SymboltoElementName` <- function(els){
  els <- sub("_\\w+","",els)
  els <- sub("\\d+","",els)
  elTable <- data.frame(Element = c("Boron","Sodium","Magnesium","Aluminum","Phosphorus","Sulfur",
                                    "Potassium","Calcium","Manganese","Iron","Cobalt","Nickel",
                                    "Copper","Zinc","Arsenic","Selenium","Rubidium","Strontium",
                                    "Molybdenum","Cadmium","Indium","Yttrium","Lead"),
                        Symbol = c("B","Na","Mg","Al","P","S","K","Ca","Mn","Fe","Co","Ni","Cu","Zn","As","Se","Rb","Sr","Mo","Cd","In","Y","Pb"), stringsAsFactors = FALSE)
  out <- elTable[match(els,elTable$Symbol),"Element"]
  if(length(which(is.na(out)))>length(which(is.na(els)))){warning("Some values from input not found in symbol lookup table. Original value returned")}
  out[which(is.na(out))] <- els[which(is.na(out))]
  out
}