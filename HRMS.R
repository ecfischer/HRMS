library(DescTools)

readHRMS<-function(x){
  data<-read.delim2(x, header=F, col.names=c("mass","int")) #This reads the table in. Tab delim.
  data[,1]<-as.numeric(as.character(data[,1])) #These lines changes the strings to numerics... they have to go thru characters first... dunno why
  data[,2]<-as.numeric(as.character(data[,2]))
  return(data)
}

relativeHRMS<-function(data){
  #This functions takes a data table from readHRMStxt and converts the MS int (ie. intensity) to relative numbers in percentage.
  data_rel<-data #these lines converts the data to relative intensity
  data_rel$int<-data_rel$int/max(data_rel$int)*100
  return(data_rel)
}

subsetHRMS<-function(data, from, to){
  to<-as.numeric(to)
  from<-as.numeric(from)
  a<-data[data$mass <= to & data$mass >= from, ]
  return(a)
}

findmaxHRMS<-function(data, from=min(data[,1]), to=max(data[,1])){
  #finds the maximum coordinates in a given range
  sub<-subsetHRMS(data, from, to)
  return(sub[which(sub$int==max(sub$int)),])
}

findpeakHRMS<-function(data){
  spl<-smooth.spline(data$int~data$mass) #first it makes a smoothened line of the data
  peak_x<-findmaxHRMS(data)[,1] #finds the highest values mass 
  pred1<-predict(spl, x=peak_x, deriv=1) # finds the 1st derivate using the smoothened line for x
}


addslopesHRMS<-function(data){
  #This function creates a smoothened approximation for the curve in the given data set ... and gives back the data set with slopes for each data point.
  spl<-smooth.spline(data$int~data$mass, spar=0.01) #first it makes a smoothened line of the data
  slopes=c()
  for (masses in data$mass){
    pred1<-predict(spl, x=masses, deriv=1)
    slopes<-c(slopes,pred1$y)
  }
  result<-cbind(data,slopes)
  return(result)
}


integrateHRMS<-function(data, thres, thres2=0.01, from=min(data[,1]), to=max(data[,1]), plottype="l", output_text=TRUE, output_plot=TRUE){
  #This function takes a data frame and optionally from/to intervals. It will integrate the highest peak in that frame and give an output of the integral as well as coordinates for the cutoff points defining the peak. 
  #You can specify T/F if you want a text output and/or a plot output for inspecting the cutoff points.
  #The threshhold for peak cutoffs is optimized as 5% (change using thres2=0.01) of the highest slope in the data. A threshhold slope can also be defined using <thres>. 
  
  #preparing the input data frame
  data1<-subsetHRMS(data, from, to)
  data1<-addslopesHRMS(data1)
  peak_x<-findmaxHRMS(data1)[,1] #finds the highest values mass 
  upper_sub<-subsetHRMS(data1, from=peak_x, to=max(data1[,1])) #makes a subset of the upper half of the peak.
  lower_sub<-subsetHRMS(data1, to=peak_x, from=min(data1[,1]))
  lower_sub<-lower_sub[order(-lower_sub$mass),]
  
  #Determining appropiate threshold if no thres was given
  if (missing(thres)){
    thres<-max(abs(data1$slopes))*thres2 #threshhold is defined as (default) 5% of the maximum slope in the peak. 
  }
  
  #The upper part of the max peak
  upper_limit_x=0
  for (i in seq(1,length(upper_sub[,1])-1)){
    if (upper_sub[i+1,3] >= -thres & upper_sub[i+1,3] >= upper_sub[i,3]){ #If the slope is higher than -threshhold AND the slope is higher than the slope for the previous data point...It will determine i+1 as the limit. 
      upper_limit_x<-upper_sub[i+1,1]
      break}
  }
  #The lower part of the max peak
  lower_limit_x=0
  for (i in seq(1,length(lower_sub[,1])-1)){
    if (lower_sub[i+1,3] <= thres & lower_sub[i+1,3] <= lower_sub[i,3]){ #If the slope is lower than threshhold AND the slope is lower than the slope for the previous data point...It will determine i+1 as the limit. 
      lower_limit_x<-lower_sub[i+1,1]
      break}
  }
  if (lower_limit_x == 0 | upper_limit_x == 0){
    print("Threshhold was not met. Please expand data set.")
  }
  
  #Outputs
  peak_sub<-subsetHRMS(data, from=lower_limit_x, to=upper_limit_x)
  int<-AUC(x=peak_sub$mass,y=peak_sub$int)
  result<-list("integrale"=int, 
               "mass"=peak_x,
               "lower"=data1[which(data1$mass==lower_limit_x),], 
               "upper"=data1[which(data1$mass==upper_limit_x),])
  if (output_text==TRUE){ #text output
    print(paste0("The determined peak interval is: ", toString(lower_limit_x),"-",toString(upper_limit_x)," Da"))
    print(paste0("The integrale of the peak is: ", toString(int), " Da"))
    print(paste0("Threshhold +/-: ", toString(thres)))
  }
  if (output_plot==TRUE){ #visual output
    spl<-smooth.spline(data1$int~data1$mass, spar=0.01)
    plot(x=data1$mass,y=data1$int, type=plottype)
    lines(spl)
    abline(v=upper_limit_x, col="Red", lty=2)
    abline(v=lower_limit_x, col="Red",lty=2)
    arrows(angle=90, code=3, length=0,
           x0=result$lower[,1],
           y0=result$lower[,2],
           x1=result$upper[,1],
           y1=result$upper[,2])
    mtext(paste("Peak",toString(peak_x)," Da", "    Integrale", toString(int)), side=3)
  }
  return(result)
}

runpeakHRMS<-function(filenames, peak, thres2=0.05, size=10, output_text=F, output_plot=T){
  #This functions runs thru multiple files and collects the integral of a given peak. The threshold for the peak can be specified by thres2 (as in integrateHRMS). 
  #The function ask you to press any key every time it integrates a peak to verify the integration ... but only if any of the output types are set to TRUE.
  #It will return a list of filenames with integrals. 
  peakto<-as.numeric(peak+size)
  peakfrom<-as.numeric(peak-size)
  result=list()
  for (files in filenames){
    dat<-readHRMS(files)
    int<-integrateHRMS(dat, to=peakto, from=peakfrom, output_text=output_text, thres2=thres2, output_plot=output_plot)
    result<-append(result,list(int$integrale))
    if (sum(output_text,output_plot)>0){
      readline(prompt="Press [enter] to continue")
    }
  }
  names(result)<-filenames
  return(result)
}


findmaxintHRMS<-function(filenames){
  a<-c()
  for (file in filenames){
    new<-findmaxHRMS(readHRMS(file))
    a<-c(a,new$int)}
  return(max(a))
}

findmassrangeHRMS<-function(filenames){
  top<-c()
  bottom<-c()
  for (file in filenames){
    data<-readHRMS(file)
    data$mass
    top<-c(top, max(data$mass))
    bottom<-c(bottom, min(data$mass))
  }
  return(list(
    "min"=min(bottom),
    "max"=max(top)
  ))
}

offsetHRMS<-function(data, reference=findmaxHRMS(data)$mass){
  data2<-data
  data2$mass<-data$mass-reference
  return(data2)
}

detectpeaksHRMS<-function(file, from=lower, to=upper, rel=T, grains=5, minpeak=7, int_size=10, thres2=0.15){
  ### This function returns masses, intensities, and integrals for peaks above a value (minpeak which is a percentage).
  ### Is first loads and simplifies a plot (with grains meaning every <grains> oberservation only). 
  spec<-readHRMS(file)
  lower<-min(spec$mass)
  upper<-max(spec$mass)
  spec<-subsetHRMS(spec, to=upper, from=lower)
  if (rel==TRUE){spec<-relativeHRMS(spec)}
  else {minpeak<-max(spec$int)*0.01*minpeak}
  
  ### Now it simplifies the plot and cuts off everything below minpeak. 
  grains<-5
  a<-seq(from=1,to=length(spec$mass), grains)
  spec2<-spec[a,]  
  spec2<-subset(spec2, spec2$int > minpeak)
  
  ### Then is collapses the adjacent (mass) to only one potential peak per area. 
  peaks<-c()
  local<-c()
  diffs <- diff(spec2[,1])
  for (i in seq(1, length(spec2[,1])-1)){
    local<-c(local, spec2[i,1])
    if (diffs[i] > 5 | i == length(spec2[,1])-1){
      maxpeak<-round(median(local))
      peaks<-c(peaks,maxpeak)
      local<-c()}}
  ### To get the actual peaks ... and not the potential peaks from the simplified plot... we use integrateHRMS. 
  masses<-c()
  ints<-c()
  integrals<-c()
  for (element in peaks){
    info<-integrateHRMS(spec, from=element-int_size, to=element+int_size,thres2=thres2, output_text=F, output_plot=F)
    masses<-c(masses, round(info$mass, digits=1))
    integrals<-c(integrals, round(info$integrale, digits=1))
    index_int<-which(spec$mass==info$mass)
    ints<-c(ints, round(spec[index_int,2], digits=1))}
  results<-list("mass"=masses,"intensity"=ints,"integrale"=integrals)
  return(results)}

plotHRMS<-function(filenames, from, to, rel=TRUE, int.from=0, int.to=upperlimit, ticks.mass=100,
                    style.line=c("solid","solid","solid", "solid","dotted","dotted","dotted","dotted"),
                    style.col=c("black", "blue","red","orange","black", "blue","red","orange"),
                    labels.add=T, labels.ref=1, labels.minpeak=7, labels.thres2=0.15){
  #This functions takes one or more file names and plots the spectra together. By default rel is TRUE which normalizes all spectra to their maximum intensity (rel int). 
  #If rel is false it will take the max intensity value in all of the spectra and use that as an upper limit. 
  #The y-limits can be manipulated manually set with int.to and int.from.
  #Set the tick interval for mass by ticks.mass
  upperlimit<-105
  if (rel == F){upperlimit<-findmaxintHRMS(filenames)*1.05}
  howmany<-length(filenames)
  intrange<-findmassrangeHRMS(filenames)
  if (missing(from)){from<-intrange$min}
  if (missing(to)){to<-intrange$max}
  plot(x=c(from,to), y=c(int.from,int.to), type="l", col="white", axes=F, xlab="mass (Da)", ylab="Intensity")
  axis(side=1, at=seq(from,to,by=ticks.mass), pos=-0.02*int.to, las=2)
  axis(side=2, las=2)
  legend("topright", col=style.col[1:howmany], lty=style.line[1:howmany], bty="n", legend=filenames)
  for (i in 1:howmany){
    if (rel==TRUE){dataset<-subsetHRMS(relativeHRMS(readHRMS(filenames[i])),from=from, to=to)}
    else {dataset<-subsetHRMS(readHRMS(filenames[i]),from=from, to=to)}
    points(dataset, col=style.col[i], type="l", lty=style.line[i])
  }
  if (labels.add == TRUE){
    file<-filenames[labels.ref]
    detect<-detectpeaksHRMS(file, from=from, to=to, rel=rel, minpeak=labels.minpeak, thres2=labels.thres2)
    masses<-detect$mass
    ints<-detect$intensity+(0.025*upperlimit)
    text(x=masses, y=ints, labels=masses, cex=0.8)
  }
}

majorminorHRMS<-function(majorlist, minorlist){
  #This function computes the relative integrants of a list of peaks from runpeakHRMS... 
  #It computes the proportion of minorlist peaks. 
  a<-c()
  names<-c()
  for (file in names(majorlist)){
    names<-c(names,file)
    majorpeak<-as.numeric(majorlist[file])
    minorpeak<-as.numeric(minorlist[file])
    a<-c(a,minorpeak/(minorpeak+majorpeak))
  }
  names(a)<-names
  return(a)
}

plotoffsetHRMS<-function(filenames, from, to, rel=TRUE, int.from=0, int.to=upperlimit, ticks.mass=100,
                         style.line=c("solid","solid","solid", "solid","dotted","dotted","dotted","dotted"),
                         style.col=c("black", "blue","red","orange","black", "blue","red","orange"),
                         reference){
  #This functions takes one or more file names and plots the spectra together. By default rel is TRUE which normalizes all spectra to their maximum intensity (rel int). 
  #If rel is false it will take the max intensity value in all of the spectra and use that as an upper limit. 
  #The y-limits can be manipulated manually set with int.to and int.from.
  #Set the tick interval for mass by ticks.mass
  upperlimit<-100
  if (rel == F){upperlimit<-findmaxintHRMS(filenames)}
  howmany<-length(filenames)
  
  massrange<-findmassrangeHRMS(filenames)
  if (missing(from)){from<-massrange$min}
  if (missing(to)){to<-massrange$max}
  
  if (missing(reference)){reference<-findmaxHRMS(readHRMS(filenames[1]))$mass}
  offset_from<-round_any(from-reference, ticks.mass, f=ceiling)
  offset_to<-round_any(to-reference, ticks.mass, f=floor)
  xvalues<-c(seq(offset_from, 0, by=ticks.mass), seq(ticks.mass,offset_to, by=ticks.mass))
  
  plot(x=c(offset_from,offset_to), y=c(int.from,int.to), type="l", col="white", axes=F, xlab="mass (Da)", ylab="Intensity")
  axis(side=3, at=seq(from,to,by=ticks.mass), pos=0.02*int.to, las=2)
  axis(side=2, las=2)
  axis(side=1, at=xvalues, pos=-0.02*int.to, las=2)
  legend("topright", col=style.col[1:howmany], lty=style.line[1:howmany], bty="n", legend=filenames)
  
  for (i in 1:howmany){
    if (rel==TRUE){dataset<-subsetHRMS(relativeHRMS(readHRMS(filenames[i])),from=from, to=to)}
    else {dataset<-subsetHRMS(readHRMS(filenames[i]),from=from, to=to)}
    dataset<-offsetHRMS(dataset, reference=reference)
    points(dataset, col=style.col[i], type="l", lty=style.line[i])
  }
}

#Define notin function for good measure
`%notin%` <- Negate(`%in%`)





