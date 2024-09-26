distanceDate<-function(nowday,nowmonth,nowyear, startday=1, startmonth=1, startyear=2013){
  if((is.na(nowday))||(is.na(nowmonth))||(is.na(nowyear)))
  {
    return(NA)
  }
  if(nowday<startday){
    endstartmonth<-monthLength(startmonth)-startday
    startmonth<-startmonth+1
    distday<-endstartmonth+nowday
  }
  else{
    distday<-nowday-startday
  }
  if (nowmonth<startmonth)
  {
    endstartyear<-sum(monthLength(startmonth:12))
    startyear<-startyear+1
    distmonth<-endstartyear+sum(monthLength(1:(nowmonth-1)))
  }
  else{
    if(startmonth==nowmonth){
      distmonth=0
    }
    else{
    distmonth<-sum(monthLength(startmonth:(nowmonth-1)))
  }}
  distyear<-(nowyear-startyear)*365
  
  Dist<-distday+distmonth+distyear
  return(Dist)
}

monthLength<-function(month){
  hash<-c()
  hash[12]<-31
  hash[11]<-30
  hash[10]<-31
  hash[9]<-30
  hash[8]<-31
  hash[7]<-31
  hash[6]<-30
  hash[5]<-31
  hash[4]<-30
  hash[3]<-31
  hash[2]<-28
  hash[1]<-31
monthL<-hash[month]
return(monthL)
}


MLSTclean<-function(MLST){
  if (is.na(MLST) )
  {
    return(paste("-","-",sep="\\"))
  }
  if ( str_detect(MLST,"^[0-9]+$"))
  {return(paste(MLST,"-",sep="\\"))}
  
  if ((str_detect(MLST,"Novel"))||(str_detect(MLST,"NOVEL")))
  {
    if (str_detect(MLST,"cannot"))
    {
      return(paste("-","Novel",sep="\\"))
    } 
    if (str_detect(MLST,"Closest"))
    {
      tmp<-str_split(MLST,"\\:")
      tmptmp<-str_split(tmp[[1]][2]," ")
      return(paste(tmptmp[[1]][1],"Novel",sep="\\"))
    }
  }
  if ( str_detect(MLST,"\\*"))
  {
    pippo<-str_replace_all(MLST,"\\*","")
    return(paste(pippo,"*",sep="\\"))
  }
  if (str_detect(MLST,"Failed"))
  {
    return(paste("-","Failed",sep="\\"))
  }
  if (str_detect(MLST,""))
  {
    return(paste("-","-",sep="\\"))
  }
}
