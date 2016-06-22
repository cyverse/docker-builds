pie_compare<-function(gene.name,n1,estimates.file="estimates.txt",geometry="Euclidean",adjust.weight=1E-300,output.file=paste("Pieplot_",gene.name,".pdf",sep=""),group.name=c("1","2")){
  estimates<-read.delim(estimates.file,stringsAsFactors=FALSE,comment.char="#")
  estimates.gene<-estimates[estimates[,1]==gene.name,,drop=FALSE]
  
  if (nrow(estimates.gene)==0 | sum(is.na(estimates.gene[1,-(1:2)]))==(ncol(estimates.gene)-2)){
    try("No data for the input gene! \n")
  }
  else{
    isoforms<-estimates.gene[,2]
    
    usage.data<-t(estimates.gene[,-(1:2)])
    usage.1<-usage.data[1:n1,,drop=FALSE]; usage.1<-usage.1[complete.cases(usage.1),,drop=FALSE]
    usage.2<-usage.data[-(1:n1),,drop=FALSE]; usage.2<-usage.2[complete.cases(usage.2),,drop=FALSE]
    
    if (nrow(usage.1)==0){
      if (geometry=="Euclidean"){
        temp<-matrix(colMeans(usage.2),nrow=1)
      }
      else {
        usage.2.m<-apply(usage.2,1:2,function(j){ifelse(j==0,adjust.weight,j)})
        temp<-matrix(ilrInv(matrix(colMeans(ilr(usage.2.m)),nrow=1)),nrow=1)
      }
      pie_plot(temp,isoforms,gene.name,output.file,group=group.name[2])
    }
    else if (nrow(usage.2)==0){
      if (geometry=="Euclidean"){
        temp<-matrix(colMeans(usage.1),nrow=1)
      }
      else {
        usage.1.m<-apply(usage.1,1:2,function(j){ifelse(j==0,adjust.weight,j)})
        temp<-matrix(ilrInv(matrix(colMeans(ilr(usage.1.m)),nrow=1)),nrow=1)
      }
      pie_plot(temp,isoforms,gene.name,output.file,group=group.name[1])
    }
    else{
      if (geometry=="Euclidean"){
        temp<-rbind(matrix(colMeans(usage.1),nrow=1),matrix(colMeans(usage.2),nrow=1))
      }
      else{
        usage.1.m<-apply(usage.1,1:2,function(j){ifelse(j==0,adjust.weight,j)}); usage.2.m<-apply(usage.2,1:2,function(j){ifelse(j==0,adjust.weight,j)});
        temp<-rbind(matrix(ilrInv(matrix(colMeans(ilr(usage.1.m)),nrow=1)),nrow=1),matrix(ilrInv(matrix(colMeans(ilr(usage.2.m)),nrow=1)),nrow=1))
      }
      pie_plot(temp,isoforms,gene.name,output.file,group=group.name)
    }
  }
}

