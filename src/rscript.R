

do_plot<-function(data,id){
	ti=paste(data[id,1], paste(", score: ",data[id,2],sep=""),sep="")
        plot((as.numeric(data[id+1,3:length(data[1,])-1])),type='n',main=ti,ylab="coverage",xlab="position",col='green',xlim=c(1,400))
        lines((as.matrix(data[id,3:length(data[1,])-1])),col='green')
	lines((as.matrix(data[id+1,3:length(data[1,])-1])),col='red')
	abline(v=200)
}


setEPS()
tmp=as.matrix(read.table('plots.txt',sep='\t',header=F))

i=0
for(id in seq(from=1,to=length(tmp[,1]), by=2)){
  name=paste("VShape",i,sep="")
  name=paste(name,".eps",sep="")
  postscript(name,paper="a4")
  i=i+1
  print(id)
  do_plot(tmp,id)
  dev.off()
}
