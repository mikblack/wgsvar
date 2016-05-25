rm(list=ls())
library(DNAcopy)
dt<-format(Sys.time(), "%Y%m%d")
cf<-paste("calledFiles-",dt,".txt",sep='')
sysCmd<-paste("ls -1 *called > ",cf,sep='')
system(sysCmd)
ff<-readLines(cf)
sn<-ff
pltdir<-paste('CNVplots-',dt,sep='')
rdatdir<-paste('CNVrdata-',dt,sep='')
system(paste('mkdir -p',pltdir))
system(paste('mkdir -p',rdatdir))

cn<-c()
for(i in 1:length(sn)){
    aa<-read.table(ff[i],header=T)
    aa$adjusted_log_ratio<- aa$adjusted_log_ratio - median(aa$adjusted_log_ratio)
    cn <- rbind(cn,aa)
}
CNA.object <-CNA( genomdat = cn[,"adjusted_log_ratio"], 
                 chrom = gsub("chr","",as.vector(cn[,"chrom"])), maploc = cn[,"chr_start"], 
                 data.type = 'logratio', sampleid = sn[i])
CNA.smoothed <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(CNA.smoothed)

segs2 = segment.smoothed.CNA.object$output
gg<-paste("daisy5-allchrom-CNV-CBS-",dt,".txt",sep='')
ll<-segs2$loc.end-segs2$loc.start
segs2$length<-ll
write.table(segs2[,c(2:4,7,6)], file=gg, row.names=F, col.names=T, quote=F, sep="\t")
cnvr<-seq(1,nrow(segs2))[segs2$length>1000 & 
                           (segs2$seg.mean>log2(3/2) | segs2$seg.mean< -1)]
gg<-paste("daisy5-allchrom-CNVR-",dt,".txt",sep='')
write.table(segs2[cnvr,c(2:4,7,6)], file=gg, row.names=F, col.names=T, quote=F, sep="\t")

pdf(paste(pltdir,"/daisy-cnvplot-allchrom-",dt,".pdf",sep=''),height=6,width=12)
plot(segment.smoothed.CNA.object, plot.type="w",main='')
dev.off()

pdf(paste(pltdir,"/daisy-cnvplot-bychrom-",dt,".pdf",sep=''),height=18,width=18)
plot(segment.smoothed.CNA.object, plot.type="s",main='')
dev.off()

save.image(paste(rdatdir,'/D5SP_D6TL-HiSeq_dnaCopy.RData',sep=''))
