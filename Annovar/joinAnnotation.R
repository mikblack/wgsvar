## Script to combine annotation data from multiple chromosomes into a single file (one each for SNV and InDel).

## NB - This script is specific to the analysis performed for the Sleeping Beauty mouse model for leukemia.
##	Use in other settings will require modification of the code below to match the characteristics of the
##      data being analysed.  

## Create variable to hold date information for output filenames
dt <- format(Sys.time(), "%Y%m%d")

##############################################################
## SNV data
##############################################################

## Define chromosome data to be read in.
## Y is excluded here because no somatic exonic variants were found
samples<-paste("chr",c(1:19,"X"),".dedup",sep='')

out<-list()
for(i in seq(samples)){

## Read in SNP (i.e., SNV) data for this chromosome
    aa<-read.table(paste('../Varscan/',samples[i],'.snp',sep=''),header=T)

## Select only SNVs called as somatic
    bb<-aa[aa$somatic_status=="Somatic",]

## Read in annotation data for somatic SNVs
    cc<-read.table(paste(samples[i],'.snp.SOMATIC.anno.variant_function',sep=''))
    dd<-read.table(paste(samples[i],'.snp.SOMATIC.anno.exonic_variant_function',sep=''),sep='\t')
    ee<-matrix("",nrow(cc),2)
    ln<-as.numeric(gsub("line","",dd[,1]))
    ee[ln,1]<-as.vector(dd[,2])
    ee[ln,2]<-as.vector(dd[,3])
    out[[i]]<-cbind(rep(samples[i],nrow(cc)),cc[,1:2],ee,cc[,c(3:8,10)],bb[,c(5:15)],p.adjust(bb[,15],method='holm'),p.adjust(bb[,15],method='BH'),bb[,16:23])
    colnames(out[[i]])[1:12]<-c("sample","type","gene/region","func","info","chr","start","end","ref","var","zygosity","readdepth")
    colnames(out[[i]])[24:25]<-c("fwer_adj_p","fdr_adj_p")
## Output progress per chromosome
## - Total SNVs
## - Non-synomynous SNVs
## - Significant non-synomynous SNVs (FDR adjusted: chromosome-wide)
    cat(i,samples[i],":",nrow(out[[i]]),":",sum(out[[i]]$func=="nonsynonymous SNV"),":",
        sum(out[[i]]$func=="nonsynonymous SNV" & p.adjust(out[[i]]$somatic_p_value)<0.05),"\n")
}

## Generate summary information
aa<-c()
for(i in seq(samples)){
    aa<-rbind(aa,
              c(samples[i],nrow(out[[i]]),
                sum(out[[i]]$somatic_p_value<0.05),
                sum(p.adjust(out[[i]]$somatic_p_value,method='BH')<0.05),
                sum(p.adjust(out[[i]]$somatic_p_value,method='holm')<0.05),
                sum(out[[i]]$func=="nonsynonymous SNV"),
                sum(out[[i]]$func=="nonsynonymous SNV" & p.adjust(out[[i]]$somatic_p_value,method='BH')<0.05),
                sum(out[[i]]$func=="nonsynonymous SNV" & p.adjust(out[[i]]$somatic_p_value,method='holm')<0.05)))

}
aa<-as.data.frame(aa)
colnames(aa)<-c("sample","allSNV","raw_p","fdr_adj_p","fwer_adj_p","NonSynSNV","fdr_adj_p","fwer_adj_p")
## Print summary information to screen
print(aa)

## Combine per-chromosome data into single output file
outx<-c()
for(i in 1:length(out)) outx<-rbind(outx,out[[i]])

## Write out combined (all chromosomes) SNV data
write.csv(outx,file=paste0('annovar-snv-',dt,'.csv'),row.names=F)

##############################################################
## InDel data
##############################################################

outIndel<-list()
## No exonic indels for chr 2 3 6 7 9 12 14 15 16 17 18 19 X Y
## Only want chr 1 4 5 8 10 11 13
kp<-c(1,4,5,8,10,11,13)
## Check
cat(samples[kp])

## Loop through chromosomes
for(i in kp){

## Read in InDel data for this chromosome
    aa<-read.table(paste('../Varscan/',samples[i],'.indel',sep=''),header=T)

## Select only InDels called as somatic  
    bb<-aa[aa$somatic_status=="Somatic",]

## Read in annotation data for somatic InDels
    cc<-read.table(paste(samples[i],'.indel.SOMATIC.anno.variant_function',sep=''))
    dd<-read.table(paste(samples[i],'.indel.SOMATIC.anno.exonic_variant_function',sep=''),sep='\t')
    ee<-matrix("",nrow(cc),2)
    ln<-as.numeric(gsub("line","",dd[,1]))
    ee[ln,1]<-as.vector(dd[,2])
    ee[ln,2]<-as.vector(dd[,3])
    outIndel[[i]]<-cbind(rep(samples[i],nrow(cc)),cc[,1:2],ee,cc[,c(3:8,10)],bb[,c(5:15)],p.adjust(bb[,15],method='holm'),p.adjust(bb[,15],method='BH'),bb[,16:23])
    colnames(outIndel[[i]])[1:12]<-c("sample","type","gene/region","func","info","chr","start","end","ref","var","zygosity","readdepth")
    colnames(outIndel[[i]])[24:25]<-c("fwer_adj_p","fdr_adj_p")

## Output progress per chromosome
## - Total InDels
## - Non-synomynous InDels
## - Significant non-synomynous InDels (FDR adjusted: chromosome-wide) 
    cat(i,samples[i],":",nrow(outIndel[[i]]),":",sum(p.adjust(outIndel[[i]]$somatic_p_value,method="BH",)<0.05),
        ":",sum(p.adjust(outIndel[[i]]$somatic_p_value,method="holm")<0.05),"\n")
}

## Generate summary information
aa<-c()
for(i in kp){
    aa<-rbind(aa,
              c(samples[i],nrow(outIndel[[i]]),
                sum(out[[i]]$somatic_p_value<0.05),
                sum(p.adjust(outIndel[[i]]$somatic_p_value,method='BH')<0.05),
                sum(p.adjust(outIndel[[i]]$somatic_p_value,method='holm')<0.05)))
}
aa<-as.data.frame(aa)
colnames(aa)<-c("sample","allIndel","raw_p","fdr_adj_p","fwer_adj_p")

## print summary information to screen
print(aa)

## Write out combined (all chromosomes) InDel data
write.csv(outx,file=paste0('annovar-indel-',dt,'.csv'),row.names=F)

## Save session as an R data file
save.image(paste0('joinAnnotation-',dt,'.RData'))

