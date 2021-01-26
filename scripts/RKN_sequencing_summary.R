#install.packages('dplyr')
library(dplyr)
getwd()

path = '../data/libraries/RKN_lib2'

seq_summary=read.csv(paste(path,'/sequencing_summary.txt',sep=''), header = T,sep ='\t')

barcodes = unique(sort(seq_summary$barcode_arrangement[seq_summary$barcode_arrangement != 'unclassified']))
shades = c('grey95','grey80','lightskyblue3')

png('read_hist_plot.png',width=1600,height=2000,units='px',pointsize=30)

LO <- matrix(c(1,1,1,2,2,2,3,4,5,6,7,8,9,10,11,12,13,14), nrow=6, ncol=3, byrow=TRUE) # for RKN_lib2

WIDTHS <- c(1,1,1) #widths of each figure in layout (i.e. column widths)
HEIGHTS <- c(3,1,3,3,3,3)  #heights of each figure in layout (i.e. row heights)

layout(LO, heights=HEIGHTS, widths=WIDTHS)

par(mar=c(4,5,4,1))

seq_len=seq_summary$sequence_length_template
filt_len=seq_summary$sequence_length_template[seq_summary$passes_filtering == TRUE]
barcoded=seq_summary$sequence_length_template[seq_summary$barcode_arrangement != 'unclassified']
gbp=round(sum(seq_summary$sequence_length_template)/1e9,2)

hist(seq_len,xlim=c(0,10000),
     breaks=((round((max(seq_len,na.rm=T)-min(seq_len,na.rm=T)+500),-3)/200)+1),
     xlab='Read length (kb)',ylab=NA,axes=F,cex.lab=1.5,main=paste('Sequenced reads (',gbp,' Gbp)',sep=''),cex.main=2,col=shades[1])
axis(1,seq(0,14000,1000),labels=seq(0,14,1),cex.axis=1.2,tck=-0.02)
axis(2,seq(0,16000,2000),labels=seq(0,16000,2000),cex.axis=1.2,tck=-0.02,las=1)
hist(filt_len,xlim=c(0,10000),
     breaks=((round((max(filt_len,na.rm=T)-min(filt_len,na.rm=T)+500),-3)/200)+1),
     xlab=NA,ylab=NA,axes=F,cex.lab=1.2,main=NA,cex.main=2,add=T,col=shades[2])
hist(barcoded,xlim=c(0,10000),
     breaks=((round((max(barcoded,na.rm=T)-min(barcoded,na.rm=T)+500),-3)/200)+1),
     xlab=NA,ylab=NA,axes=F,cex.lab=1.2,main=NA,cex.main=2,add=T,col=shades[3])
legend('topright',bty='n',cex=1.2,
       legend=c(length(seq_len),
                length(filt_len),
                length(barcoded)),fill=shades[1:3],title='reads')

par(mar=c(0,3,0,3))
plot(1,1,cex=0,axes=F)
legend('center', legend=c('raw reads', 'pass filter', 'barcoded'),
       fill = shades, cex=1.5, bty='n',
       text.font = 3, pt.cex = 10,
       pt.lwd=4, ncol=3)

par(mar=c(4,5,4,1))

for (i in barcodes){diff=(seq_summary[seq_summary$barcode_arrangement==i,])

  
  dat=diff$sequence_length_template
  hist(dat,xlim=c(0,10000),
       breaks=((round((max(dat,na.rm=T)-min(dat,na.rm=T)+500),-3)/500)+1),
       xlab='Read length (kb)',ylab=NA,axes=F,cex.lab=1.5,main=i,cex.main=2,col=shades[3])
  axis(1,seq(0,14000,2000),labels=seq(0,14,2),cex.axis=1.2,tck=-0.02)
  axis(2,seq(0,5000,200),labels=seq(0,5000,200),cex.axis=1.2,tck=-0.02,las=1)
  legend('topright',bty='n',cex=1.2,
         legend=length(dat),
         fill=shades[3],title='reads')}

dev.off()

############################
qc_reads = read.csv(paste('../RKN_lib2_readlengths.tsv',sep=''), header = T,row.names=1,sep ='\t')

shades = colorRampPalette(c('grey95','grey80','lightskyblue3','dodgerblue','purple','magenta'))(6)
non_spp = c(654580,189290)

krk_reads=read.csv(paste('../results/kraken2/outputs/', i, '.krk',sep=''), sep='\t', header = F)
krk_class = krk_reads[krk_reads[1] == 'C',]
mel_class = krk_class[(grep('Meloidogyne', krk_class[,3])),]
spp_class = mel_class[!grepl(paste(non_spp,collapse="|"),mel_class[,3]),]
inc_class= spp_class[(grep('6306',spp_class[,3])),]


qc_dat=(qc_reads[,i])

hist(qc_dat,xlim=c(0,10000),
     breaks=((round((max(qc_dat,na.rm=T)-min(qc_dat,na.rm=T)+500),-3)/500)+1),
     xlab=NA,ylab=NA,axes=F,cex.lab=1.2,main=NA,add=T,col=shades[4])
hist(mel_class[,4],xlim=c(0,10000),
     breaks=((round((max(mel_class[,4],na.rm=T)-min(mel_class[,4],na.rm=T)+500),-3)/500)+1),
     xlab=NA,ylab=NA,axes=F,cex.lab=1.2,main=NA,add=T,col=shades[5])
hist(inc_class[,4],xlim=c(0,10000),
     breaks=((round((max(inc_class[,4],na.rm=T)-min(inc_class[,4],na.rm=T)+500),-3)/500)+1),
     xlab=NA,ylab=NA,axes=F,cex.lab=1.2,main=NA,add=T,col=shades[6])
legend('topright',bty='n',cex=1.2,
       legend=c(length(dat),sum(qc_dat != 'NA', na.rm=T),length(mel_class[,4]),length(inc_class[,4])),
       fill=shades[3:6],title='reads')


par(mar=c(0,3,0,3))
plot(1,1,cex=0,axes=F)
legend('center', legend=c('raw reads', 'pass filter', 'barcoded', 'pass qc','Meloidogyne','M. incognita'),
       fill = shades, cex=1.5, bty='n',
       text.font = 3, pt.cex = 10,
       pt.lwd=4, ncol=3)
round(sum(seq_summary$sequence_length_template)/1000000000,2)
names(seq_summary)
