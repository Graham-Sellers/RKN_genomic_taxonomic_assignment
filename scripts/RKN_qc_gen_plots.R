#install.packages('circlize')
library('circlize')


path = '../qc_gen/kraken2/reports'

file_names = dir(path, pattern ='.txt')
file_names = gsub('.txt','',file_names[file_names != 'unclassified.txt'])
samples = list()
for(i in (file_names)){
  df = read.csv(paste(path,'/',i,'.txt',sep = ''),header=F, sep='\t')
  df = df[(grep('Meloidogyne', df[,6])),]
  df$taxlen = lengths(regmatches(df[,6], gregexpr("  ", df[,6])))
  df$prop = rev(cumsum(rev(df[,3]))/sum(df[,3]))
  df$taxlen= 1/(df$taxlen - (min(df$taxlen) -1))
  df$otus = trimws(df[,6])
  samples[[i]]=df
}
names(samples[4])

otu_list = c('Meloidogyne','Meloidogyne chitwoodi','Meloidogyne enterolobii','Meloidogyne floridensis','Meloidogyne graminicola', 'Meloidogyne hapla', 'Meloidogyne luci', 'Meloidogyne incognita group',
             'Meloidogyne arenaria', 'Meloidogyne javanica', 'Meloidogyne incognita')

shades = colorRampPalette(c('grey95','grey85','grey70','lightblue4','lightskyblue2','dodgerblue','royalblue3','purple','magenta','maroon','grey10'))(length(otu_list))
names(shades) = otu_list

#barplot(rep(1,length(otu_list)),space=0,col=shades)

# PLOT #

png('RKN_qc_gen_javanica03_plot.png',width=1600,height=1700,units='px',pointsize=30) # for RKN_lib2
#png('RKN_lib3_plot.png',width=1600,height=1600,units='px',pointsize=30) # for RKN_lib3

#LO <- matrix(c(1,2,3,4,5,6,7,8,13,13,9,10,13,13,11,12), nrow=4, ncol=4, byrow=TRUE) # for RKN_lib2
#LO <- matrix(c(1,2,3,4,11,11,5,6,11,11,7,8,11,11,9,10), nrow=4, ncol=4, byrow=TRUE) # for RKN_lib3

lo=c(5,10,15,20,
4,9,14,19,
3,8,13,18,
2,7,12,17,
1,6,11,16)

LO <- matrix(c(21,5,10,15,20,23,
               21,4,9,14,19,24,
               22,3,8,13,18,25,
               22,2,7,12,17,26,
               22,1,6,11,16,27,
               28,28,28,28,28,28), nrow=6, ncol=6, byrow=TRUE) # for RKN_lib2

WIDTHS <- c(1,2,2,2,2,3) #widths of each figure in layout (i.e. column widths)
HEIGHTS <- c(2,2,2,2,2,2.5)  #heights of each figure in layout (i.e. row heights)

layout(LO, heights=HEIGHTS, widths=WIDTHS)

par(mar=c(0.5,0,0.5,0))

used_otus = list()
for (s in 1:length(samples)){
  dfp=samples[[s]]
  used_otus = c(used_otus, dfp$otus)
  used_otus = unique(used_otus)
  plot(0,0,axes = F,cex = 0, ylab=NA,xlab=NA,asp = 1, main = NA, cex.main = 2.5)
  for (i in 1:dim(dfp)[1]){
    draw.sector(
      start.degree = 270,
      end.degree = -(360 * dfp[i,]$prop)+270,
      rou1 = 1,
      rou2 = ifelse(1 - dfp[i,]$taxlen == 0, 0.25, 1- dfp[i,]$taxlen),
      center = c(0, 0),
      clock.wise = T,
      col = shades[dfp[i,]$otus],
      border = "black",
      lwd = par("lwd"),
      lty = par("lty"))
  }
  text(-1.1,-1.1,paste('n =',sum(dfp[,3]),sep=' '),pos=4,font=3,cex=1.2,xpd=T)
}

#par(mar=c(0,0,0,0))

plot(0,0,xlim=c(0,1),ylim=c(0,2), cex=0,axes=F)
text(0.5,1, 'Kraken 2',cex=3, srt = 90)

lines(c(0.9,0.9),c(0.2,1.8),lwd=3)
plot(0,0,xlim=c(0,1),ylim=c(0,3), cex=0,axes=F)
text(0.5,1.5, 'fastp QC',cex=3, srt = 90)
lines(c(0.9,0.9),c(0.2,2.8),lwd=3)

plot(0,0,xlim=c(0,1),ylim=c(0,1), cex=0,axes=F)
text(0.5,0.5, 'Kraken 2:\nconfidence score 0.01\nminimum base quality 10\n+ full fastp QC',cex=1.5,font=3)
plot(0,0,xlim=c(0,1),ylim=c(0,1), cex=0,axes=F)
text(0.5,0.5, 'Kraken 2:\nminimum base quality 10\n+ full fastp QC',cex=1.5,font=3)
plot(0,0,xlim=c(0,1),ylim=c(0,1), cex=0,axes=F)
text(0.5,0.5, 'Full fastp QC:\ntrim leading 30 bases\nminimum base quality 15',cex=1.5,font=3)
plot(0,0,xlim=c(0,1),ylim=c(0,1), cex=0,axes=F)
text(0.5,0.5, 'fastp QC:\nminimum base quality 15',cex=1.5,font=3)
plot(0,0,xlim=c(0,1),ylim=c(0,1), cex=0,axes=F)
text(0.5,0.5, 'fastp QC:\nminimum base quality 0',cex=1.5,font=3)
plot(3,0,xlim=c(0.5,11.5),ylim=c(0,1), cex=0,axes=F)
text(2,0.9, '≥ 0 bp',cex=2)
text(4,0.9, '≥ 1000 bp',cex=2)
text(6,0.9, '≥ 2000 bp',cex=2)
text(8,0.9, '≥ 4000 bp',cex=2)
text(5,0.6, 'Read length',cex=3)

used_otus = otu_list[otu_list %in% used_otus]

shades[used_otus]
par(mar=c(0,0,0,0))

legend('bottom', legend=(gsub('Meloidogyne ','M. ', used_otus)),
       fill = shades[used_otus], cex=1.5, bty='n',
       text.font = 3, pt.cex = 10, ncol=4, y.intersp =1, x.intersp =1.5,
       pt.lwd=4)

dev.off()

#############

used_otus = otu_list[otu_list %in% used_otus]

shades[used_otus]
par(mar=c(0,0,0,0))

plot(0,0,xlim=c(0,2),ylim=c(0,3), cex=0, axes = F)
legend('left', legend=(gsub('Meloidogyne ','M. ', used_otus)),
       fill = shades[used_otus], cex=2.5, bty='n',
       text.font = 3, pt.cex = 10, y.intersp = 1.5,
       pt.lwd=4)

barplot(rep(1,length(used_otus)),space=0,col=shades[used_otus])
