args = commandArgs(T)
inputdir = args[1]
outdir = args[2]
sumdir = args[3]
runlogfile = args[4]
builddir = args[5]
chrheld = as.numeric(args[6])

old.repos <- getOption("repos")
on.exit(options(repos = old.repos))
new.repos <- old.repos
new.repos["CRAN"] <- "http://cran.stat.ucla.edu"
options(repos = new.repos)

install.packages('Rcpp')
install.packages('inline')
install.packages('bitops')
install.packages('gplots')
install.packages('RColorBrewer')

require('Rcpp');require('inline')
uniconv ='
Rcpp::NumericVector cvadd(vadd);
Rcpp::NumericVector cvret(cvadd.size());
int cksize = as<int>(ksize);
int halfsz = (cksize-1)/2;
double temp=0;
for(int i=0; i < (halfsz);i++){
 temp+=cvadd[i];
}
for(int i=0; i < cvadd.size();i++){
 if((i+halfsz)< cvadd.size()){
  temp+=cvadd[i+halfsz];
 }
 cvret[i]=temp/cksize;
 if((i-halfsz)>=0){
  temp-=cvadd[i-halfsz];
 }
}
return cvret;
'
uconv <- cxxfunction(signature(vadd="numeric",ksize="integer"),uniconv,plugin="Rcpp",includes="#include <numeric>")

runset=sort(as.double(sapply(strsplit(list.files(outdir,pattern='yall'),'[_.]'),function(i){i[2]})))
metapars=sapply(runset,function(i){
        print(i)
        c(readBin(file.path(outdir,paste0('eta_',i,'.bin')),double(),size=4,n=1),readBin(file.path(outdir,paste0('x0_',i,'.bin')),double(),size=4,n=2),readBin(file.path(outdir,paste0('heldout_',i,'.bin')),double(),size=4,n=1))#,activebase,activek)
    })

bn = which.min(metapars[4,])
 bestreg=runset[bn]
bx0nt = metapars[2,bn]
bx0 = metapars[3,bn]

runlog=readLines(runlogfile)
print(runlog[1:20])

kmax=as.integer(strsplit(runlog[c(grep('Kmax',runlog),grep('KBIG',runlog))],' ',fixed=T)[[1]][2])
k2=as.integer(strsplit(runlog[grep('K2',runlog)],' ',fixed=T)[[1]][2])
resol=as.integer(strsplit(runlog[grep('RESOL',runlog)],' ',fixed=T)[[1]][2])

xszall=(sum(4^(0:kmax)))*k2/resol

dir.create(sumdir)
pdf(paste0(sumdir,'/heldout.pdf'))
plot(log(metapars[1,],10),metapars[4,])
dev.off()

allbins=t(matrix(readBin(file.path(outdir,paste0('yall_',bestreg,'.bin')),double(),size=4,n=xszall),nrow=k2/resol))[-1,]

require(bitops)
bases=c('A','T','G','C')
ntokmer<-function(num,kl){
    kmer=''
    for(i in 1:kl){
        kmer=paste0(kmer,bases[bitAnd(num,3)+1])
        num=bitShiftR(num,2)
    }
    kmer
}

kmernames=do.call(c,lapply(1:kmax,function(i){
    sapply((1:(4^(i)))-1,ntokmer,kl=i)
}))

subs=1
of=order(-rowSums(abs(allbins)))[1:min(1500,nrow(allbins))]
qtl=quantile(unlist(allbins[of,]),c(0.001,0.999))*2
pdf(paste0(sumdir,'/profile.pdf'))
kseq = seq((-k2/2),(k2/2-1),by=resol)
for(i in of){
    plot(kseq,unlist(allbins[i,]),type='l',main=kmernames[i],ylim=qtl)
}
dev.off()

source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings",ask=F,suppressUpdates=T,suppressAutoUpdate=T)

require('Biostrings')
matall = allbins

sd=as.matrix(stringDist(kmernames[of]))
sd2=as.matrix(stringDist(c(kmernames[of],as.character(reverseComplement(DNAStringSet(kmernames[of]))))))[1:(nrow(sd)),nrow(sd)+(1:(nrow(sd)))]
sd[sd2<sd] = sd2[sd2<sd]
abd = as.matrix(dist(allbins[of,]))

hcc=hclust(as.dist(sd + abd/mean(abd[sd==1])),method='complete')

hccut=cutree(hcc,h=2)
range.y = range(matall[of,])
pdf(paste0(sumdir,'/hctest-dist2.pdf'))
for(i in unique(hccut)){
    wi=which(hccut==i)
    rb=rainbow(min(length(wi),10))
    #range.y = range(matall[of,][wi,])
    plot(unlist(matall[of,][wi[1],]),type='l',ylim=range.y)
    if(length(wi)>1){
        for(j in 2:min(length(wi),10))
            points(unlist(matall[of,][wi[j],]),type='l',col=rb[j])
    }
    legend('topleft',legend=kmernames[of][wi],lwd=1,col=c('black',rb[2:length(wi)]))
}
dev.off()

## heatmappy thing
require('gplots')
require('RColorBrewer')
ofsub = 1:min(500,nrow(sd))
hccsub=hclust(as.dist(sd[ofsub,ofsub] + abd[ofsub,ofsub]/mean(abd[sd==1])),method='complete')

pal=rev(colorRampPalette(brewer.pal(11,'RdBu'))(100))

pdf(paste0(sumdir,'/heatplot.pdf'),7,35)
heatmap.2(allbins[of[ofsub],],Rowv=as.dendrogram(hccsub),Colv=FALSE,scale='none',trace='none',dendrogram='row',labRow=kmernames[of[ofsub]],key=F, col=pal, symbreaks=T,labCol='', lhei = c(0.1,10))
dev.off()



offsets = as.double(readLines(file.path(inputdir,'offsets.txt')))
chrlen=diff(offsets)

i=18
unlink('/tmp/*.done')
for(i in 1:length(chrlen)){
    print(i)
    stc = offsets[i]
    edc = offsets[i+1]
    system(paste0(file.path(builddir,'validation'),' --output=',sumdir,'/fitted.',i,'.bin --genome=',file.path(inputdir,'genome.in'),' --xopt=\'',outdir,'/yall_',bestreg,'.bin\' --start=',stc+1,' --end=',edc,' --x0=',bx0,' --x0nt=',bx0nt, ' --covariates=',file.path(inputdir,'covariates.in'),' --beta=',file.path(outdir,paste0('beta_',bestreg,'.bin')),' && touch /tmp/',i,'.done'))
}
while(length(list.files('/tmp','done'))<length(chrlen)){Sys.sleep(5)}

readsize=chrlen[chrheld]

con <- file(file.path(inputdir,'reads.in'),open='rb')
seek(con,where = offsets[chrheld] * 4)
rbequiv=readBin(con,double(),size=4,n=readsize)
rbequiv[rbequiv>5]=5
close(con)

#rbequiv[rbequiv>1]=1
rbhat=readBin(file.path(sumdir,paste0('fitted.',chrheld,'.bin')),double(),size=4,n=readsize)
rbhat[rbhat>5]=5

subr = 1:(readsize-1)

wsz=2000
ea=uconv(rbhat[subr],wsz)
eq=uconv(rbequiv[subr],wsz)
c.signal.raw=cor(ea,eq)
c.signal=cor(log(ea+1),log(eq+1))
c.r.signal=cor(sqrt(ea),sqrt(eq))
corvec=c(c.signal.raw,c.signal,c.r.signal)
print(corvec)

writeLines(format(corvec),paste0(sumdir,'/cor.txt'))

ssub = sample(1:length(subr),200000)
irs = isoreg(ea[ssub],sqrt(eq[ssub]))
afun=approxfun(sort(irs$x),irs$yf,rule=0)
cor(afun(ea),sqrt(eq))
cor(ea,eq)


set.seed(1)
yl=c(0,mean(eq)+6*sd(eq))
pdf(paste0(sumdir,'/detail.pdf'),20,7)
for(j in 1:200){
        ssub=sample(50000:(length(subr)-50000),1)
            subhat=sort(sample(ssub+(-50000:50000),10000))
            plot(subhat,ea[subhat],type='l',ylim=yl,col=rgb(1,0,0),lwd=1,xlab='coordinate',ylab='reads')
            points(subhat,eq[subhat],type='l',col=rgb(0,1,0),lwd=1)
    }
dev.off()

####
# call stuff

set.seed(1)
eqsub = sort(sample(which(eq > mean(eq)+sd(eq)*2),size=10000,replace=T))

pdf(paste0(sumdir,'/err.pdf'))
plot(ea[eqsub],eq[eqsub]+1e-5,pch='.',xlab='predicted',ylab='observed',log='xy')
dev.off()

sdfact=6
eqcall=which(eq > mean(eq)+sd(eq)*sdfact)

require(IRanges)
callrs=reduce(flank(IRanges(start=eqcall,width=1),width=1000,both=T))
callord=order(-ea[(start(callrs)+end(callrs))/2])
yl = c(0,mean(eq)+sd(eq)*sdfact*2)#range(c(ea[eqcall],eq[eqcall]))
set.seed(1)
pdf(paste0(sumdir,'/callplots.pdf'),20,7)
par(mfrow=c(1,3))
for(j in callord){
        subhat=sort(sample((start(callrs[j])-1000):(end(callrs[j])+1000),size=2000))
            plot(subhat,ea[subhat],type='l',ylim=yl,col=rgb(1,0,0),lwd=1,xlab='coordinate',ylab='reads')
            points(subhat,eq[subhat],type='l',col=rgb(0,1,0),lwd=1)
            points(subhat,rbhat[subr][subhat],type='l',col=rgb(1,0,0),lwd=1,lty=2)
            abline(h=mean(ea)+sd(ea)*sdfact,col='red')
            rug(subhat[which(rbequiv[subr[subhat]]>0)])
    }
dev.off()

eacall=which(ea > mean(ea)+sd(ea)*sdfact)
callea=reduce(flank(IRanges(start=eacall,width=1),width=1000,both=T))
fos=findOverlaps(callea,callrs)
coverage.kmer=length(fos)/length(callrs)
ppv.kmer=length(fos)/length(callea)

hits=c(coverage.kmer,ppv.kmer)
print(hits)
writeLines(format(hits),paste0(sumdir,'/hits.txt'))
