args <- commandArgs(T)
paramdir = args[1]
source(paste0(paramdir,'/preprocess/kmer_mat/prepare_kmer_mat_params.R'))
require(bitops)
library(R.matlab)
bases=c('A','T','G','C')
ntokmer<-function(num,kl){
    kmer=''
    for(i in 1:kl){
        kmer=paste0(kmer,bases[bitAnd(num,3)+1])
        num=bitShiftR(num,2)
    }
    kmer
}	

    procfiles=list.files(OUT_DIR)
    procdirs=file.info(paste0(OUT_DIR,procfiles))$isdir
    outnames=procfiles[procdirs]
     
    outdir = paste0(OUT_DIR,ename)
     
    readsize=19000000
    xszall=sum(4^(0:kmax))*K2/RESOL
     
    runset=sort(as.double(sapply(strsplit(list.files(outdir,pattern=param1),'[_.]'),function(i){i[2]})))
    metapars=sapply(runset,function(i){
        print(i)
        c(readBin(paste0(outdir,'/eta_',i,'.bin'),double(),size=4,n=1),readBin(paste0(outdir,'/',param2,'_',i,'.bin'),double(),size=4,n=1),readBin(paste0(outdir,'/heldout_',i,'.bin'),double(),size=4,n=1))#,activebase,activek)
    })
     
    bestreg=runset[which.min(metapars[3,])]
     
    allbins=lapply(bestreg:bestreg,function(j){
        readBin(paste0(outdir,'/',param1,'_',j,'.bin'),double(),size=4,n=xszall)
    })
     
    kmernames=do.call(c,lapply(1:kmax,function(i){
        sapply((1:(4^(i)))-1,ntokmer,kl=i)
    }))
     
    dmall=lapply(allbins,function(rb2){
        rbmat=t(matrix(rb2,nrow=(K2/RESOL)))[-1,]
        data.frame(kmernames,rbmat)
    })

	mat = dmall[[1]];
	mat = mat[,2:ncol(mat)];
	writeMat(outfile,mat=as.matrix(mat));
