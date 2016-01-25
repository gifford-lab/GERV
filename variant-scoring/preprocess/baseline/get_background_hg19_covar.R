args = commandArgs(T)
paramdir = args[1]
source(paste0(paramdir,"/preprocess/baseline/get_background_params.R"))
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

    outdir = file.path(OUT_DIR,ename);
	system(paste0('lzop -d ',outdir,'/*.lzo'))
	system(paste0('lzop -d ',indir,'/covariates.in.lzo'))
     
    readsize=19000000
    xszall=sum(4^(0:kmax))*K2
     
    runset=sort(as.double(sapply(strsplit(list.files(outdir,pattern=param1),'[_.]'),function(i){i[2]})))
    metapars=sapply(runset,function(i){
        print(i)
        c(readBin(paste0(outdir,'/eta_',i,'.bin'),double(),size=4,n=1),readBin(paste0(outdir,'/',param2,'_',i,'.bin'),double(),size=4,n=1),readBin(paste0(outdir,'/heldout_',i,'.bin'),double(),size=4,n=1))#,activebase,activek)
    })
     
    bestreg=runset[which.min(metapars[3,])]
	system(paste0('cp ',outdir,'/',param2,'_',bestreg,'.bin ',x0filename))
	x0data = readBin(paste0(outdir,'/',param2,'_',bestreg,'.bin'),double(),size=4,n=2);
	if (length(x0data)==2){
		x0nt = x0data[1];
		x0 = x0data[2];
	}else{
		x0nt = x0data;
		x0 = x0nt;
	}
	indir = paste0(OUT_DIR,ename_in)
	covfile = paste0(indir,'/covariates.in')
	system(paste0('lzop -d ',covfile,'.lzo'))
	betafile = paste0(outdir,'/beta_',bestreg,'.bin')

setwd('/kmm/delete_later/build/')
cmd = (paste0('sh -c \"make -j CXX_DEFINES=\\\"-DK=',K,' -DRESOL=',RESOLUTION,' -DKBIG=',kmax,' -DK_BETA=',DK_BETA,' -DNUM_COV=',DNUM_COV,' -DLINK=',DLINK,'\\\" \"'))
print(cmd)
system(cmd)

cmd =(paste0('sh -c \"./validation  --output=',bg.file,' --genome=',genomefile,' --xopt=',outdir,'/',param1,'_',bestreg,'.bin --start=',genomestart,' --end=',genomeend,' --x0=',x0,' --x0nt=',x0nt,' --covariates=',covfile,' --beta=',betafile,'\"'))
print(cmd)
system(cmd)
#system(paste0(validationfile,' --output=',bg.file,' --genome=',genomefile,' --xopt=',outdir,'/yall_',bestreg,'.bin',' --start=0 --end=',genomesize_midpoint,' --x0=',0),intern=T)

