#!/usr/bin/Rscript

options(echo=TRUE)
cms <- commandArgs(trailingOnly = TRUE)
print(cms)

load(paste0('state',cms[1],'.RData'))


print('Starting experiment:')
print(exptname)
print(bamlist)
outdir = file.path(outdir,paste0(exptname,postfix))

##organism specific config
genomedir=file.path(gbase,paste0(genome,'.in'))
offset.file = file.path(gbase,paste0(genome,'.offsets.txt'))
chroffs = as.double(readLines(offset.file))
chr.name = as.character(1:(length(chroffs)-1))
##
trainchr = which(chroffs > 2^31)[1]-1
testchr = which(chroffs > 2^31)[1]
##do NOT exceed 2147483647 (2^31-1), chr9 human is 1680373143
train.bases = min( chroffs[trainchr + 1], 2^31 - 1)
##heldout is chr18, this is the start of chr18 in your organism
heldout.start= chroffs[trainchr + 1]
##default 90702639 is the size of chr18 in mouse, keep it smaller than this number
test.bases = chroffs[testchr + 1] - heldout.start

#########
## process bam while we wait..

sirname = i
tmpdir = paste0('/tmp/kmer-',sirname,'/')
unlink(tmpdir,T,T)
dir.create(tmpdir)

for(bamfile in bamlist){
    if(!file.exists(file.path(bam.prefix,paste0(bamfile,'.bai')))){
        print("No bam index found, reindexing")
        x=file.path(bam.prefix,bamfile)
        y=paste0("samtools index ",shQuote(x))
        print(y)
        system(y,wait=T)
    }
}

print('Extracing reads from bam file')
dir.create(paste0(tmpdir,exptname),F)
unlink(paste0(tmpdir,exptname,'/*'),T,T)

for(bamfile in bamlist){
    print(bamfile)
    clist = sapply(1:length(chr.name),function(chr){
        paste0("(samtools view -F 788 -q ",quality," ",shQuote(file.path(bam.prefix,bamfile))," chr",chr.name[chr]," | cut -f 4 >> \'",tmpdir,exptname,"/allreads-",chr,".csv\'; touch ",tmpdir,exptname,"/chr",chr,".done)")
    })
    writeLines(clist,paste0(tmpdir,'commlist',cms[1],'.txt'))
    system(paste0('cat ',tmpdir,'commlist',cms[1],'.txt | parallel --progress -j 4'))
}

chrin = testchr
readsin=scan(paste0(tmpdir,exptname,'/allreads-',chrin,'.csv'),list(0))
rrle=rle(sort(readsin[[1]]))
rcoord=rrle$value
rnum=rrle$length
rbequiv=rep(0,test.bases)
print(max(rcoord))
rbequiv[rcoord[rcoord<test.bases]]=rnum[rcoord<test.bases]

con=file(paste0(tmpdir,exptname,'/heldout.in'),'wb')
writeBin(rbequiv,con,4)
close(con)

## If we have covariates, proc them here.
if(covariate!='none'){
    covbamlist = strsplit(covariate,'[,]')[[1]]
    tmpdir.cov = paste0('/tmp/kmer-',sirname,'-cov/')
    unlink(tmpdir.cov,T,T)
    dir.create(tmpdir.cov)
    for(bamfile in covbamlist){
        if(!file.exists(file.path(bam.prefix,paste0(bamfile,'.bai')))){
            print("No bam index found, reindexing")
            x=file.path(bam.prefix,bamfile)
            y=paste0("samtools index ",shQuote(x))
            print(y)
            system(y,wait=T)
        }
    }
    print('Extracing reads from bam file')
    dir.create(paste0(tmpdir.cov,exptname),F)
    unlink(paste0(tmpdir.cov,exptname,'/*'),T,T)
    for(bamfile in covbamlist){
        print(bamfile)
        clist = sapply(1:length(chr.name),function(chr){
            paste0("(samtools view -F 788 -q ",quality," ",shQuote(file.path(bam.prefix,bamfile))," chr",chr.name[chr]," | cut -f 4 >> \'",tmpdir.cov,exptname,"/allreads-",chr,".csv\'; touch ",tmpdir.cov,exptname,"/chr",chr,".done)")
        })
        writeLines(clist,paste0(tmpdir.cov,'commlist',cms[1],'.txt'))
        system(paste0('cat ',tmpdir.cov,'commlist',cms[1],'.txt | parallel --progress -j 4'))
    }
    chrin = testchr
    readsin=scan(paste0(tmpdir.cov,exptname,'/allreads-',chrin,'.csv'),list(0))
    rrle=rle(sort(readsin[[1]]))
    rcoord=rrle$value
    rnum=rrle$length
    rbequiv=rep(0,test.bases)
    print(max(rcoord))
    rbequiv[rcoord[rcoord<test.bases]]=rnum[rcoord<test.bases]
    con=file(paste0(tmpdir.cov,exptname,'/heldout.covariates.in'),'wb')
    writeBin(rbequiv,con,4)
    close(con)
}


####
# Transfer files.
inputdir = file.path(outdir,'input')
dir.create(inputdir,showWarnings=FALSE,recursive =TRUE)

scptoclus(genomedir,file.path(inputdir,'genome.in'))
scptoclus(offset.file,file.path(inputdir,'offsets.txt'))

clist=sapply(1:(length(chr.name)),function(chr){
    scpstring(paste0(tmpdir,exptname,'/allreads-',chr,'.csv'),inputdir)
})
writeLines(clist,paste0(tmpdir,'commlist',cms[1],'.txt'))
system(paste0('cat ',tmpdir,'commlist',cms[1],'.txt | parallel --progress -j 4'))
scptoclus(paste0(tmpdir,exptname,'/heldout.in'),inputdir)
## If we have covariates, transfer covar files here:
if(covariate!='none'){
	covardir = file.path(outdir,'covar')
	dir.create(covardir,showWarnings=FALSE,recursive =TRUE)
    clist=sapply(1:(length(chr.name)),function(chr){
        scpstring(paste0(tmpdir.cov,exptname,'/allreads-',chr,'.csv'),covardir)
    })
    writeLines(clist,paste0(tmpdir.cov,'commlist',cms[1],'.txt'))
    system(paste0('cat ',tmpdir.cov,'commlist',cms[1],'.txt | parallel --progress -j 4'))
    scptoclus(paste0(tmpdir.cov,exptname,'/heldout.covariates.in'),inputdir)
}

unlink(tmpdir,T,T)
if(covariate!='none')
    unlink(tmpdir.cov,T,T)

### run
kmmdir = file.path(outdir,'delete_later')
kmmbuild_dir = file.path(outdir,'delete_later/build')
real_outdir = file.path(outdir,'output')

dir.create(real_outdir,showWarnings=FALSE,recursive=TRUE)
#system('echo -e "Host github.mit.edu\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config')
if (!file.exists(kmmdir)){
	system(paste('git clone https://thashim-ro:*SybT2X9@bitbucket.org/ddkang/ccm_testing.git', kmmdir, sep=' '))
}
dir.create(kmmbuild_dir,showWarnings=FALSE)
system(paste0('rm -rf ',kmmbuild_dir,'/*'))

if(covariate!='none'){
    ncov = 1
    branch = 'glm_v2'
}else{
    ncov = 0
    kbeta = 0
    branch = 'no91'
}
system(paste('cd',kmmdir,'; git pull; git checkout ',branch,sep=' '))
cmakestr = paste0('-DK=',k,' -DRESOL=',resol,' -DKBIG=',maxk,' -DLINK=',link,' -DNUM_COV=',ncov,' -DK_BETA=',kbeta)
print(cmakestr)
system(paste('cd ',kmmbuild_dir,'; cmake .. ; make clean; make -j CXX_DEFINES=\"',cmakestr,'\"',sep=''))
##make reads
read.options=paste0(
    '--num_chr=',rev(chr.name)[1],' ',
    '--num_bases=',rev(chroffs)[1],' ',
    '--max_ct=100 ',
    '--offsets_file=',file.path(inputdir,'offsets.txt'),' ')

readstr=paste('cd ',kmmbuild_dir,'; ./reads ',read.options,' --out_file=',file.path(inputdir,'reads.in'),' --reads_dir=',inputdir,sep='')

if(covariate!='none'){
    covstr=paste('cd ',kmmbuild_dir,'; ./reads ',read.options,' --out_file=',file.path(inputdir,'covariates.in'),' --reads_dir=',covardir,sep='')
}else{
    covstr=''
}

scptoclus(input.list,outdir)

if(covariate!='none'){
    coption=paste('--covariates=',file.path(inputdir,'covariates.in'),' --heldout_covariates=',file.path(inputdir,'heldout.covariates.in'),' --cov_max=',cov.max,sep='')
	postproc_script='/kmm/postproc.r'
}else{
    coption='--cov_max=1'
	postproc_script='/kmm/postproc_covar.r'
}

runstr=paste(file.path(kmmbuild_dir,'mpi_motif'),' --out_dir=',real_outdir,' --genome=',file.path(inputdir,'genome.in'),' --reads=',file.path(inputdir,'reads.in'),' --num_bases=',train.bases,' --read_max=',read.max,' --smooth_window_size=',smooth.window,' --heldout_start=',heldout.start,' --heldout_size=',test.bases,' --heldout_reads=',file.path(inputdir,'heldout.in'),' ', coption,' 2>&1 | tee ',file.path(outdir,'runlog.txt'),sep='')
print(runstr)

rl=readLines('/kmm/standalone.template.local.txt')
rl=gsub('READ_STR',readstr,rl)
rl=gsub('COV_STR',covstr,rl)
rl=gsub('RUN_STR',runstr,rl)
rl=gsub('TEST_CHR',testchr,rl)
rl=gsub('KMMDIR',kmmdir,rl)
rl=gsub('BUILD_DIR',kmmbuild_dir,rl)
rl=gsub('INPUTDIR',inputdir,rl)
rl=gsub('REALOUTDIR',real_outdir,rl)
rl=gsub('OUTDIR',outdir,rl)
rl=gsub('POSTPROC_SCRIPT',postproc_script,rl)

runall = file.path(outdir,'runall.sh')
system(paste('printf \'',paste0(rl,collapse='\n'),'\' >',runall,sep=' '))
system(paste('chmod +x',runall,sep=' '))
system(paste('nohup',runall,' `</dev/null` >',file.path(outdir,'nohup.txt'),' 2>&1 ',sep=' '))
#print('Launch finished in:')
#print(Sys.time()-t1)
