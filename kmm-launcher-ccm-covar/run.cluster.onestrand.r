#!/usr/bin/Rscript

options(echo=TRUE)
cms <- commandArgs(trailingOnly = TRUE)
print(cms)

load(paste0('state',cms[1],'.RData'))


print('Starting experiment:')
print(exptname)
print(bamlist)

##organism specific config
genomedir=paste0(gbase,genome,'.in')
offset.file = paste0(gbase,genome,'.offsets.txt')
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

##launch instance
userdatablob=paste0(system('cat /kmm/user-data.txt | base64',intern=T),collapse='')
lspec = paste0("\'{\"UserData\":\"",userdatablob,"\",\"ImageId\":\"",ami,"\",\"KeyName\":\"",keyname,"\",\"InstanceType\":\"",itype,"\"}\'")

awscomm="aws"
if(price != "Inf"){
    launch=system(paste0(awscomm,' --region ',realm,' --output text ec2 request-spot-instances --spot-price ',price,' --launch-specification ',lspec),intern=T)
    sirname = strsplit(launch,'\t')[[1]][4]
    sistatus = ''
}else{
    launch=system(paste0(awscomm,' --region ',realm,' --output text ec2 run-instances --image-id ',ami,' --key-name ',keyname,' --instance-type ',itype,' --user-data file:///kmm/user-data.txt'),intern=T)
    iname=strsplit(launch,'\t')[[3]][8]
    sirname=iname
}

#########
## process bam while we wait..

tmpdir = paste0('/tmp/kmer-',sirname,'/')
unlink(tmpdir,T,T)
dir.create(tmpdir)

for(bamfile in bamlist){
    if(!file.exists(paste0(bam.prefix,bamfile,'.bai'))){
        print("No bam index found, reindexing")
        x=paste0(bam.prefix,bamfile)
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
        paste0("(samtools view -F 788 -q ",quality," ",shQuote(paste0(bam.prefix,bamfile))," chr",chr.name[chr]," | cut -f 4 >> \'",tmpdir,exptname,"/allreads-",chr,".csv\'; touch ",tmpdir,exptname,"/chr",chr,".done)")
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
        if(!file.exists(paste0(bam.prefix,bamfile,'.bai'))){
            print("No bam index found, reindexing")
            x=paste0(bam.prefix,bamfile)
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
            paste0("(samtools view -F 788 -q ",quality," ",shQuote(paste0(bam.prefix,bamfile))," chr",chr.name[chr]," | cut -f 4 >> \'",tmpdir.cov,exptname,"/allreads-",chr,".csv\'; touch ",tmpdir.cov,exptname,"/chr",chr,".done)")
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

#####
## check if spot is up

if(price!="Inf"){
print('wait for spot fulfilment')
while(sistatus!='fulfilled'){
    cat('.')
    tryCatch({
        sitest=system(paste0(awscomm,' --region ',realm,' --output text ec2 describe-spot-instance-requests --spot-instance-request-ids ',sirname),intern=T)
        sistatus=strsplit(sitest,'\t')[[grep('STATUS',sitest)]][2]
    }, error = function(e){
        print(e)
        Sys.sleep(10)
    })
    Sys.sleep(5)
}
iname = strsplit(sitest,'\t')[[1]][3]
}

rname=paste0(exptname,postfix)
system(paste0(awscomm,' --region ',realm,' --output text ec2 create-tags --resources ',iname,' --tags Key=Name,Value=',rname))

istatus = 'initializing'

checks.passed = 0

print('wait for status checks')
while(checks.passed < 2){
    cat('.')
    tryCatch({
        itest=system(paste0(awscomm,' --region ',realm,' --output text ec2 describe-instance-status --instance-ids ',iname),intern=T)
        if(length(itest)>=3){
            checks.passed = length(grep('passed',itest))
        }
    }, error = function(e){
        print(e)
        Sys.sleep(10)
    })
    Sys.sleep(5)
}

istat=system(paste0(awscomm,' --region ',realm,' --output text ec2 describe-instances --instance-ids ',iname),intern=T)

isub=strsplit(istat,'\t')[[grep('INSTANCES',istat)]]
INSTANCE_NAME = isub[grep('amazonaws',isub)]

while(length(grep('done',rsystem('ls /mnt',intern=T)))==0) { Sys.sleep(5) }


####
# Transfer files.

scptoclus(genomedir,'/mnt/input/genome.in')

scptoclus(offset.file,'/mnt/input/offsets.txt')

clist=sapply(1:(length(chr.name)),function(chr){
    scpstring(paste0(tmpdir,exptname,'/allreads-',chr,'.csv'),'/mnt/input/')
})
writeLines(clist,paste0(tmpdir,'commlist',cms[1],'.txt'))
system(paste0('cat ',tmpdir,'commlist',cms[1],'.txt | parallel --progress -j 4'))
scptoclus(paste0(tmpdir,exptname,'/heldout.in'),'/mnt/input/heldout.in')
## If we have covariates, transfer covar files here:
if(covariate!='none'){
    clist=sapply(1:(length(chr.name)),function(chr){
        scpstring(paste0(tmpdir.cov,exptname,'/allreads-',chr,'.csv'),'/mnt/covar/')
    })
    writeLines(clist,paste0(tmpdir.cov,'commlist',cms[1],'.txt'))
    system(paste0('cat ',tmpdir.cov,'commlist',cms[1],'.txt | parallel --progress -j 4'))
    scptoclus(paste0(tmpdir.cov,exptname,'/heldout.covariates.in'),'/mnt/input/heldout.covariates.in')
}

while(length(grep('setup',rsystem('ls /home/ubuntu/',intern=T)))==0) {
    Sys.sleep(5)
}

unlink(tmpdir,T,T)
if(covariate!='none')
    unlink(tmpdir.cov,T,T)

### run
scptoclus(rsa_key,'~/.ssh/id_rsa')
rsystem('echo -e "Host github.mit.edu\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config')
#rsystem('git clone git@github.mit.edu:thashim/ccm-devel.git /home/ubuntu/delete_later')
rsystem('git clone https://thashim-ro:*SybT2X9@bitbucket.org/ddkang/ccm_testing.git /home/ubuntu/delete_later')
rsystem('mkdir /home/ubuntu/delete_later/build')
rsystem('rm -rf ~/delete_later/build/*')
rsystem(paste0('cd ~/delete_later/; git pull; git checkout ',branch))

if(covariate!='none'){
    ncov = 1
}else{
    ncov = 0
    kbeta = 0
}
cmakestr = paste0('-DK=',k,' -DRESOL=',resol,' -DKBIG=',maxk,' -DLINK=',link,' -DNUM_COV=',ncov,' -DK_BETA=',kbeta)
rsystem(paste0('cd ~/delete_later/build/; cmake .. ; make clean; make -j CXX_DEFINES=\"',cmakestr,'\"'))
##make reads
read.options=paste0(
    '--num_chr=',rev(chr.name)[1],' ',
    '--num_bases=',rev(chroffs)[1],' ',
    '--max_ct=100 ',
    '--offsets_file=/mnt/input/offsets.txt ')

readstr=paste0('cd ~/delete_later/build/; ./reads ',read.options,' --out_file=/mnt/input/reads.in --reads_dir=/mnt/input/')

if(covariate!='none'){
    covstr=paste0('cd ~/delete_later/build/; ./reads ',read.options,' --out_file=/mnt/input/covariates.in --reads_dir=/mnt/covar/')
}else{
    covstr=''
}

scptoclus(input.list,'~/input.list')
rsystem(paste0('printf \'',paste0('i=',i,'\n'),'\' > ~/params.txt'));

##enable credentials remotely
rsystem('mkdir ~/.aws')
rsystem(paste0('printf \"[default]\naws_access_key_id=',access.key,'\naws_secret_access_key=',secret.key,'\" > ~/.aws/config'))


if(covariate!='none'){
    coption=paste0('--covariates=/mnt/input/covariates.in --heldout_covariates=/mnt/input/heldout.covariates.in --cov_max=',cov.max)
}else{
    coption='--cov_max=1'
}

runstr=paste0('~/delete_later/build/mpi_motif --out_dir=/mnt/output --genome=/mnt/input/genome.in --reads=/mnt/input/reads.in --num_bases=',train.bases,' --read_max=',read.max,' --smooth_window_size=',smooth.window,' --heldout_start=',heldout.start,' --heldout_size=',test.bases,' --heldout_reads=/mnt/input/heldout.in ',coption,' 2>&1 | tee /home/ubuntu/runlog.txt')

rl=readLines('/kmm/standalone.template.txt')
rname=paste0(exptname,postfix)
rl=gsub('READ_STR',readstr,rl)
rl=gsub('COV_STR',covstr,rl)
rl=gsub('RUN_NAME',rname,rl)
rl=gsub('REGION',realm,rl)
rl=gsub('SIRNAME',sirname,rl)
rl=gsub('INAME',iname,rl)
rl=gsub('EMAIL',mailaddr,rl)
rl=gsub('RUN_STR',runstr,rl)
rl=gsub('TEST_CHR',testchr,rl)
rl=gsub('BUCKET_NAME',bucket_name,rl)

rsystem(paste0('printf \'',paste0(rl,collapse='\n'),'\' > ~/runall.sh'))
rsystem('chmod +x ~/runall.sh')
rsystem('nohup ~/runall.sh `</dev/null` >nohup.txt 2>&1 &')
print('Launch finished in:')
print(Sys.time()-t1)
