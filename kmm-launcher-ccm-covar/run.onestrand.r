#!/usr/bin/Rscript

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

input.list = args[1]
cred.file = args[2]


#ec2 default settings
price = '3.0'
realm = 'us-east-1'
region = paste0(realm,'d')
ami='ami-d05e75b8'
itype='c3.8xlarge'

print('parse credential file')
cf=readLines(cred.file)
cf=cf[-grep('#',cf)]
for(sp in strsplit(cf,':')){
    print(sp)
    assign(sp[1],sp[2])
}
if(!file.exists(rsa_key)){print('check rsa key is readable')}
keyname=rev(strsplit(rsa_key,'[/.]')[[1]])[2]
print('setting key name to:')
print(keyname)

tmp = paste0(tempfile(),'.rsa')
file.copy(rsa_key,tmp)
Sys.chmod(tmp,mode='600')
rsa_key = tmp


#set up credentials
starcluster.rsa = rsa_key
#mailaddr = 'thashim@csail.mit.edu' # MAKE SURE THIS EMAIL IS SET UP ON AMAZON SES
key.name = key_name
access.key = access_key
secret.key = secret_key
bucket.name = bucket_name

system('mkdir ~/.aws')
system(paste0('printf \"[default]\naws_access_key_id=',access_key,'\naws_secret_access_key=',secret_key,'\" > ~/.aws/config'))

############
# default parameters (can be overwirtten in .list files

#organism specific
genome='hg19'

#runtime params
maxk=8
k=200
resol=1
cov.max = 1
covariate='none'
link='POISSON'

#runtime params (used)
read.max=50
smooth.window=1
require('utils')
cov.num = 0

#probably dont need changing..
branch = 'master'
gbase = '/genome/'


### load input
input.lines=readLines(input.list)
input.lines=input.lines[nchar(input.lines)>0]
input.options=input.lines[grep('#',input.lines)]
bam.prefix=''
out.prefix=''
postfix=''

epointers = (1:length(input.lines))[-grep('#|//',input.lines)]
ep2 = c(0,epointers)

############
# functions

rsystem <- function(sh,intern=F,wait=T){
    system(paste0('ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -i ',starcluster.rsa,' ubuntu@',INSTANCE_NAME,' ',shQuote(sh)),intern=intern,wait=wait)
}

scptoclus <- function(infile,out,intern=F){
    system(paste0('scp -C -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -c arcfour -i ',starcluster.rsa,' ',shQuote(infile),' ubuntu@',INSTANCE_NAME,':',shQuote(out)))
}

scpstring <- function(infile,out,intern=F){
    paste0('scp -C -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -c arcfour -i ',starcluster.rsa,' ',shQuote(infile),' ubuntu@',INSTANCE_NAME,':',shQuote(out))
}

scpfromclus <- function(infile,out,intern=F){
    system(paste0('scp -C -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -c arcfour -i ',starcluster.rsa,' -r ubuntu@',INSTANCE_NAME,':',shQuote(infile),' ',shQuote(out)))
}

getFilename<-function(x){
    rev(strsplit(x,'/')[[1]])[1]
}

############
# script


for( i in 1:length(epointers) ) {
    t1=Sys.time()

    if(epointers[i] > (ep2[i]+1)){
        for(option in (ep2[i]+1):(epointers[i]-1)){
            ss=strsplit(input.lines[option],'[# ]')[[1]]
            ss=ss[nchar(ss)>0]
            print(paste0('set:',ss[1],' to:',ss[2]))
            assign(ss[1],ss[2])
        }
    }
    #
    input.bams = input.lines[epointers[i]]
    rl = strsplit(input.bams,'[,]')[[1]]
    exptname = rl[1]
    bamlist = rl[-1]
    save.image(paste0('state',i,'.RData'))
}

runlist=sapply(1:length(epointers),function(i){
    paste0('/kmm/run.cluster.onestrand.r ',i)
})
writeLines(runlist,'/kmm/runlist.txt')
system('cat /kmm/runlist.txt | parallel')
