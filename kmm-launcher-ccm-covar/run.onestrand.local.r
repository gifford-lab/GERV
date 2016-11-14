#!/usr/bin/Rscript

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

input.list = args[1]

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
scptoclus <- function(infile,out,intern=F){
	system(paste('cp',infile,out,sep=' '))
}

scpstring <- function(infile,out,intern=F){
	paste('cp',infile,out,sep=' ')
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
    paste0('/kmm/run.cluster.onestrand.local.r ',i,' > log_',format(Sys.time(), "%h_%d_%R_%Y"),'_',i,' 2>&1')
})
writeLines(runlist,'/kmm/runlist.txt')
system('cat /kmm/runlist.txt | parallel')
