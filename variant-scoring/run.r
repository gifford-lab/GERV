#!/usr/bin/Rscript
runModel <- function(order){
	cmd = paste('./kmersnp_module',order,outputdir,seq_dir,snp_dir,model_parentdir,expt_name,genomefile,blfile_prefix,blfile_outfile_suffix,VT,VT_prefix,VT_suffix,kss_suffix,maxchr,maxk,K2,rc,resol,km_ver,SGE,win_resol,DLINK,DNUM_COV,kbeta,read.max,covmax,sep=" ")
	print(cmd)
	system(cmd)
}

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args)<2){
	stop('Not enough input arguments!')
}

paramfile = args[1]
order = args[2]

if (order=='score' || order=='STR_score'){
	if (length(args)!=3){
		stop('Please input the chr to score!')
	}
	chr_to_score = as.numeric(args[3])
}


#### Default Prams
rc=0;
km_ver=2;
win_resol=20;
DLINK = 'POISSON';
DNUM_COV = 0;
kbeta = 0;
covmax = 5;
SGE = 0;
covariate =c();

### load params
param.lines=readLines(paramfile)
param.lines=param.lines[nchar(param.lines)>0]
param.options=param.lines[grep('#',param.lines)]

for (option in 1:length(param.options)){
	ss=strsplit(param.options[option],'[# ]')[[1]]
	ss=ss[nchar(ss)>0]
	print(paste0('set:',ss[1],' to:',ss[2]))
	assign(ss[1],ss[2])
}

input.list = paste0(model_parentdir,'/',expt_name,'/input.list')

### load kmm input.list
input.lines=readLines(input.list)
input.lines=input.lines[nchar(input.lines)>0]
input.options=input.lines[grep('#',input.lines)]

epointers = (1:length(input.lines))[-grep('#|//',input.lines)]
ep2 = c(0,epointers)

flag = F
for( i in 1:length(epointers) ) {
	if(epointers[i] > (ep2[i]+1)){
		for(option in (ep2[i]+1):(epointers[i]-1)){
			ss=strsplit(input.lines[option],'[# ]')[[1]]
	    	ss=ss[nchar(ss)>0]
			print(paste0('set:',ss[1],' to:',ss[2]))
			assign(ss[1],ss[2])
		}
	}
	input.bams = input.lines[epointers[i]]
	rl = strsplit(input.bams,'[,]')[[1]]
	exptname = rl[1]
	if (exptname == expt_name){
		flag = T
		break
	}
}

### Induced params
blfile_tag = paste0('.',expt_name);
kss_suffix = paste0(blfile_tag,'.',snp_name);
blfile_prefix = paste0('fitbg',blfile_tag);
blfile_outfile_suffix = paste0(blfile_tag,'.baseline.bin');
#blfile_outfile_suffix = paste0(blfile_tag,'.reads.bin');
seq_dir = paste0('/seq_data/');
K2 = 2*as.numeric(k)
if (length(covariate)>0){
	if (covariate[1]!='none'){
		DNUM_COV = 1
		km_ver = 3
	}
}
if (km_ver==3){
	system('ln -s /kmm/master /kmm/delete_later')
}else{
	system('ln -s /kmm/nocov /kmm/delete_later')
}
### Run the script
if (order=='preprocess'){
	runModel(2)
	runModel(3)

}else{
	if (order=='score'){
		maxchr = chr_to_score;
		runModel(4)
	}else{
		if (order=='combine.result'){
			runModel(5)
		}else{
			if (order == 'STR_score'){
				maxchr = chr_to_score;
			 	runModel(7)
			}else{
				print('Order not recognized!')
			}
		}
	}
}








