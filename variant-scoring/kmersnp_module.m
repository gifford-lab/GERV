function [] = kmersnp_module(order,topdir,seq_dir,snp_dir,kmer_binary_dir,kmer_binary_folder,genomefile,blfile_prefix,blfile_outfile_suffix,VT,VT_prefix,VT_suffix,kss_suffix,maxchr,kmax,K2,rc,RESOL,km_ver,SGE,win_resol,DLINK,DNUM_COV,DK_BETA,readmax,covmax)

order = str2num(order);
VT= str2num(VT);
maxchr= str2num(maxchr);
kmax= str2num(kmax);
K2= str2num(K2);
rc= str2num(rc);
RESOL= str2num(RESOL);
km_ver= str2num(km_ver);
SGE= str2num(SGE);
win_resol= str2num(win_resol);
DNUM_COV= str2num(DNUM_COV);
DK_BETA= str2num(DK_BETA);
readmax= str2num(readmax);
covmax= str2num(covmax);
%%%%%%%%%%%%%%%%%%%%%
% Hardcoding params %
%%%%%%%%%%%%%%%%%%%%%

valid_limit = 2147483647;
blfile_suffix = '.bin'; 
snp_suffix = '.vcf';
seq_suffix = '.fa';
combine_prefix = 'allChr';

CountBase_cmd = 1;
GetBaseline_cmd = 2;
GetKmermat_cmd = 3;
KSS_cmd = 4;
combine_KSS_cmd = 5;
combine_win_KSS_cmd = 6;
KSS_STR_cmd = 7;

%%% Params depending on input params
K=K2/2;
baseline_dir = horzcat(topdir,'/baseline/');
kss_raw_dir = horzcat(topdir,'/snp-ranking-raw-result//');
kss_raw_SGE_dir = horzcat(kss_raw_dir,'SGE_log/');
kmer_mat_file = horzcat(topdir,'/',kmer_binary_folder,'.kmer.mat');
x0file = horzcat(topdir,'/',kmer_binary_folder,'.bin');

if (km_ver==0)
	param1 = 'xopt';
	param2 = 'copt';
else
	param1 = 'yall'; 
	param2 = 'x0';
	if (km_ver>=2)
		kmer_binary_folder_in = horzcat([kmer_binary_folder,'/input/']);
		kmer_binary_folder = horzcat([kmer_binary_folder,'/output/']);
	end
end
if order==KSS_STR_cmd
	snp_suffix = '_estr.tab';
end

%%% Checking directory
if (exist(baseline_dir)~=7)
	mkdir(baseline_dir);
end
if (exist(kss_raw_dir)~=7)
	mkdir(kss_raw_dir);
end
if (SGE==1 && exist(kss_raw_SGE_dir)~=7)
	mkdir(kss_raw_SGE_dir);
end

%%% Main functions

if (~isempty(find(order==CountBase_cmd)))
	display('Start Counting Base');
	cmd = horzcat('python preprocess/seq_data/countBase.py ', num2str(maxchr),' ',seq_dir);
    system(cmd);
end


if (~isempty(find(order==GetBaseline_cmd)))
	if ~isdeployed
		addpath('preprocess/seq_data');
	end
	[chrchunk,genomesize] = findValidationCutoff(maxchr,seq_dir,valid_limit);

	[nrow,~] = size(chrchunk);
	for i=1:nrow
		display(horzcat('Getting Basline for chr',num2str(chrchunk(i,1)),'-',num2str(chrchunk(i,2))));
		out = fopen('preprocess/baseline/get_background_params.R','w');
		fprintf(out,'%s\n',horzcat('kmax=',num2str(kmax)));
		fprintf(out,'%s\n',horzcat('K=',num2str(K)));
		fprintf(out,'%s\n',horzcat('K2=',num2str(K2)));
		fprintf(out,'%s\n',horzcat('genomestart=',num2str(genomesize(i,1))));
		fprintf(out,'%s\n',horzcat('genomeend=',num2str(genomesize(i,2))));
		fprintf(out,'%s\n',horzcat('RESOLUTION=',num2str(RESOL)));
		blfileName = horzcat(blfile_prefix,'.',num2str(chrchunk(i,1)),'to',num2str(chrchunk(i,2)),blfile_suffix);
		fprintf(out,'%s\n',horzcat('bg.file="',baseline_dir,blfileName,'"'));
		fprintf(out,'%s\n',horzcat('genomefile="',genomefile,'"'));
		fprintf(out,'%s\n',horzcat('OUT_DIR="',kmer_binary_dir,'"'));
		fprintf(out,'%s\n',horzcat('ename="',kmer_binary_folder,'"'));
		fprintf(out,'%s\n',horzcat('ename_in="',kmer_binary_folder_in,'"'));
		fprintf(out,'%s\n',horzcat('param1="',param1,'"'));
		fprintf(out,'%s\n',horzcat('param2="',param2,'"'));
		fprintf(out,'%s\n',horzcat('x0filename="',x0file,'"'));
		fprintf(out,'%s\n',horzcat('DLINK="',DLINK,'"'));
		fprintf(out,'%s\n',horzcat('DNUM_COV=',num2str(DNUM_COV)));
		fprintf(out,'%s\n',horzcat('DK_BETA=',num2str(DK_BETA)));
		fprintf(out,'%s\n',horzcat('readmax=',num2str(readmax)));
		fprintf(out,'%s\n',horzcat('covmax=',num2str(covmax)));
		
        fclose(out);
		%cd preprocess/baseline/;
        if (km_ver==3)
            r= system(horzcat('Rscript preprocess/baseline/get_background_hg19_covar.R /scripts/'));
        else
		    r= system(horzcat('Rscript preprocess/baseline/get_background_hg19.R /scripts/'));
        end
		if ~isdeployed
			addpath('preprocess/baseline');
		end

		display('start parse baseline');
		parseBaseline(chrchunk(i,1),chrchunk(i,2),baseline_dir,seq_dir,blfileName,blfile_outfile_suffix,rc);
	
		%cd('../../');
	end
end

if (~isempty(find(order==GetKmermat_cmd)))
	display('Getting Kmermat');
	out = fopen('preprocess/kmer_mat/prepare_kmer_mat_params.R','w');
	fprintf(out,'%s\n',horzcat('OUT_DIR="',kmer_binary_dir,'"'));
	fprintf(out,'%s\n',horzcat('ename="',kmer_binary_folder,'"'));
	fprintf(out,'%s\n',horzcat('K2=',num2str(K2)));
	fprintf(out,'%s\n',horzcat('RESOL=',num2str(RESOL)));
	fprintf(out,'%s\n',horzcat('kmax=',num2str(kmax)));
	fprintf(out,'%s\n',horzcat('outfile="',kmer_mat_file,'"'));
	fprintf(out,'%s\n',horzcat('param1="',param1,'"'));
	fprintf(out,'%s\n',horzcat('param2="',param2,'"'));
	
	fclose(out);
	%cd preprocess/kmer_mat/;
	display('Start to run R file')
    r= system(horzcat('Rscript preprocess/kmer_mat/prepare_kmer_mat.R /scripts/'));
    %cd ../../
	if ~isdeployed
		addpath('preprocess/kmer_mat');
	end
	display('Start to correct for resol')
	correct_for_resol(kmer_mat_file,RESOL);
end


if (~isempty(find(order==KSS_cmd)))
	display('Start Kmer-SNP-Scoring');
	main_baseline(num2str(maxchr),baseline_dir,snp_dir,seq_dir,kss_raw_dir,kmax,K,kmer_mat_file,blfile_outfile_suffix,snp_suffix,seq_suffix,kss_suffix,VT,VT_prefix,VT_suffix,genomefile,0,-1,win_resol);
	
	%out = fopen('snp-kmer-scoring/parellel_run_params.R','w');
	%fprintf(out,'%s\n',horzcat('kmermat_file="',kmer_mat_file,'"'));
	%fprintf(out,'%s\n',horzcat('snpdir="',snp_dir,'"'));
	%fprintf(out,'%s\n',horzcat('seqdir="',seq_dir,'"'));
	%fprintf(out,'%s\n',horzcat('bldir="',baseline_dir,'"'));
	%fprintf(out,'%s\n',horzcat('outdir="',kss_raw_dir,'"'));
	%fprintf(out,'%s\n',horzcat('VT=',num2str(VT)));
	%fprintf(out,'%s\n',horzcat('VT_prefix="',VT_prefix,'"'));
	%fprintf(out,'%s\n',horzcat('VT_suffix="',VT_suffix,'"'));
	%fprintf(out,'%s\n',horzcat('windowsize=',num2str(K)));
	%fprintf(out,'%s\n',horzcat('ksize=',num2str(kmax)));
	%fprintf(out,'%s\n',horzcat('SGE=',num2str(SGE)));
	%fprintf(out,'%s\n',horzcat('SGE_dir="',kss_raw_SGE_dir,'"'));
	%fprintf(out,'%s\n',horzcat('maxchr=',num2str(maxchr)));
	%fprintf(out,'%s\n',horzcat('bl_suffix="',blfile_outfile_suffix,'"'));
	%fprintf(out,'%s\n',horzcat('snp_suffix="',snp_suffix,'"'));
	%fprintf(out,'%s\n',horzcat('seq_suffix="',seq_suffix,'"'));
	%fprintf(out,'%s\n',horzcat('out_suffix="',kss_suffix,'"'));
	%fprintf(out,'%s\n',horzcat('genomefile="',genomefile,'"'));
	%fprintf(out,'%s\n',horzcat('TOPKMER=',num2str(0)));
	%fprintf(out,'%s\n',horzcat('topkmer_thresh=',num2str(-1)));
	%fprintf(out,'%s\n',horzcat('win_resol=',num2str(win_resol)));
	%fclose(out);
	
	%cd snp-kmer-scoring/;
	%r= system('Rscript snp-kmer-scoring/parellel_run.R /scripts/');
	%cd ('../');
end

if (~isempty(find(order==combine_KSS_cmd)))
	if ~isdeployed
		addpath('postprocess/combineKSS/');
	end
	combine_snp(kss_raw_dir,kss_suffix,maxchr,combine_prefix)
end

if (~isempty(find(order==combine_win_KSS_cmd)))
	addpath('postprocess/combineKSS/');
	combine_snp_win(kss_raw_dir,kss_suffix,maxchr,combine_prefix)
end

if (~isempty(find(order==KSS_STR_cmd)))
	display('Start STR-HIT-Scoring');
    str_main_baseline_brutal(num2str(maxchr),baseline_dir,snp_dir,seq_dir,kss_raw_dir,kmax,K,kmer_mat_file,blfile_outfile_suffix,'_estr.tab',seq_suffix,kss_suffix,genomefile,x0file)	
	%out = fopen('snp-kmer-scoring/str_parellel_run_params.R','w');
	%fprintf(out,'%s\n',horzcat('kmermat_file="',kmer_mat_file,'"'));
	%fprintf(out,'%s\n',horzcat('snpdir="',snp_dir,'"'));
	%fprintf(out,'%s\n',horzcat('seqdir="',seq_dir,'"'));
	%fprintf(out,'%s\n',horzcat('bldir="',baseline_dir,'"'));
	%fprintf(out,'%s\n',horzcat('outdir="',kss_raw_dir,'"'));
	%fprintf(out,'%s\n',horzcat('windowsize=',num2str(K)));
	%fprintf(out,'%s\n',horzcat('ksize=',num2str(kmax)));
	%fprintf(out,'%s\n',horzcat('SGE=',num2str(SGE)));
	%fprintf(out,'%s\n',horzcat('SGE_dir="',kss_raw_SGE_dir,'"'));
	%fprintf(out,'%s\n',horzcat('maxchr=',num2str(maxchr)));
	%fprintf(out,'%s\n',horzcat('bl_suffix="',blfile_outfile_suffix,'"'));
	%fprintf(out,'%s\n',horzcat('snp_suffix="_estr.tab','"'));
	%fprintf(out,'%s\n',horzcat('seq_suffix="',seq_suffix,'"'));
	%fprintf(out,'%s\n',horzcat('out_suffix="',kss_suffix,'"'));
	%fprintf(out,'%s\n',horzcat('genomefile="',genomefile,'"'));
	%fprintf(out,'%s\n',horzcat('x0file="',x0file,'"'));
	%fclose(out);
	%cd snp-kmer-scoring/;
	%r= system('Rscript str_parellel_run.R');
	%cd ('../');
end


