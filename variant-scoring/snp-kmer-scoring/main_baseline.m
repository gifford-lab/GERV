function [] = main_baseline(chr,bldir,snpdir,seqdir,outdir,ksize,K,kmermat_file,bl_suffix,snp_suffix,seq_suffix,out_suffix,VT,VT_prefix,VT_suffix,genomebin,TOPKMER,topkmer_thresh,resol)
diaryfile = fullfile(outdir,horzcat('scoring_log_chr',chr));
if (exist(diaryfile)==2)
	system(horzcat('rm ',diaryfile));
end
diary(diaryfile);

load(kmermat_file);
mi = ones(1,ksize);
for i=2:ksize
	mi(i) = mi(i-1) * 4;
end

K2 = K*2;
% windowsize = 200;
% ksize = 8;
allsizeFile = fullfile(seqdir,'all.size.txt');
allsize = textread(allsizeFile);

sizeFile = fullfile(seqdir,horzcat(['chr',chr,'.size.txt']));
gsize = textread(sizeFile);

genomeread = fopen(genomebin,'r');
fseek(genomeread,allsize(str2num(chr)),-1);
genome = fread(genomeread,gsize);

fclose(genomeread);

blfile = fullfile(bldir,horzcat(['chr',chr,bl_suffix]));
bl = fopen(blfile);
baseline = fread(bl,'*single');
if (length(baseline)~=gsize)
    error('baseline size mismatches with size file!');
end
fclose(bl);

snpFile = fullfile(snpdir,horzcat(['chr',chr,snp_suffix]));
snp = fopen(snpFile);

outputfile =  fullfile(outdir,horzcat(['chr',chr,out_suffix]));
out = fopen(outputfile, 'w');

distr_range = 2*K +ksize-1;

keys = {'A','T','G','C','a','t','g','c'};
values = {0,1,2,3,0,1,2,3};
dic = containers.Map(keys, values);

keys = {0,1,2,3,4};
values = {'A','T','G','C','N'};
revdic = containers.Map(keys, values);

kmermap = zeros((1+ksize)*ksize/2,2);
cnt = 0;
for i = 1:ksize
	for j = ksize:(i+ksize-1)
		cnt = cnt + 1;
		kmermap(cnt,:) = [i,j];
	end
end


oneLine = fgetl(snp);
while ischar(oneLine)
    myVector = regexp(oneLine, '\t','split');        %%% if it is not VCF, this may change
    first = myVector{1};
	
	%% Filter for header    
    if (first(1) == '#')
        oneLine = fgetl(snp);
        continue;
    end
    
	%% Parse SNP information and Filter for non-SNP
    if (VT==1)
        info = myVector{8};
        vt_split = regexp(info,VT_prefix,'split');
        if (length(vt_split{2})<3)
            oneLine = fgetl(snp);
            continue;
        end
        variant_type = vt_split{2}(1:3);
        if (strcmp(variant_type,VT_suffix)~=1)
            oneLine = fgetl(snp);
            continue;
        end

    end
    
    chr = first;
    pos = str2double(myVector{2});
    ref = myVector{4};
    if (length(ref)>1)
        oneLine = fgetl(snp);
        continue;
    end
        
    alt = myVector{5};
    alt_s = regexp(alt, '[,]','split');
    flag = 0;
    for alt_iter = 1:length(alt_s)     
        if (length(alt_s{alt_iter})>1)
            flag = 1;
            break;
        end
    end
    if (flag==1)
        oneLine = fgetl(snp);
        continue;
    end
    alt_num = length(alt_s);
	
	%%% Find the sequence around the SNP
	genomeleft = max(1,pos - ksize + 1);
	genomeright = min(gsize,pos + ksize - 1);
	snppos = min(pos,ksize);

	seq = genome(genomeleft:genomeright);
    
	%%% Find baseline 
	blleft = max(1,pos-K-(ksize-1));
	leftcut = blleft - (pos-K-(ksize-1)) + 1;
	blright = min(pos+K-1,gsize);
	rightcut = distr_range - (pos+K-1 - blright);

	base = baseline(blleft:blright);

    %sig = sum(base);
    %all_sig_str = num2str(sig);
    %fprintf(out, '%s\t%d\t%s\t%s\t%f\t%s\n', chr,pos,ref,alt,sig,all_sig_str);
    %oneLine = fgetl(snp);
    %continue;

	%%% Filter for too many Ns
	ns = find(seq==4);
    n = length(ns);

	if (n>3)
		if (TOPKMER)
			error('TOPKMER while having too many Ns in seq');
		end
		fprintf(out, '%s\t%d\t%s\t%s\t%f\t%s\n', chr,pos,ref,alt,-1,'-1');
		oneLine = fgetl(snp);
		continue;
	end
	
	%%% Aggregate the influence vector of all kmers overlapping the SNP

	if (TOPKMER)
		distr = zeros(1+alt_num,distr_range,(1+ksize)*ksize/2);
	else
		distr = zeros(1+alt_num,distr_range);
		%distr2 = zeros(1+alt_num,1);
		kmersupport = zeros(1+alt_num,(1+ksize)*ksize/2);
	end

	if (n>0)
		num = 4^n;
    	nmer = generateKmer(n);
	else
		num = 1;
	end

	for multi_var = 1:num
		seq_to_score = seq;
	    
		%%% If there are several 'N' in the sequence, we need to enumerate all possible values of them and take average
		for multi_pos = 1:n
			seq_to_score(ns(multi_pos)) =nmer(multi_var,multi_pos);
		end

		%%%% Check if ref allele matches the base in ref genome at the SNP location
		if (dic(ref)~=seq_to_score(snppos))
			display(chr);
			display(pos);
			display(ref);
			display(seq_to_score(snppos));
			display(seq(snppos));
			error('ref doesn''t match!')
		end

		%%%% Start aggregating
		if (TOPKMER)
			%%%% When SNP = ref
			distr(1,:,:) = distr(1,:,:) + aggregateKmer_topkmer(seq_to_score,snppos,K2,ksize,mat,distr_range);
			
			%%%% When SNP = alt
			for i = 1:alt_num
				seq_to_score(snppos) = dic(alt_s{i});
				distr(i+1,:,:) = distr(i+1,:,:) + aggregateKmer_topkmer(seq_to_score,snppos,K2,ksize,mat,distr_range);
			end
		else
			%%%% When SNP = ref
			%distr2(1,:) = distr2(1,:) +aggregateKmer_logspace(seq_to_score,snppos,K2,ksize,mat,distr_range,mi);
			%distr(1,:) = distr(1,:) +aggregateKmer(seq_to_score,snppos,K2,ksize,mat,distr_range,mi);
			[a,b] = aggregateKmer(seq_to_score,snppos,K2,ksize,mat,distr_range,mi);
			distr(1,:) = distr(1,:) + a;
			kmersupport(1,:) = b;
			%%%% When SNP = alt
			for i = 1:alt_num
				seq_to_score(snppos) = dic(alt_s{i});
				%distr2(i+1,:) = distr2(i+1,:) + aggregateKmer_logspace(seq_to_score,snppos,K2,ksize,mat,distr_range,mi);
				%distr(i+1,:) = distr(i+1,:) + aggregateKmer(seq_to_score,snppos,K2,ksize,mat,distr_range,mi);
				[a,b] = aggregateKmer(seq_to_score,snppos,K2,ksize,mat,distr_range,mi);
				distr(i+1,:) = distr(i+1,:) + a;
				kmersupport(i+1,:) = b;
			end
		end
	end
	distr = distr ./ num;
    
	%%% Truancate distr if nessacary
	if (TOPKMER)
    	distr = distr(:,leftcut:rightcut,:);
	else
    	distr = distr(:,leftcut:rightcut);	
	end
	
	%%% Score the SNP based on changes in aggreagted influence vector
	if (TOPKMER)
    	distr_score = zeros(1+alt_num,distr_range);
		kmer_num = (1+ksize)*ksize/2;
		part = zeros(distr_range,kmer_num);
		for i = 1 : (1+alt_num)
			part(:,:) = distr(i,:,:);
			distr_score(i,:) = part * ones(kmer_num,1);
		end
   		%%% Find the SNP score as reference 
		[sig,all_sig_str ] = exp_errorL2_BL(base,distr_score',alt_num);

		%%% Look for kmers contribute more than threshold to SNP score
		[topkmer_score,topkmer_seq,topkmer_choose] = exp_errorL2_BL_topkmer(base,distr,alt_num,snppos,seq,ksize,topkmer_thresh,sig,revdic);
		
		%%% Output
		for i = 1:length(topkmer_score)
			fprintf(out, '%s\t%d\t%s\t%s\t%f\t%f\t%f\t%s\t%s\n', chr,pos,ref,alt,sig,topkmer_score(i),abs(topkmer_score(i)/sig),alt_s{topkmer_choose(i)},topkmer_seq{i});
		end
	else
		%%% Score the SNP based on the change between ref and alt
		[sig,all_sig_str ] = exp_errorL2_BL(base,distr',alt_num);
		%[sig2,all_sig_str2 ] = logspace_diff(distr2,alt_num);
		%all_sig_str = num2str(sig2);
		refseq = seq;
		altseq = seq;
		altseq(snppos) = dic(alt_s{1});
		all_sig_str = getKeyKmer(kmersupport,refseq,altseq,ksize,revdic,kmermap);
		%%% Output
		fprintf(out, '%s\t%d\t%s\t%s\t%f\t%s\n', chr,pos,ref,alt,sig,all_sig_str);
	end

    oneLine = fgetl(snp);
end

fclose(out);
fclose(snp);
