function [] = combine_snp(dir,snpfile_suffix,MaxChr,combine_prefix)

outfile = horzcat([dir,combine_prefix,snpfile_suffix,'.mat']);
sizefile =  horzcat([dir,combine_prefix,snpfile_suffix,'.size.mat']);
chr = cell(1,MaxChr);
snploc = cell(1,MaxChr);
score = cell(1,MaxChr);
snpnrow = zeros(1,MaxChr);
numl = 0;
for i = 1:MaxChr
    display(horzcat('combining chr',num2str(i)));
	snpfile = horzcat([dir,'chr',num2str(i),snpfile_suffix]);
    if exist(snpfile)~=2
		continue
	end
    [~,snploc{i},~,~,score{i},~] = textread(snpfile,'%s %f %s %s %f %s');
    [snpnrow(i),~] = size(snploc{i});
    chr{i} = ones(snpnrow(i),1).*i;
    numl = numl + snpnrow(i);   
end

mat = zeros(numl,3);
idx = 1;
for i = 1:MaxChr
    mat(idx:idx+snpnrow(i)-1,:) = [chr{i},snploc{i},score{i}];
    idx = idx + snpnrow(i);
end

meaningful = mat(:,1)~=-1;
mat = mat(meaningful,:);
save(outfile,'mat');
save(sizefile,'snpnrow');


