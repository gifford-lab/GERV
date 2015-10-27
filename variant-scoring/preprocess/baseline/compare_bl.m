close all;clear;clc;
olddir = '/cluster/zeng/research/kmer/hg19/K562/baseline/';
newdir = '/cluster/zeng/research/kmer/hg19/TEST_K562/baseline/';

oldfile = horzcat([olddir,'fitbg.hg19.bin']);
newfile = horzcat([newdir,'fitbg.K562.1to12.bin']);

old = fopen(oldfile,'r');
new = fopen(newfile,'r');
oldbl = fread(old,'*single');
newbl = fread(new,'*single');

old_padding = 0;
new_padding = old_padding;

[nrow,~] = size(newbl);

flag = -1;
tolerance = 1e-3;
badcnt = 0;
maxdiff = 0;
for i = 1:nrow
    %if (abs(oldbl(i+old_padding)-newbl(i+new_padding))>tolerance)
    %    flag = 1;
	%	pdisplay(horzcat('diff at bit ',num2str(i)));
    %    badcnt = badcnt + 1;
	%	thisdiff = abs(oldbl(i+old_padding)-newbl(i+new_padding));
	%	if thisdiff>diff
	%		diff = thisdiff;
	%	end
    %end
end
if (flag==-1)
	display('Match completely!');
end


diff = oldbl(1:nrow)-newbl;
sdiff = sort(diff);
binsize = 0.0001;
MIN = sdiff(1);
MAX = sdiff(nrow);
bin = MIN+binsize/2;
idx = 1;
x = [];
y = [];
while(bin-binsize/2<=MAX)
	newidx=idx;
	while(sdiff(newidx)<bin+binsize/2)
		newidx = newidx + 1;
		if (newidx>nrow)
			break;
		end
	end
	x	= [x,bin];
	y = [y,newidx - idx];
	bin = bin + binsize;
	idx = newidx;
end
h=figure;
bar(x,log(y+ones(1,length(y))));
xlabel('Diff');
ylabel('Log of (Cnt+1)');
saveas(h,'/cluster/zeng/research/kmer/hg19/TEST_K562/compare_result_1to12','png');

h=figure;
bar(x,y/nrow*100);
xlabel('Diff');
ylabel('Percentage');
saveas(h,'/cluster/zeng/research/kmer/hg19/TEST_K562/compare_result_1to12_pct','png');
