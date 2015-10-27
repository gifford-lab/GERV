function [chrchunk,genomesize] = findValidationCutoff(maxchr,dir,limit)
accum = 0;
cutcnt = 0;
chrchunk = zeros(maxchr,2);
genomesize = zeros(maxchr,2);
all = 0;
for i = 1:maxchr
	file = horzcat(dir,'chr',num2str(i),'.size.txt');
	if (exist(file)==0)
		error('Fail to open chr size file')
	end

	cnt = dlmread(file);
	if (accum+cnt>limit)
		cutcnt = cutcnt +1;
		chrchunk(cutcnt,2) = i-1;
		chrchunk(cutcnt+1,1) = i;
		genomesize(cutcnt,2) = all;
		genomesize(cutcnt+1,1) = all;
		accum = cnt;
	else
		accum = accum+cnt;
	end
	all = all + cnt;

end
chrchunk(1,1) = 1;
chrchunk(cutcnt+1,2) = maxchr;
genomesize(1,1) = 0;
genomesize(cutcnt+1,2) = all;

chrchunk = chrchunk(1:cutcnt+1,:);
genomesize = genomesize(1:cutcnt+1,:);
