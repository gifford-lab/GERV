function [] = correct_for_resol(file,resol)
load(file);
[nrow,ncol] = size(mat);
out = zeros(nrow,ncol*resol);
pick = 1:ncol;
pick = pick*resol - (resol-1);

for i = 1:resol
	out(:,pick) = mat;
	pick = pick + 1;
end
mat = out;
system(horzcat('rm ',file));
save(file,'mat','-v7.3');

