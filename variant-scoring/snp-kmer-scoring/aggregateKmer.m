function [dist,kmersupport] = aggregateKmer(seq,snppos,K2,ksize,kmermat,range,mi)
dist = zeros(1,range);
seq = seq + 1;
kmersupport = [];
for i=1:snppos
	minj = min(length(seq),i+ksize-1);
	kmernum = zeros(1,minj-snppos+1);
	for j=snppos:minj
		kmernum(j-snppos+1) = mi(1:(j-i+1)) *  seq(i:j);
	end
	kmermat_part = kmermat(kmernum,:);
	dist(1,i:(i+K2-1)) = dist(1,i:(i+K2-1)) + ones(1,minj-snppos+1) * kmermat_part;
	kmersupport = [kmersupport;kmermat_part * ones(K2,1)];
end
kmersupport = kmersupport';
