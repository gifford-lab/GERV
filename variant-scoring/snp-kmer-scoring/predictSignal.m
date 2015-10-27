function [dist] = predictSignal(seq,effectrange,kmermat,ksize,K,mi,x0)
dist = zeros(1,length(seq)+2*K-1);
seq = seq + 1;
seq(find(seq==5),:) = 1;
for i=1:length(seq)
	minj = min(length(seq),i+ksize-1);
	knum_count = minj-i+1;
	kmernum = zeros(1,knum_count);
	for j=i:minj
		kmernum(j-i+1) = mi(1:(j-i+1)) *  seq(i:j);
	end
	%display(kmernum);
	%display(seq(i:minj));
	kmermat_part = kmermat(kmernum,:);
	
	eff_start = i;
	eff_end = 2*K +i - 1;
	dist(1,eff_start:eff_end) = dist(1,eff_start:eff_end) + ones(1,knum_count)*kmermat_part;
end
goodpick_start = K+effectrange(1);
goodpick_end = K + effectrange(2);
dist = exp(dist(goodpick_start:goodpick_end)-x0(1));

