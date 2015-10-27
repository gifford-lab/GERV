function [kmer]= generateKmer(k)
kmer = zeros(4^k,k);
for i = 0:(4^k-1)
    str = '';
    cur = i;
    for j = k:-1:1
        num = cur/(4^(j-1));
		num = floor(num);
        kmer((i+1),j) = num;        
		cur = cur-(num)*(4^(j-1));
    end
end
end
