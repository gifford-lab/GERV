function [out] = getKeyKmer(support,ref,alt,ksize,revdic,kmermap)
topnum = 5;
diff = abs(support(1,:) - support(2,:));
[~,b] = sort(diff,'descend');
b = b(1:topnum);

out = '';
for num = 1:topnum
	spos = kmermap(b(num),1);
	epos = kmermap(b(num),2);
	for k = spos :epos
		out = horzcat(out,revdic(ref(k)));
	end
	out = horzcat(out,',');
	for k = spos:epos
		out = horzcat(out,revdic(alt(k)));
	end
	out = horzcat(out,';');
end
