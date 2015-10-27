function [] = parseBaseline(startchr,endchr,dir,sizedir,blfileName,outfile_suffix,rc)

blfile = horzcat([dir,blfileName]);
bl = fopen(blfile,'r');

offset = 0;
for i=startchr:endchr
% if (i==Maxchr)
%     sizefile = horzcat([dir,'chr','X','.size.txt']);
%     outfile = horzcat([dir,'chr','X','.baseline.bin']);
% else
    sizefile = horzcat([sizedir,'chr',num2str(i),'.size.txt']);
    outfile = horzcat([dir,'chr',num2str(i),outfile_suffix]);
% end
s = textread(sizefile);
out = fopen(outfile,'w');

baseline = fread(bl,[s,1],'*single');
fwrite(out,baseline,'*single');
if (rc==1)
	offset = offset + 2*s;
	fread(bl,[s,1],'*single');
else
	offset = offset + s;
end

fclose(out);
end
left = fread(bl,'*single');
display(horzcat(['Bytes in tail: ',num2str(length(left)*4)]));
display(horzcat(['Bytes in data: ',num2str(offset*4)]));
fclose(bl);
