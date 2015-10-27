function [sig,all_sig_str] = exp_errorL2_BL(base,distr,alt_num)
[distr_range,~] = size(distr);
base_expand = log(base)*ones(1,alt_num);
gain =  - distr(:,1)*ones(1,alt_num) + distr(:,2:(alt_num+1));
new = gain + base_expand;
diff = -exp(base_expand) + exp(new);
all_sig = sqrt(sum(diff.^2,1));
all_sig_str = num2str(all_sig(1));
if (length(all_sig)>=2)
    for i = 2:length(all_sig)
        all_sig_str = horzcat([all_sig_str,',',num2str(all_sig(i))]);
    end
end
[~,pick ] =max(abs(all_sig));
sig = all_sig(pick);
end
