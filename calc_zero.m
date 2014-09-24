function [stim_stat_zero] = calc_zero(dwt, dwf, fwidth, stim_stat);
% Zero outs the modulation spectrum outside the rectangle.
% k is the constant that relates f_width to the upper frequency of the
% sampled rectangle
k = 2.57;
wt_max = k*fwidth;

% tp1 and tp2 are the index points for dwt boundaries
nwt = length(dwt);
tp1 = 1;
tp2 = nwt;
for i=1:nwt-1
    if (dwt(i) <= -wt_max & dwt(i+1) > -wt_max )
        tp1 = i;
    end
    if (dwt(i) < wt_max & dwt(i+1)>= wt_max )
        tp2 = i+1;
    end
end

% fp1 and fp2 are the index points for the dwf boundaries
wf_max = 1./wt_max;
nwf = length(dwf);
fp1 = 1;
fp2 = nwf;
for i=1:nwf-1
    if (dwf(i) <= -wf_max & dwf(i+1) > -wf_max )
        fp1 = i;
    end
    if (dwf(i) < wf_max & dwf(i+1)>= wf_max )
        fp2 = i+1;
    end
end

% Zero out areas outside rectangle
stim_stat_zero = stim_stat;
stim_stat_zero(:,1:tp1) = 0;
stim_stat_zero(:,tp2:nwt) = 0;
stim_stat_zero(1:fp1,:) = 0;
stim_stat_zero(fp2:nwf,:) = 0;
