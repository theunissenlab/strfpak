function est_spike_k = predict_one_tol_std_combo_2(forwardJN_s,nlen,end_window,fidx,nband,stim_avg)
%  This routine calculates the predicted PSTH in a way that's easy to cache.
global predDS

stim_env = Check_And_Load(predDS{fidx}.stimfiles);
stimval = zeros(nband, nlen);
for tgood = 1:nlen
    stimval(:, tgood) = stim_env(1:nband, tgood) - stim_avg(1:nband);
end

% NEW algorithm
est_spike_k = zeros(nlen+end_window,1);
clear tempResult


for ib = 1:nband

    tempResult = conv(stimval(ib, :), forwardJN_s(ib, :));
    chopsize = size(forwardJN_s, 2);

    % Chop to our wanted length
    chopsize = (chopsize +1) /2;
    est_spike_k = est_spike_k +...
        tempResult(1,chopsize:size(tempResult, 2) - chopsize +1)';
end     % END of ib

