function est_spike_k = predict_one_tol_std_combo(forward,forwardJNstd,filecount,fidx,stdfilt,nlen,end_window,nband,stimval,matchflg);
%  This routine calculates the predicted PSTH in a way that's easy to cache.
if matchflg == 1
    % Compute the filtered JN STRF
    forwardJN_s(:,:) = fast_filter_filter(mean(forward(:,:,[1:(fidx-1) (fidx+1):filecount]),3), forwardJNstd(:,:,fidx), stdfilt);
else
    forwardJN_s(:,:) = fast_filter_filter(forward, squeeze(mean(forwardJNstd,3)), stdfilt);
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

