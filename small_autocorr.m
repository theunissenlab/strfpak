function CS = small_autocorr(stimval,nband,ntrials,size_CS,twin)
xb = 1;
CS = zeros(size_CS);
for ib1 = 1:nband

    if 0%running_flag == 1
      waitbar(ib1/nband, tempWait);
    end

    for ib2 = ib1:nband

        % NEW version of algorithm by using xcorr
        CS(xb, :) = CS(xb, :) + xcorr(stimval(ib2,:), stimval(ib1, :),...
            twin)*ntrials;

        xb = xb +1;
    end            % END of ib2
end                 % END of ib1
