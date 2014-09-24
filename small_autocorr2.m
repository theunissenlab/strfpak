function CS = small_autocorr2(stimval,nband,size_CS,twin)
%  Same as small_autocorr, but without the ntrials multiplication, which
%  can be done outside.
xb = 1;
CS =   zeros(size_CS);
for ib1 = 1:nband

    if 0==1%running_flag == 1
      waitbar(ib1/nband, tempWait);
    end

    for ib2 = ib1:nband

        % NEW version of algorithm by using xcorr
        CS(xb, :) = CS(xb, :) + xcorr(stimval(ib2,:), stimval(ib1, :),...
            twin);

        xb = xb +1;
    end            % end of ib2
end                 % end of ib1 
