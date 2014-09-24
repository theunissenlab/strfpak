function CS = small_autocorr3(stimval,nband,size_CS,twin)
%  Same as small_autocorr, but without the ntrials multiplication, which
%  can be done outside.
xb = 1;
CS =  zeros(nband/2*(nband+1),2*twin+1);% zeros(size_CS);

N = size(stimval,2);
savedStimFft = fft(stimval',2^nextpow2(2*N-1))';

for ib1 = 1:nband

    if 0==1%running_flag == 1
      waitbar(ib1/nband, tempWait);
    end

    for ib2 = ib1:nband
        c = real(ifft(conj(savedStimFft(ib2,:)).*(savedStimFft(ib1, :))));
        

        % NEW version of algorithm by using xcorr
        %size(c(end-twin+1:end))
        %size(c(1:twin+1))
        CS(xb, :) = [c(end-twin+1:end) c(1:twin+1)];

        xb = xb +1;
    end            % end of ib2
end                 % end of ib1 
