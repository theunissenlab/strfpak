function [JN,y,stP] = make_slepian(x,WinLength,NW,nTapers,nChannels,nFFTChunks,nFFT,Detrend,winstep)
% allocate memory now to avoid nasty surprises later
stP = zeros(nFFT, nChannels, nChannels);
varP = zeros(nFFT, nChannels, nChannels);
[Tapers V]=dpss(WinLength,NW,nTapers, 'calc');
Periodogram = complex(zeros(nFFT, nTapers, nChannels)); % intermediate FFTs
Temp1 = complex(zeros(nFFT, nTapers)); %Temps are particular psd or csd values for a frequency and taper
Temp2 = complex(zeros(nFFT, nTapers));
Temp3 = complex(zeros(nFFT, nTapers));
eJ = complex(zeros(nFFT,1));
JN = complex(zeros(nFFTChunks,nFFT, nChannels, nChannels));  %jackknifed csd
y=complex(zeros(nFFT, nChannels, nChannels)); % output array for csd
Py=zeros(nFFT, nChannels, nChannels); % output array for psd's


% New super duper vectorized alogirthm
% compute tapered periodogram with FFT 
% This involves lots of wrangling with multidimensional arrays.

TaperingArray = repmat(Tapers, [1 1 nChannels]);
for j=1:nFFTChunks
	Segment = x((j-1)*winstep+[1:WinLength], :);
	if (~isempty(Detrend))
		Segment = detrend(Segment, Detrend);
	end;
	SegmentsArray = permute(repmat(Segment, [1 1 nTapers]), [1 3 2]);
	TaperedSegments = TaperingArray .* SegmentsArray;

	Periodogram(:,:,:) = fft(TaperedSegments,nFFT);

	% Now make cross-products of them to fill cross-spectrum matrix
	for Ch1 = 1:nChannels
		for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
			Temp1 = squeeze(Periodogram(:,:,Ch1));
			Temp2 = squeeze(Periodogram(:,:,Ch2));	
			Temp2 = conj(Temp2);
			Temp3 = Temp1 .* Temp2;

            %eJ and eJ2 are the sum over all the tapers.
            eJ=sum(Temp3, 2)/nTapers;
            JN(j,:,Ch1, Ch2) = eJ;
			y(:,Ch1, Ch2)= y(:,Ch1,Ch2) + eJ;
		end
	end
end 

% now fill other half of matrix with complex conjugate
for Ch1 = 1:nChannels
	for Ch2 = (Ch1+1):nChannels % don't compute cross-spectra twice
		y(:, Ch2, Ch1) = y(:,Ch1,Ch2);
        Py(:, Ch1, Ch2) = atanh(abs(y(:,Ch1,Ch2)./sqrt(abs(y(:,Ch1,Ch1)).*abs(y(:,Ch2,Ch2)))));
	end
end


for j = 1:nFFTChunks
    JN(j,:, :, :) = abs(y - squeeze(JN(j,:, :,:)));
    for Ch1 = 1:nChannels
        for Ch2 = (Ch1+1):nChannels  
            % Calculate the transformed coherence
            JN(j,:, Ch1, Ch2) =atanh(abs(JN(j,:,Ch1,Ch2))./sqrt(abs(JN(j,:,Ch1,Ch1)).*abs(JN(j,:,Ch2,Ch2))));
            % Obtain the pseudo values
            JN(j,:, Ch1, Ch2) = nFFTChunks*Py(:, Ch1, Ch2)' - (nFFTChunks-1)*squeeze(JN(j,:, Ch1,Ch2));
        end
    end  
end

meanP=squeeze(mean(JN,1));
for Ch1=1:nChannels
    for Ch2=Ch1:nChannels
        varP(:,Ch1, Ch2) = (1/nFFTChunks)*var(JN(:,:,Ch1, Ch2),1);
    end
end

%upper and lower bounds will be 2 standard deviations away.
stP=sqrt(varP);
