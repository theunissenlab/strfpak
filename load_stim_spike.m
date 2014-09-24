% Get dimensions for the time, f-bands and JN 
nstd_val = 0.5;
name = pwd;
load stim_init_count.avg;
nJN = stim_init_count(1);
nb = stim_init_count(2);
nt = stim_init_count(3);
flow = stim_init_count(4);
fhigh = stim_init_count(5);
fstep = stim_init_count(6);
nt = 2*nt +1;

asize(1) = nt;
asize(2) = nb;
w = hanning(nt);

% Read the overall cross-correlation
fid = fopen('stim_spike.avg','r');
stim_spike = fread(fid,asize,'double');
fclose(fid);
stim_spike = fliplr(stim_spike');
for ib=1:nb
   stim_spike(ib,:)=stim_spike(ib,:).*w';
end

% Read the JN cross-correlation
stim_spike_JN = zeros(nb, nt, nJN);
w = hanning(nt);
for iJN=1:nJN
	fname = sprintf('stim_spike_JN%d.avg',iJN);
	fid = fopen(fname,'r');
	for ib=1:nb
		stim_spike_JN(ib,:,iJN) = wrev(fread(fid,asize(1),'double')).*w;
	end
	fclose(fid);
end

% Calculate fft for all JN values.
stim_spike_JNf = fft(stim_spike_JN,[],2);
stim_spike_JNmf = mean(stim_spike_JNf,3);
stim_spikef = fft(stim_spike,[],2);

JNv = (nJN-1)*(nJN-1)/nJN;
j = sqrt(-1);
hcuttoff=0;
nf = (nt-1)/2 +1;
for ib=1:nb
	itstart = 1;
	itend = nf;
	below = 0;
	for it=1:nf
		stim_spike_JNvf(ib,it) = JNv*cov(permute(real(stim_spike_JNf(ib,it,:)),[3 2 1]))+ j*JNv*cov(permute(imag(stim_spike_JNf(ib,it,:)),[3 2 1]));
%		rmean = real(stim_spike_JNmf(ib,it));
		rmean = real(stim_spikef(ib,it));
      rstd = sqrt(real(stim_spike_JNvf(ib,it)));
%		imean = imag(stim_spike_JNmf(ib,it));
      imean = imag(stim_spikef(ib,it));
		istd = sqrt(imag(stim_spike_JNvf(ib,it)));
		if abs(rmean) < nstd_val*rstd & abs(imean) < nstd_val*istd
			if itstart == 1 
				itstart = it;
			end
			below = below + 1;
			if below == 3
				itend = it;
			end
		end
		stim_spike_sf(ib,it)= rmean + j*imean;
	end
	for it=1:nf
		if it > itstart 
            expval = exp(-0.5*(it-itstart)^2/(itend-itstart)^2);
            stim_spike_sf(ib,it) = stim_spike_sf(ib,it)*expval;
            stim_spike_JNf(ib,it,:) = stim_spike_JNf(ib,it,:).*expval;
% old way
			%if it < itend
				%stim_spike_sf(ib,it) = stim_spike_sf(ib,it)*(itend-it)/(itend-itstart);
				%stim_spike_JNf(ib,it,:) = stim_spike_JNf(ib,it,:).*((itend-it)/(itend-itstart));
            %else
				%stim_spike_sf(ib,it) = 0.0;
				%for iJN=1:nJN
					%stim_spike_JNf(ib,it,iJN) = 0.0;
                %end
            %end
		end
		if it ~= 1
			stim_spike_sf(ib,nt+2-it) = conj(stim_spike_sf(ib,it));
			stim_spike_JNf(ib,nt+2-it,:) = conj(stim_spike_JNf(ib,it,:));
		end
	end
end

% Reverse fft
stim_spike_JNm = real(ifft(stim_spike_JNmf,[],2));
stim_spike_s = real(ifft(stim_spike_sf,[],2));
stim_spike_JNs = real(ifft(stim_spike_JNf,[],2));

