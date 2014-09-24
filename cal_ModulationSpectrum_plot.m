function  cal_ModulationSpectrum_plot(...
               fwidthHz, amp_samplerate, fstep)
%
%             STRFPAK: STRF Estimation Software
% Copyright ©2003. The Regents of the University of California (Regents).
% All Rights Reserved.
% Created by Theunissen Lab and Gallant Lab, Department of Psychology, Un
% -iversity of California, Berkeley.
%
% Permission to use, copy, and modify this software and its documentation
% for educational, research, and not-for-profit purposes, without fee and
% without a signed licensing agreement, is hereby granted, provided that
% the above copyright notice, this paragraph and the following two paragr
% -aphs appear in all copies and modifications. Contact The Office of Tec
% -hnology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510,
% Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing
% opportunities.
%
%IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
%SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
%ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
%REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
%LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
%PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
%PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PRO
%-VIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

% Crated and Modified by FET and JXZ, 2003.

% ===========================================================
% Load stim auto-correlation matrix
% ===========================================================
global outputPath
if isempty(outputPath)
    errordlg('Where do you want to load your autocorr matrix?')
    return
end

autoCorrfile = fullfile(outputPath, 'Stim_autocorr.mat');
stim_autocorr = Check_And_Load(autoCorrfile);

%
% Get dimension size of modulation spectrum from stim_autocorr 
%  
ncorr = size(stim_autocorr, 1);
nt = size(stim_autocorr, 2);
nb = (-1 + sqrt(1+8*ncorr))/2;

% stim_modspec has size 2*nb-1 because of positive and negative db */
%stim_modspec_temp=zeros(2*nb-1,nt);
stim_modspec =zeros(2*nb-1,nt);

% Take average of stim along space
stimdiag = zeros(nb, nt);
for xboff=1:nb
    xb=xboff;
    for i=1:nb-xboff+1     
	    %stim_mean(xboff,i) = sqrt(mean(stim(xb,:)));
        stimdiag(xboff,:) = stimdiag(xboff,:)+stim_autocorr(xb,:);
	    xb = xb + (nb-i+1);
    end
    stimdiag(xboff, :) = stimdiag(xboff, :)/(nb-xboff+1);
end

% Stuff a symmetric correlation by repeating the frequency.
for xb=1:nb
    stim_modspec(nb-1+xb,:) = stimdiag(xb,:); 
    
    if ( xb ~= 1)
        stim_modspec(nb+1-xb,:)=fliplr(stimdiag(xb,:));
        
    end   
end

% plot intermediate results
%
% figure;
% imagesc(stim_modspec); axis xy;

% Smoothing purpose, define hanning window for windowing
% Get a hanning window for windowing in time and frequency;
wht = hanning(nt);
whf = hanning(2*nb-1);

% Make a 2-d window
whtf = zeros(size(stim_modspec));
for it=1:nt
    whtf(:,it)=whf;
end
for ib=1:2*nb-1
    whtf(ib,:) = whtf(ib,:).*wht';
end

% Subtract DC value and window
stim_modspec0 = stim_modspec-mean(mean(stim_modspec));
stim_modspec0w = stim_modspec0.*whtf;
stim_modspecw = stim_modspec.*whtf;

% Next take 2d fft of stim_modspec
stim_modspec_f   =fft2(stim_modspec);
stim_modspec0w_f =fft2(stim_modspec0w);
stim_modspec0_f  =fft2(stim_modspec0);
stim_modspecw_f  =fft2(stim_modspecw);

% Prepare for plotting
stim_modspec_f = fftshift(abs(stim_modspec_f));
stim_modspec0w_f = fftshift(abs(stim_modspec0w_f));
stim_modspec0_f = fftshift(abs(stim_modspec0_f));
stim_modspecw_f = fftshift(abs(stim_modspecw_f));

% Find labels for x and y axis
% fstep is the separation between frequency bands
for i=1:2*nb-1
   dwf(i)= (i-nb)*(1/(2*fstep*(nb-1)));
end
% 1 ms (1000 on the numerator) is the sampling rate
ntt=(nt+1)/2;
for i=1:2*ntt-1
   dwt(i) = (i-ntt)*(amp_samplerate/(2*(ntt-1)));
end


% Ready for showing modulation spectrum at zero mean and windowed 
% in log coordinates
%

% Ask user what fraction of the total power they want to draw
% junli - 5/7/04

prompt = {'Enter the fraction of the total power to draw contours:'};
def = {'0.5 0.6 0.7 0.8 0.9'};
dlgTitle='Input for fraction of the power';
lineNo = 1;
fracP = inputdlg(prompt, dlgTitle,lineNo,def);

% parse the input string to a set of data files
% delimited by ' '
[stoken, rem] = strtok(char(fracP), ' \t,');
if isempty(stoken)
    errordlg('Please enter valid fraction_power.', 'Input Error', 'modal')
    return;
end
tolval = str2double(stoken);
if isnan(tolval) | tolval < 0
    errordlg('Please enter valid fraction_power.', 'Input Error', 'modal')
    return;
end

fracpower = [tolval];

while  ~isempty(rem)
    [stoken, rem] = strtok(rem, ' \t,');
    if isempty(stoken)
        break
    end
    tolval = str2double(stoken);
    if isnan(tolval) | tolval < 0
        errordlg('Please enter valid fraction_power.', 'Input Error', 'modal')
        return;
    end
     fracpower = [fracpower tolval];
    end

figure;
stim_modspec0w_f = calc_zero(dwt, dwf, fwidthHz, stim_modspec0w_f);
imagesc(dwt,dwf*1000,log(stim_modspec0w_f+1));
%shading interp;
%lighting phong;
axis square;
colormap(jet);
axis xy;
hold on;

%fracpower = [0.5 0.6 0.7 0.8 0.9];
fracvalues = calc_contour_values((stim_modspec0w_f), fracpower);

if length(fracvalues) == 1
    [C h] = contour(dwt,dwf*1000,log(stim_modspec0w_f+1),[fracvalues fracvalues],'k-');
else
    [C h] = contour(dwt,dwf*1000,log(stim_modspec0w_f+1),fracvalues,'k-');
end

hold off;
axis([-100 100 0 2]);
title('Modulation Spectrum');
xlabel('{\omega}_{t}(Hz)');
ylabel('{\omega}_{x}(Cycles/kHz)');

