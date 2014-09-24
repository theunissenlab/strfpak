function [infopre, infouppre, infodownpre, info, infoup, infodown,...
          cc_spike_pre_best, cc_spike_pre_constant, cc_two_halves_tmax,...
          cc_two_halves_constant, tmax_pre, cc_ratio_max, cc_ratio_const, tmax_ratio,...
          fpxypre, cxy, cxypre, cxyup, cxydown, cxyuppre, cxydownpre,...
	      cc_spike_pre_corrval, cc_spike_pre_corrval_low, cc_spike_pre_corrval_high,...
          cc_two_halves_corrval, cc_two_halves_corrval_low, cc_two_halves_corrval_high,raw_r]...
        = cal_Validate(this_est,spike_est1, spike_est2,ntol,istd,sample_rate,numTimeBin,running_in_script_mode) %        = cal_Validate(spike_est1, spike_est2,ntol,istd,sample_rate,numTimeBin,smoothVector,smoothConstant)
% function [infopre, infouppre, infodownpre, info, infoup, infodown,...
%          cc_spike_pre_best, cc_spike_pre_constant, cc_two_halves_tmax,...
%          cc_two_halves_constant, tmax_pre, cc_ratio_max, tmax_ratio,...
%          fpxy, cxy, cxypre, cxyup, cxydown, cxyuppre, cxydownpre,...
%          cc_spike_pre_corrval, cc_two_halves_corrval]...
%        = cal_Validate(ntol,sample_rate, smoothVector, smoothConstant)
%
%    -- Calculates the information between one trial and the predicted, 
%       infopre, and between one trial and the mean rate, info, 
%    Input:
%       ntol: index of the list of tol values
%       sample_rate: sampling rate in Hz
%      smoothVector: the vector of filter width
%      smoothConstant: constant filterwidth
%    Output:
%      infopre:   predicted information value
%      infouppre: upbound of predicted information value
%      infodownpre: lowbound of predicted information value
%      info: origianl information value
%      infoup: upbound of origianl information value
%      infodown: lowbound of origianl information value
%
%             STRFPAK: STRF Estimation Software
% Copyright �2003. The Regents of the University of California (Regents).
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

% Created by Anne and FET. 
% Modified by JXZ, Dec. 8, 2002
% Modified by JXZ, 7/14/2005
%    Add one more parameter 'istd'.
% Modified by JXZ, 8/22/2005
%    Load avg_spike1 and avg_spike2 one time.


% ==========================================================
%  Check if we have valid input
% ==========================================================
if ~exist('smoothVector')
    global smoothVect
    if isempty(smoothVect)
        smoothVector = [5 8 29];
    else
        smoothVector = smoothVect;
    end
end
global psth_smoothconst

if ~exist('smoothConstant','var')
    if isempty(psth_smoothconst)
        smoothConstant = 21;
    else
        smoothConstant = psth_smoothconst;
    end
end
% % find the indices of the ones that have the correct number of trials;
global outputPath
spike_psth = Check_And_Load(fullfile(outputPath,'spike_psth_count.mat'));

ntrials_proper=spike_psth(end);
spike_psth=spike_psth(1:end-1);
ntrials_index=find(spike_psth>=ntrials_proper);
ntrials = spike_psth;
%nrec = length(ntrials_index);
nrec = length(spike_psth>=ntrials_proper);
clear spike_psth;

% Read all the psth first to calculate the overall mean rate
spike_pre=[];

% Load predicted PSTH for tol val( ntol )
% spredresult = fullfile(outputPath,...
%            sprintf('predResult_EstSpike_Tol%d.mat', ntol));
% estSpike = Check_And_Load(spredresult);

   
% JXZ 7/18/2005
%  get time-varying mean rate
global timevary_PSTH
if timevary_PSTH ==1
    loadpsth = load(fullfile(outputPath, 'stim_avg.mat'), 'Avg_psth','constmeanrate');
    mean_rate = loadpsth.Avg_psth;
    const_meanrate = loadpsth.constmeanrate;
else
    loadpsth = load(fullfile(outputPath, 'stim_avg.mat'), 'constmeanrate');
    const_meanrate = loadpsth.constmeanrate;
end

% 8/10/05. Fet. The mean rate or time-varying mean rate is now added to the prediction
for k=1:nrec
    irec = k; %ntrials_index(k);
    % Load spike_pre
    spike_pretemp = this_est{k};   
    if timevary_PSTH == 1
        timedur1 = max(size(spike_pretemp));
        spike_pretemp = spike_pretemp + mean_rate(irec, 1:timedur1)';
        
    else
        spike_pretemp = spike_pretemp + const_meanrate;
    end   
    spike_pre=[spike_pre; spike_pretemp];
end

% Rectify the prediction for values less than 0
global  allow_negative_rates
if isempty(allow_negative_rates)
    allow_negative_rates = 0;
end
if ~allow_negative_rates
    spike_pre(find(spike_pre<0.0)) = 0.0;
end
% Call Signal Noise Ratio function to calculate predinfo value
% Set numTimeBin as 2% of time length
if ~exist('numTimeBin')
    numTimeBin =  round(0.02 * size(spike_est2, 1));

    if mod(numTimeBin, 2) == 1
        numTimeBin = numTimeBin +1;
    end

    if numTimeBin < 6
        numTimeBin = 8;
    end
end

% Modified by JXZ based on calc_est.m
%[infopre, infouppre, infodownpre, fpxypre,cxypre,cxyuppre, cxydownpre, cxy_notnormalizedpre,...
%    cxy_notnormalizeduppre,cxy_notnormalizedlopre]=SNRinfo(ntrials_proper,1, [(spike_est1+spike_est2)/2 spike_pre],...
%    numTimeBin, sample_rate);
%
%[info, infoup, infodown, fpxy,cxy,cxyup, cxydown, cxy_notnormalized,cxy_notnormalizedup,cxy_notnormalizedlo]=...
%    SNRinfo(ntrials_proper,2, [spike_est1 spike_est2], numTimeBin, sample_rate);

% New calculation of information value, April 09, 2003
% Call Signal Noise Ratio function to calculate info value

% Junli: 9/23/2004
%     Get correct ntrials from global variable assigned in read_SpikeTime.m
%
% global num_trials
% if isempty(num_trials)
%     num_trials = ntrials_proper;
% end

% JXZ 7/18/2005
% Avoid cutoff in SNRinfo calculation: Information is estimated for all
% frequencies where the coherence is significantly different from zero. The
% older routine SNR_info zero out all high frequencies after reaching a
% frequency with non significant coherence.  Note that the mean is
% subracted because the absolute value of information obtained from the
% multi-taper method is affected by the DC.
% size(spike_pre-mean(spike_pre))
% size((spike_est1+spike_est2)/2-const_meanrate)
% size(spike_est1-const_meanrate)
% size(spike_est2-const_meanrate)
lsp = length(spike_pre);
lse1 = length(spike_est1);
lse2 = length(spike_est2);
the_max = max([lsp,lse1,lse2]);
the_min = min([lsp,lse1,lse2]);
if (the_max-the_min) == 1  %If there's a 1 ms fault, it can be because of even/odd window sizes, and we don't care.
    spike_pre = spike_pre(1:the_min);
    spike_est1 = spike_est1(1:the_min);
    spike_est2 = spike_est2(1:the_min);
    %disp('trimming spike_pre to make it fit');
end
global running_in_script_mode
if 0%strcmp(running_in_script_mode,'yes')
    [infopre, infouppre, infodownpre, fpxypre,cxypre,cxyuppre, cxydownpre, cxy_notnormalizedpre,cxy_notnormalizeduppre,cxy_notnormalizedlopre,cxy,cxyup, cxydown,info, infoup, infodown]=...
    SNRinfo_no_cutoff(ntrials_proper,3, [spike_pre-mean(spike_pre) (spike_est1+spike_est2)/2-const_meanrate spike_est1-const_meanrate spike_est2-const_meanrate], numTimeBin, sample_rate);
else
[infopre, infouppre, infodownpre, fpxypre,cxypre,cxyuppre, cxydownpre, cxy_notnormalizedpre,cxy_notnormalizeduppre,cxy_notnormalizedlopre,cxy,cxyup, cxydown,info, infoup, infodown]=...
    do_locally_cached_calc(get_local_cache_dir,'SNRinfo_no_cutoff',ntrials_proper,3, [spike_pre-mean(spike_pre) (spike_est1+spike_est2)/2-const_meanrate spike_est1-const_meanrate spike_est2-const_meanrate], numTimeBin, sample_rate);
end
% JXZ 9/7/05: cutoff due to statistical error
if info < infopre
    if infodownpre <= infoup
        infopre = info;
    else
        %warndlg('Warning: Predicted information value is higher than input info value.');
        disp('Warning: infopre is greater than info in cal_Validate.m due to too noisy data.')
    end
end
 
% Smoothing parameters when calculating CC
widthvector=smoothVector;
width=widthvector(1):widthvector(2):widthvector(3);

% Call cc_plot_est to calculate Correlation Coefficient and then normalize it
%
if 0 %strcmp(running_in_script_mode,'yes')
    [cc_two_halves_best1, cc_two_halves_best1_low, cc_two_halves_best1_high, cc_two_halves_constant1,...
       cc_two_halves_constant1_low, cc_two_halves_constant1_high, cc_two_halves_corrval1,...
       cc_two_halves_corrval1_low, cc_two_halves_corrval1_high, tmax_two_halves]...
    = cc_plot_est(spike_est1, spike_est2, widthvector, smoothConstant, 1);

else
[cc_two_halves_best1, cc_two_halves_best1_low, cc_two_halves_best1_high, cc_two_halves_constant1,...
       cc_two_halves_constant1_low, cc_two_halves_constant1_high, cc_two_halves_corrval1,...
       cc_two_halves_corrval1_low, cc_two_halves_corrval1_high, tmax_two_halves]...
    = do_locally_cached_calc(get_local_cache_dir,'cc_plot_est',spike_est1, spike_est2, widthvector, smoothConstant, 1);
end

% Deal with zero cc's also and add estimates of lower and upper bound...
if ( cc_two_halves_constant1 ~= 0.0 )
    temp=(-ntrials_proper+ntrials_proper*sqrt(1/cc_two_halves_constant1))/2;
    cc_two_halves_constant=1./(temp+1);
else
    cc_two_halves_constant=0;
end

if ( cc_two_halves_constant1_low ~= 0.0 )
    temp=(-ntrials_proper+ntrials_proper*sqrt(1/cc_two_halves_constant1_low))/2;
    cc_two_halves_constant_low=1./(temp+1);
else
    cc_two_halves_constant_low=0;
end

if ( cc_two_halves_constant1_high ~= 0.0 )
    temp=(-ntrials_proper+ntrials_proper*sqrt(1/cc_two_halves_constant1_high))/2;
    cc_two_halves_constant_high=1./(temp+1);
else
    cc_two_halves_constant_high=0;
end


if ( cc_two_halves_best1 ~= 0.0 )
    temp=(-ntrials_proper+ntrials_proper*sqrt(1/cc_two_halves_best1))/2;
    cc_two_halves_best=1./(temp+1);
else
    cc_two_halves_best=0;
end

notzeroind = find(cc_two_halves_corrval1 ~= 0);
temp = zeros(size(cc_two_halves_corrval1));
temp(notzeroind) =(-ntrials_proper+ntrials_proper*sqrt(1./cc_two_halves_corrval1(notzeroind)))/2;
cc_two_halves_corrval = zeros(size(cc_two_halves_corrval1));
cc_two_halves_corrval(notzeroind)=1./(temp(notzeroind)+1);

clear notzeroind temp
notzeroind = find(cc_two_halves_corrval1_low ~= 0);
temp = zeros(size(cc_two_halves_corrval1_low));
temp(notzeroind) =(-ntrials_proper+ntrials_proper*sqrt(1./cc_two_halves_corrval1_low(notzeroind)))/2;
cc_two_halves_corrval_low = zeros(size(cc_two_halves_corrval1_low));
cc_two_halves_corrval_low(notzeroind)=1./(temp(notzeroind)+1);

clear notzeroind temp
notzeroind = find(cc_two_halves_corrval1_high ~= 0);
temp = zeros(size(cc_two_halves_corrval1_high));
temp(notzeroind) =(-ntrials_proper+ntrials_proper*sqrt(1./cc_two_halves_corrval1_high(notzeroind)))/2;
cc_two_halves_corrval_high = zeros(size(cc_two_halves_corrval1_high));
cc_two_halves_corrval_high(notzeroind)=1./(temp(notzeroind)+1);

cc_two_halves_corrval=sqrt(cc_two_halves_corrval);
cc_two_halves_corrval_low=sqrt(cc_two_halves_corrval_low);
cc_two_halves_corrval_high=sqrt(cc_two_halves_corrval_high);
cc_two_halves_best=sqrt(cc_two_halves_best);
cc_two_halves_constant=sqrt(cc_two_halves_constant);

% Calculate raw r between raw psth and predicted psth
raw_r = diag(corrcoef((spike_est1+spike_est2)/2,spike_pre), 1);

% Call cc_plot_est to calculate predicted Correlation Coefficient 
% and then normalize it
% Calculate upper and lower bounds using the expected value obtained from
% cc_two_halves_best1
if 0%strcmp(running_in_script_mode,'yes')
    [cc_spike_pre_best1, cc_spike_pre_best1_low, cc_spike_pre_best1_high, cc_spike_pre_constant1,...
       cc_spike_pre_constant1_low, cc_spike_pre_constant1_high, cc_spike_pre_corrval1,...
       cc_spike_pre_corrval1_low, cc_spike_pre_corrval1_high, tmax_pre]...
    =cc_plot_est((spike_est1+spike_est2)/2,spike_pre,widthvector, smoothConstant , 0);

else
[cc_spike_pre_best1, cc_spike_pre_best1_low, cc_spike_pre_best1_high, cc_spike_pre_constant1,...
       cc_spike_pre_constant1_low, cc_spike_pre_constant1_high, cc_spike_pre_corrval1,...
       cc_spike_pre_corrval1_low, cc_spike_pre_corrval1_high, tmax_pre]...
    =do_locally_cached_calc(get_local_cache_dir,'cc_plot_est',(spike_est1+spike_est2)/2,spike_pre,widthvector, smoothConstant , 0);
end

if ( cc_two_halves_best1 == 0 )
    cc_spike_pre_best = 0;
else
    cc_spike_pre_best=cc_spike_pre_best1*(1+sqrt(1/cc_two_halves_best1))/(-ntrials_proper+ntrials_proper*sqrt(1/cc_two_halves_best1)+2);
end

if ( cc_two_halves_constant1 == 0 )
    cc_spike_pre_constant = 0;
else
    cc_spike_pre_constant=cc_spike_pre_constant1*(1+sqrt(1/cc_two_halves_constant1))/(-ntrials_proper+ntrials_proper*sqrt(1/cc_two_halves_constant1)+2);
end
if ( cc_two_halves_constant1_low == 0 )
    cc_spike_pre_constant_low = 0;
else
    cc_spike_pre_constant_low=cc_spike_pre_constant1_low*(1+sqrt(1/cc_two_halves_constant1_low))/(-ntrials_proper+ntrials_proper*sqrt(1/cc_two_halves_constant1_low)+2);
end
if ( cc_two_halves_constant1_high == 0 )
    cc_spike_pre_constant_high = 0;
else
    cc_spike_pre_constant_high=cc_spike_pre_constant1_high*(1+sqrt(1/cc_two_halves_constant1_high))/(-ntrials_proper+ntrials_proper*sqrt(1/cc_two_halves_constant1_high)+2);
end

clear notzeroind
notzeroind = find(cc_two_halves_corrval1 ~=0 );
cc_spike_pre_corrval = zeros(size(cc_spike_pre_corrval1));
cc_spike_pre_corrval(notzeroind)=cc_spike_pre_corrval1(notzeroind).*(1+sqrt(1./cc_two_halves_corrval1(notzeroind)))./(-ntrials_proper+ntrials_proper*sqrt(1./cc_two_halves_corrval1(notzeroind))+2);

clear notzeroind
notzeroind = find(cc_two_halves_corrval1_low ~=0 );
cc_spike_pre_corrval_low = zeros(size(cc_spike_pre_corrval1_low));
cc_spike_pre_corrval_low(notzeroind)=cc_spike_pre_corrval1_low(notzeroind).*(1+sqrt(1./cc_two_halves_corrval1_low(notzeroind)))./(-ntrials_proper+ntrials_proper*sqrt(1./cc_two_halves_corrval1_low(notzeroind))+2);

clear notzeroind
notzeroind = find(cc_two_halves_corrval1_high ~=0 );
cc_spike_pre_corrval_high = zeros(size(cc_spike_pre_corrval1_high));
cc_spike_pre_corrval_high(notzeroind)=cc_spike_pre_corrval1_high(notzeroind).*(1+sqrt(1./cc_two_halves_corrval1_high(notzeroind)))./(-ntrials_proper+ntrials_proper*sqrt(1./cc_two_halves_corrval1_high(notzeroind))+2);

cc_spike_pre_corrval=sqrt(cc_spike_pre_corrval);
cc_spike_pre_corrval_low=sqrt(cc_spike_pre_corrval_low);
cc_spike_pre_corrval_high=sqrt(cc_spike_pre_corrval_high);

cc_spike_pre_best=sqrt(cc_spike_pre_best);
cc_spike_pre_constant=sqrt(cc_spike_pre_constant);


tindex=find(width==tmax_pre);
cc_two_halves_tmax=cc_two_halves_corrval(tindex);

% Find the best ratio of predicted_CC and CC
cc_ratio = zeros(size(cc_spike_pre_corrval));
notzeroind = find(cc_two_halves_corrval~=0);
cc_ratio(notzeroind) = cc_spike_pre_corrval(notzeroind)./cc_two_halves_corrval(notzeroind);
for icc=1:length(notzeroind)
   if (cc_ratio(icc) > 1.0) 
       if (cc_spike_pre_corrval_low(icc) < cc_two_halves_corrval_high(icc) )
           cc_ratio(icc) = 1.0;
       else
           %print a warning of cc ratio greater than one...
           disp('Warning: cc_ratio is greater than one in cal_Validate.m due to too noisy data.');
       end
   end
end

% Find the constant ratio of cc_predicted and cc
if cc_two_halves_constant==0
    cc_ratio_const = 0.0;
else
    cc_ratio_const = cc_spike_pre_constant ./cc_two_halves_constant;
    if cc_ratio_const > 1.0
        if (cc_spike_pre_constant_low < cc_two_halves_constant_high)
            cc_ratio_const = 1.0;
        else
            %warndlg('Warning: cc_ratio_const is greater than one.');
            disp('Warning: cc_ratio_const is greater than one in cal_Validate.m due to too noisy data.')
        end
    end
end

cc_ratio_max = max(cc_ratio);
tindex = find(cc_ratio == cc_ratio_max);
tmax_ratio = width(tindex);

% ======================================================================
% END of cal_Validate
% ======================================================================


