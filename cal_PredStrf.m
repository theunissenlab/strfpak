function [ErrorFlag, est_spike, avg_spike1, avg_spike2] =...
    cal_PredStrf(running_flag,predDS, stim_avg, strfFiles, ntol, twindow,...
    nband, running_in_script_mode)
%
% [est_spike, avg_spike1, avg_spike2] = cal_PreStrf(predDS,
%         stim_avg, twindow,nband,matchflg,end_window)
%    -- calculate pred neuron response based on stimulus data files
%       and estimated STRF. We use the same calc. parameters as
%       before and convolve stimulus and strf and get pred. spike.
%       In order to calculate how good we fit, we also calculate
%       avg_spike1 and avg_spike2 based on neuron spike data files.
%   Input:
%       predDS:  the cell of each data struct that contains four fields:
%               stimfiles  - stimulus file name
%               respfiles  - response file name
%               nlength    - length of time domain
%               ntrials    - num of trials
%              e.g. predDS{1} = struct('stimfiles', 'stim1.dat', 'respfiles',
%                    'resp1.dat', 'nlength', 1723, 'ntrials', 20);
%       stim_avg: avg stimulus that used to smooth the noise
%                 If stim_avg is empty, we will call cal_AVG to get it.
%       twindow: the variable to set the time interval to calculate
%                autocorrelation. e.g. twindow=[-TimeLag TimeLag]
%       nband: the size of spatio domain of the stimulus file
%       matchflg: the flag to show which STRF we used to do prediction.
%                 If matchflg is 1, we used JN-version STRF. 0 uses
%                 STRF to check the good match(overestimated).
%       end_window: the time interval which dont count for data analysis
%       ntol: the number of tol_value for computing prediction.
%   Output:
%       est_spike: predicted neuron spike
%       avg_spike1: avg of half real neuron spike.
%       avg_spike2: avg of the rest half real neuron spike.
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

% Created by JXZ
% Sep. 27, 2002
% Modified by JXZ, 7/13/2005
%  1. Add prediction on filtered STRFs.
%
% Modified by JXZ, 8/19/2005
%  1. Remove identical saved files predResult_AvgSpike1&2.mat
%  2. Remove mean_constant saved files since it can be loaded from
%     stim_avg.mat.

% ========================================================
% check if we have valid argument
% ========================================================
ErrorFlag = 0;
if ~exist('matchflg')
    global matchflg
    if isempty(matchflg)
        errordlg('Please set match_flag first to do prediction',...
            'Global Variable Error', 'modal')
        ErrorFlag = 1;
        return;
    end
end

if ~exist('end_window')
    end_window = 0;
end

if isempty(predDS)
    errordlg('There are no predicted files','Variable Error', 'modal')
    ErrorFlag = 1;
    return
end

% ========================================================
% initialize the results and get size of xy axis
% ========================================================
filecount = length(predDS);

% temporal axis range
tot_corr = diff(twindow) + 1;

% spatial axis range
%spa_corr = nband;
global outputPath
if 1%~exist(fullfile(outputPath, 'spike_psth_count.mat'),'file')
    % ========================================================
    % find the most common number of trials based on predDS
    % ========================================================
    spike_psth = zeros(1, filecount);
    for ii=1:filecount
        spike_psth(1,ii) = predDS{ii}.ntrials;
        %fprintf(fpsth, '%d\n', spike_psth(1,i));
    end

    minTrials = min(spike_psth);
    maxTrials = max(spike_psth);
    rangeTrials = minTrials:maxTrials;

    hist_trialcount = histc(spike_psth, rangeTrials);
    [maxTrial, maxIdx] = max(hist_trialcount);
    %current_max_val = rangeTrials(maxIdx);

    %     % Only predict on the same trial number data
    %     ntrials_index=find(spike_psth>=current_max_val)
    ntrials_index = 1:length(spike_psth);
    current_max_val = 'not used';
    % =======================================================
    % save the psth_count file into .mat format
    % =======================================================
    spike_psth(1, filecount+1) = rangeTrials(maxIdx);

    if isempty(outputPath)
        disp('Saving output to Output Dir.');
        stat = mkdir('Output');
        outputPath = fullfile(pwd, 'Output');
    end

    save(fullfile(outputPath, 'spike_psth_count.mat'), 'spike_psth');
end

% ========================================================
%   See if the calculation has to be done.
% ========================================================
[cached_dir,maxsize] = dir_of_caches;
cache_crosscorr = ~strcmp(cached_dir,'null');
global Std_val
use_more_memory = 1;
if use_more_memory
    preload_stims = make_ready_stims(predDS,stim_avg);
else
    make_ready_stims(predDS,stim_avg);
end
global now_do_untouched best_std_index

if strcmp(now_do_untouched,'yes')
    kk_index = best_std_index;
else
    kk_index = 1:length(Std_val);
end
if ~exist('running_in_script_mode','var')
    global running_in_script_mode
end
for kk = kk_index
    the_Std = Std_val(kk);
    strf_checksum = get_strf_checksum(strfFiles);
    if cache_crosscorr% & ~strcmp(running_in_script_mode,'yes') Use a local cache even in batch mode now.
        hashes_of_pred_stims = create_stim_cache_file(outputPath,predDS);
        the_checksum = checksum(load_function_text('cal_cache_PredStrf_2'),now_do_untouched,strf_checksum,hashes_of_pred_stims,the_Std);  %  Note: psth_option, the actual tol value, etc. are all implicit in the STRF checksum.
        if use_more_memory
            [ErrorFlag, est_spike_to_shuffle{kk}, avg_spike1, avg_spike2] =...
                do_locally_cached_calc_checksum_known(get_local_cache_dir,'cal_cache_PredStrf_2',the_checksum,cache_crosscorr,the_checksum,running_flag,predDS, stim_avg, strfFiles, ntol, twindow,...
                matchflg, end_window,the_Std,ntrials_index,current_max_val,preload_stims);
        else
            [ErrorFlag, est_spike_to_shuffle{kk}, avg_spike1, avg_spike2] =...
                do_locally_cached_calc_checksum_known(get_local_cache_dir,'cal_cache_PredStrf_2',the_checksum,cache_crosscorr,the_checksum,running_flag,predDS, stim_avg, strfFiles, ntol, twindow,...
                matchflg, end_window,the_Std,ntrials_index,current_max_val);
        end


    else
        the_checksum = 'not_used';
        %the_checksum = checksum(strf_checksum,hashes_of_pred_stims,the_Std);  %  Note: psth_option, the actual tol value, etc. are all implicit in the STRF checksum.
        if use_more_memory
            [ErrorFlag, est_spike_to_shuffle{kk}, avg_spike1, avg_spike2] =...
                cal_cache_PredStrf_2(cache_crosscorr,the_checksum,running_flag,predDS, stim_avg, strfFiles, ntol, twindow,...
                matchflg, end_window,the_Std,ntrials_index,current_max_val,preload_stims);
        else
            [ErrorFlag, est_spike_to_shuffle{kk}, avg_spike1, avg_spike2] =...
                cal_cache_PredStrf_2(cache_crosscorr,the_checksum,running_flag,predDS, stim_avg, strfFiles, ntol, twindow,...
                matchflg, end_window,the_Std,ntrials_index,current_max_val);
        end

    end
end
for ii = 1:length(est_spike_to_shuffle)
    len = length(est_spike_to_shuffle{ii});
    for jj = 1:len
        est_spike{jj}{ii} = est_spike_to_shuffle{ii}{jj};
    end
end



% ========================================================
% Save the results to the file with .mat format
% ========================================================
currentPath = pwd;

if ~isempty(outputPath)
    cd (outputPath);
else
    disp('Saving output to Output Dir.');
    stat = mkdir('Output');
    cd('Output');
end

%mean_rate = mean(spike_est_all);
global now_do_untouched
if strcmp(now_do_untouched,'yes')
    filename = sprintf('predResult_EstSpike_untouched.mat');
else
    filename = sprintf('predResult_EstSpike_Tol%d.mat', ntol);
end
save(filename, 'est_spike');

% Only save one time avg_spike1 and avg_spike2 since they are identical for
%  all tol values.
%avg_spike1
global now_do_untouched matchflg

if (ntol==1) | ((strcmp(now_do_untouched,'yes')&matchflg == 0))
    for k = 1:length(avg_spike1)     % loop through all data files
        avgspike1{k} = avg_spike1{k}';
        avgspike2{k} = avg_spike2{k}';
    end
    if strcmp(now_do_untouched,'yes')
        save('finalResult_avgSpike1.mat', 'avgspike1');
        %filename = sprintf('predResult_avgSpike2_Tol%d.mat', ntol);
        save('finalResult_avgSpike2.mat', 'avgspike2');

    else

        %filename = sprintf('predResult_avgSpike1_Tol%d.mat', ntol);
        save('predResult_avgSpike1.mat', 'avgspike1');
        %filename = sprintf('predResult_avgSpike2_Tol%d.mat', ntol);
        save('predResult_avgSpike2.mat', 'avgspike2');
    end
end
% filename = sprintf('predResult_meanrate_Tol%d.mat', ntol);
% save(filename, 'mean_rate');

cd(currentPath);

% ========================================================
%  END of CAL_PREDSTRF
% ========================================================
