function [ErrorFlag, est_spike, avg_spike1, avg_spike2] =...
    cal_cache_PredStrf_2(cache_crosscorr,the_checksum,running_flag,predDS, stim_avg, strfFiles, ntol, twindow,...
    matchflg, end_window,the_Std,ntrials_index,current_max_val,preload_stims);
%  This is a split of cal_PredStrf, with a chunk of calculations which can
%  be efficiently cached.
stdfilt = the_Std;
global outputPath now_do_untouched
if strcmp(now_do_untouched,'yes')&(matchflg == 0)

    avgspike_done = 0;
    if exist(fullfile(outputPath,'finalResult_avgSpike1.mat'),'file') & exist(fullfile(outputPath,'finalResult_avgSpike2.mat'),'file')
        load(fullfile(outputPath,'finalResult_avgSpike1.mat'));  % puts avgspike1 into workspace
        load(fullfile(outputPath,'finalResult_avgSpike2.mat'));  % puts avgspike2 into workspace
        avgspike_done = 1;
    end
else
    avgspike_done = 0;
    if exist(fullfile(outputPath,'predResult_avgSpike1.mat'),'file') & exist(fullfile(outputPath,'predResult_avgSpike2.mat'),'file')
        load(fullfile(outputPath,'predResult_avgSpike1.mat'));  % puts avgspike1 into workspace
        load(fullfile(outputPath,'predResult_avgSpike2.mat'));  % puts avgspike2 into workspace
        avgspike_done = 1;
    end
end

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

% ========================================================
%  Get the previous calculated STRF
% ========================================================
if exist(strfFiles, 'file') == 2
    % strf_result = Check_And_Load(strfFiles);
    strf_result = load(strfFiles);
else
    errordlg('There are no result files yet.',...
        'Display Error', 'modal')
    ErrorFlag = 1;
    return
end

% ========================================================
% After loading the result file, we assign the values to
% forward based on ntol input.
% ========================================================



if matchflg == 1    % read in forwardJN.filt
    forward = strf_result.STRFJN_Cell;
else                % read forward.filt directly
    forward = strf_result.STRF_Cell;
end                 % end of matchflg
forwardJNstd = strf_result.STRFJNstd_Cell;

% ========================================================
% JXZ, 7/13/2005
%    Add a loop for going through all the Std_val
% ========================================================
global Std_val Tol_val
nstd = length(Std_val);
totTol = length(Tol_val);
% %  Moved to cal_PredStrf
% % ========================================================
% % find the most common number of trials based on predDS
% % ========================================================
% spike_psth = zeros(1, filecount);
% for ii=1:filecount
%     spike_psth(1,ii) = predDS{ii}.ntrials;
%     %fprintf(fpsth, '%d\n', spike_psth(1,i));
% end
%
% minTrials = min(spike_psth);
% maxTrials = max(spike_psth);
% rangeTrials = minTrials:maxTrials;
%
% hist_trialcount = histc(spike_psth, rangeTrials);
% [maxTrial, maxIdx] = max(hist_trialcount);
% current_max_val = rangeTrials(maxIdx);
%
% % Only predict on the same trial number data
% ntrials_index=find(spike_psth>=current_max_val);
% % =======================================================
% % save the psth_count file into .mat format
% % =======================================================
% spike_psth(1, filecount+1) = current_max_val;
%
% if isempty(outputPath)
%     disp('Saving output to Output Dir.');
%     stat = mkdir('Output');
%     outputPath = fullfile(pwd, 'Output');
% end
%
% save(fullfile(outputPath, 'spike_psth_count.mat'), 'spike_psth');

% Oct 24: Junli add to calculate the overall mean rate
% spike_est_all = [];

% ========================================================
% Now caluclate the predicted Spike
% ========================================================

hashes_of_pred_stims = create_stim_cache_file(outputPath,predDS);


for k = 1:length(predDS)     % loop through all data files
    if exist('preload_stims','var')
        stimval = preload_stims{k};
    else
        filename = predDS{k}.stimfiles;
        [the_dir,filename] = fileparts(filename);
        global outputPath
        outfilename = fullfile(outputPath,[filename '_mean_removed.mat']);
        load(outfilename); %Puts "stimval" into the workspace
    end
    %nlen = predDS{fidx}.nlen;
    nlen = size(stimval,2);
    fidx = ntrials_index(k);
    if matchflg == 1
        % Compute the filtered JN STRF
        n_pred_stims = size(forward,3);
        if strcmp(now_do_untouched,'yes')
            index_to_use = sort(1 + mod([1:(n_pred_stims-2)]+fidx-1, n_pred_stims));
        else
            index_to_use = sort(1 + mod([1:(n_pred_stims-2)]+fidx, n_pred_stims));
        end

        forwardJN_s(:,:) = fast_filter_filter(mean(forward(:,:,index_to_use),3), forwardJNstd(:,:,1 + mod(fidx-1,n_pred_stims)), stdfilt);
    else
        forwardJN_s(:,:) = fast_filter_filter(forward, squeeze(mean(forwardJNstd,3)), stdfilt);
    end

    this_stim_hash = hashes_of_pred_stims{fidx};
    %   Let's not load it just yet - we might not have to load them
    %     % load stimulus files
    %     stim_env = Check_And_Load(predDS{fidx}.stimfiles);
    %


    global psth_len;   % Temp variable
    if 1%(predDS{fidx}.ntrials >= current_max_val)
        % assign avg_spike1 and avg_spike2
        % here avg_spike1 are mean value of real spike that are taken
        % from each other. e.g. trial_1, trial_2, trial_3, trial_4...
        % avg_spike1 = (trial_1 + trial_3+..)/tot/2
        % Just save avg_spike1 and avg_spike2 one time
        if ~avgspike_done %(~exist('avgspike1','var')) | (~exist('avgspike2','var'))%ntol == 1
            evenTrial = 0;
            oddTrial = 0;
            global respsamprate ampsamprate
            % load response files
            rawResp = Check_And_Load(predDS{fidx}.respfiles);

            if iscell(rawResp)
                spiketrain = zeros(predDS{fidx}.ntrials,predDS{fidx}.nlen);
                for trial_ind =1:predDS{fidx}.ntrials

                    spiketrain(trial_ind, rawResp{trial_ind}) = ones(1,length(rawResp{trial_ind}));
                end
                newpsth = resample(spiketrain', ampsamprate, respsamprate);

                newpsth = newpsth'; % make sure new response data is trials x T.
                newpsth(find(newpsth < 0)) = 0;
                psth_rec = newpsth; 
            else
                psth_rec = rawResp;
            end

            nlen = min(size(psth_rec,2), nlen);
            psth_len{fidx} = nlen;
            % get neuron spike from data
            avg_spike1{k} = zeros(1, nlen+end_window);
            avg_spike2{k} = zeros(1, nlen+end_window);
            for i=1:size(psth_rec,1)
                if mod(i, 2) == 1   % in order to consistent with dcp_forward.c
                    avg_spike2{k} = avg_spike2{k} + psth_rec(i, 1:nlen);
                    oddTrial = oddTrial +1;
                else
                    avg_spike1{k} = avg_spike1{k} + psth_rec(i, 1:nlen);
                    evenTrial = evenTrial +1;
                end
            end

            % JXZ, 9/14/2005
            % Make sure if the user only has one trial data
            %if current_max_val >1
            if i>1

                avg_spike1{k} = avg_spike1{k}/evenTrial;
                avg_spike2{k} = avg_spike2{k}/oddTrial;
            else
                avg_spike1{k} = avg_spike2{k};
            end
            avgspike1{k} = avg_spike1{k}';
            avgspike2{k} = avg_spike2{k}';
        else
            avg_spike1{k} = avgspike1{k}';
            avg_spike2{k} = avgspike2{k}';
        end
        % subtract mean of stim from each stim
        nlen = length(avg_spike1{k});

        % =====================================================
        % JXZ 7/13/2005
        %  Add std_val loop


        % Patrick: 9/26/2005
        % Move fast_filter_filter outside ib loop to improve the
        % preformance
        stdfilt = the_Std;
        est_spike_k = zeros(size(stimval,2),1);

        %  This algorithm for computing the predictions is equivalent to
        %  the one in "predict_one_tol_std_combo_2", but it's faster.
        %  Multiplication of an M X M matrix with another is NOT as hard as
        %  O(M^3); it can be done in O(M^(log2(7))), which makes a big
        %  difference (a factor of 4, we find) when the stim is big.
        twin = (size(forwardJN_s,2) - 1)/2;
        c_t = forwardJN_s(:,end:-1:1)'*stimval;  %  This makes a matrix of dot products of the filter with the stim.
        pad = zeros(2+ 2*twin + 2*twin^2,1);
        c_t = [pad; c_t(:) ; pad];
        template = (-twin:twin)*(2*twin+2) -twin + length(pad);  %  This is how the matrix adds up to get what we want.
        for jj = 1:size(stimval,2)
            %jj*(2*twin+1) + template - length(pad)
            est_spike_k(jj) = sum(c_t((jj)*(2*twin+1) + template));
        end


        %
        %         if cache_crosscorr
        %
        %             this_checksum = checksum(the_checksum,this_stim_hash,matchflg);  % sensitive on the STRF, the dataset and the STD value.
        %
        %             est_spike_k = do_cached_calc_checksum_known('predict_one_tol_std_combo_2',this_checksum,forwardJN_s,nlen,end_window,fidx,nband,stim_avg);
        %         else
        %             est_spike_k = predict_one_tol_std_combo_2(forwardJN_s,nlen,end_window,fidx,nband,stim_avg);
        %         end
        est_spike{k} = est_spike_k;
        % END of istd
        % Clear workspace
        %clear stimval;

    end                   % END of trials >= current_max_val

    % Clear workspace
    %clear stim_env;
    %clear psth_rec;
end		 % END of k
