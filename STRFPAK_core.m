%  STRFPAK_core.m
% 
%        -- 1. Preprocessing Users' input
%        -- 2. Estimation of STRF
%        -- 3. Prediction of neuron response
%        -- 4. Validation of goodness of fit
%
%   Junli Zhang, 1/23/2004
%

% Remember current raw data directory
start_cd = pwd;
preprocessName = {'SongWave->Spectrogram(STFFT)','SongWave->Scalogram(wavelet)'};

%  
% Get all the information of that subdirectory including files and i
% directory by calling dir function
%
DD = dir(pwd);
orig_outputPath = outputPath_temp

% Go through all the files the subdirectory
%     1. Get number of data sets i.e. 20
%     2. Assign my global variable, i.e. RawDS
%
for isubdir = 1:length(DD)
    
    % Go through all the scale options for each subdirectory
    for scaleoption = filteroption
        
        cd(start_cd);
        
        % Check if it is . file or .. file
        if strcmp(DD(isubdir).name, '.')
            continue;
        elseif strcmp(DD(isubdir).name, '..')
            continue;
        else
            % go to subsubdirectory, i.e. conspecfic, or flatrip or songrip
            cd (DD(isubdir).name);
            
            % Remember the current subsubdirectory
            %          
            current_Path = pwd;
            
            ddd = dir(current_Path); 
            numDir = find([ddd.isdir] == 1);
            TolNumDataSet_Estimation = (length(ddd)-length(numDir))/2;
            
            % input where you want to put your intermediate results
            % ========================================================
            % Create corresponding result directory for same bird and cell 
            % under your given output_Path
            %
            [uppath, cellname] = fileparts(start_cd);
            [upuppath, birdname] = fileparts(uppath);
            cd(orig_outputPath)
            if ~exist(fullfile(orig_outputPath,birdname))
                mkdir(birdname) 
            end
            cd(birdname)
            
            if ~exist(fullfile(pwd,cellname))
                mkdir(cellname); 
            end
            cd(cellname);
            
            if ~exist(DD(isubdir).name)
                mkdir(DD(isubdir).name); 
            end
            cd(DD(isubdir).name);
            outputPath_temp = pwd;
            
            
            % Define other global variables:
            %   1. filterwidth - in Hz
            %   2. filterwidth_ms  - in ms 
            %   3. samplerate  - Hz
            %   4. ampsamprate  - Hz
            %   5. respsamprate  - Hz
            %   6. filteroption  - i.e. linear or log 
            clear global filteroption rawDS samprate
            global filteroption  rawDS samprate 
            filteroption = scaleoption;
            
            
            % ========================================================
            %    1. PREPROCESSING STEP
            % ========================================================
            %  In the next for loop, we need assign global variable: rawDS
            %  and do more assignments
            %      1. Assign global variable DS
            %      2. Assign global variables: NBAND, initialFreq, endFre
            clear global DS NBAND
            global DS NBAND initialFreq endFreq sDim preprocessOption 
            global fwidthHz psth_smoothconst wavelength Npoints
            preprocessOption = preprocessName{preprocessoption}; 
            for ii = 1:TolNumDataSet_Estimation
                
                rawDS{ii}.respfiles = fullfile(current_Path,['spike',num2str(ii)]);
                rawDS{ii}.stimfiles = fullfile(current_Path,['stim',num2str(ii),'.wav']);
                [sp, sname] = fileparts(rawDS{ii}.stimfiles);
                [rp,rname] = fileparts(rawDS{ii}.respfiles);
                
                % 
                % PREPROCESSING stimuli files
                %
                [input, fs] = wavread(rawDS{ii}.stimfiles);
                samprate = fs;
                
                % ---------------------------------
                % Spectrogram option choosen
                % ---------------------------------
                if  preprocessoption == 1
                    
                    % Only create outputPath once
                    if ii == 1
                        clear global outputPath
                        global outputPath
                        if filteroption == 1 % linear-linear
                            if not(exist(['S_Output_lin',num2str(fwidthHz)]))
                                mkdir(['S_Output_lin',num2str(fwidthHz)]);
                            end
                            outputPath = fullfile(outputPath_temp, ...
                                ['S_Output_lin',num2str(fwidthHz)]);
                        else
                            if not(exist(['S_Output_log',num2str(fwidthHz)]))
                                mkdir(['S_Output_log',num2str(fwidthHz)]);
                            end
                            outputPath = fullfile(outputPath_temp,...
                                ['S_Output_log',num2str(fwidthHz)]);
                        end
                    end % END of creating outputPath
                    
                    %
                    % 1.1. Spectrogram calculation
                    %
                    [yy, xx] = ComplexSpectrum(input, floor(fs/ampsamprate),...
                        floor(1/(2*pi*fwidthHz)*6*fs),fs);
                    
                    samprate = fs;
                    % Working on more options -JXZ- Aug. 2003
                    if filteroption == 1 % linear-linear 
                        tmp = abs(yy);
                    else  % just temporary
                        tmp = log(abs(yy) +1);
                    end
                    if isempty(initialFreq) & isempty(endFreq)
                        initialFreq = 0;
                        endFreq = fs/2;
                    else
                        if endFreq > fs/2
                            ttt=warndlg('Max frequency limit is samprate/2.', 'high frequency warning','modal');
                            uiwait(ttt);
                            endFreq = fs/2;
                        end
                    end
                    freq_range = find(xx>=initialFreq & xx<=endFreq);
                    
                    %         if rem(size(tmp,1), 2)    % Odd
                    %             cutoff = (size(tmp,1) +1)/2;
                    %         else
                    %             cutoff = size(tmp,1)/2+1;
                    %         end
                    
                    outSpectrum = 10*tmp(freq_range,:);
                    fo = xx(freq_range);
                    NBAND = size(outSpectrum, 1);
                    
                    save(fullfile(outputPath,[sname,'_Stim_',num2str(ii),'.mat']), 'outSpectrum');
                    
                    DS{ii}.stimfiles = fullfile(outputPath,[sname,'_Stim_',num2str(ii),'.mat']);
                    DS{ii}.nlen = round(length(input)*1000/fs) +1;
                    
                    % 
                    % END OF SPECTROGRM  CAlCULATION
                    %
                    
                elseif preprocessoption == 2
                    
                    %
                    % 1.2. Wavelet transformation
                    %
                    % Only create outputPath once
                    if ii == 1
                        clear global outputPath
                        global outputPath
                        if filteroption == 1 % linear-linear
                            if not(exist(['W_Output_lin',num2str(Npoints)]))
                                mkdir(['W_Output_lin',num2str(Npoints)]);
                            end
                            outputPath = fullfile(outputPath_temp,['W_Output_lin',num2str(Npoints)]);
                        else
                            if not(exist(['W_Output_log',num2str(Npoints)]))
                                mkdir(['W_Output_log',num2str(Npoints)]);
                            end
                            outputPath = fullfile(outputPath_temp, ['W_Output_log',num2str(Npoints)]);
                        end
                    end % END of creating outputPath
                    
                    
                    % Wavelet transform
                    [tfr, ttt, fff, wt] = tfrscalo(input, time, wavelength, initialFreq/fs, endFreq/fs,Npoints);
                    % Downsample to new samplerate
                    tfr_new_sample = resample(tfr',ampsamprate, fs);
                    
                    if filteroption == 1 % linear-linear 
                        tmp = (abs(tfr_new_sample'));
                    else  % just temporary
                        tmp = log(abs(tfr_new_sample') +1);
                    end
                    
                    %save(fullfile(outputPath, ['wavelet_org_',num2str(ii),'.mat']),'ttt','fff');
                    save(fullfile(outputPath, [sname,'_Stim_',num2str(ii),'.mat']),'tmp');
                    DS{ii}.stimfiles = fullfile(outputPath,[sname,'_Stim_',num2str(ii),'.mat']);
                    DS{ii}.nlen = round(length(input)*1000/fs) +1;
                    
                    NBAND = size(tmp, 1);
                    
                    
                    %
                    % END OF WAVELET TRANSFORMATION
                    %
                end
                
                %
                %   RESPONSE DATA
                %
                
                [rawResp, trials] = read_spikeTime_2cell(rawDS{ii}.respfiles, round(length(input)*1000/fs) +1);
                
                %rawResp = read_spikeTime(rawDS{ii}.respfiles, size(outSpectrum, 2));
                
                
                % save to the file for each data pair
                save(fullfile(outputPath,[rname,'_Spike_time_',num2str(ii),'.mat']), 'rawResp');
                
                % Assign values to global variable DS
                DS{ii}.respfiles = fullfile(outputPath,[rname,'_Spike_time_',num2str(ii),'.mat']);
                DS{ii}.ntrials = trials;
                
                
            end   % END of TotalNum_Estimation
            
            % ========================================================
            %    2. ESTIMATION STEP
            % ========================================================
            % Assign the rest global variables
           
            
            
            % Start to run calculation of strf.
            disp('Start to run calculation of strf.')
            % call cal_avg to calculate avg of stim
            [stim_avg, avg_psth, psth, errFlg] = cal_AVG(DS);
            
            % Global Variable: SetSep -- used for choosing algorithm.
            clear global SetSep
            global SetSep 
            SetSep = 0;
            [CS, errFlg] = cal_AutoCorr(0,DS, stim_avg, [-TimeLag TimeLag], NBAND);
            [CSR, CSR_JN, errFlg]= cal_CrossCorr(DS,stim_avg,avg_psth,...
                psth,[-TimeLag TimeLag],NBAND);
            calStrf_script;
            % Save global variable for graphical version of STRFPAK
            save(fullfile(outputPath,'GlobalVariables.mat'),'rawDS','DS','NBAND','sDim','SetSep','Std_val',...
                'fwidthHz', 'ampsamprate','TimeLag', 'Tol_val','preprocessOption','wavelength',...
                'initialFreq', 'endFreq','psth_smoothconst','Npoints','respsamprate',...
                'smooth_rt','timevary_PSTH')
            
            
            %Done calculation Stage.
            disp('Done Calculation Stage.')
            
            %To free some memory, recommend you clear CS, CSR, CSRJN.
            clear CS CSR CSRJN;
            
            % ========================================================
            %    3. PREDICTION STEP
            % ========================================================
            
            %Variables for prediction
            clear global predinitialFreq predendFreq matchflg predDS
            global predinitialFreq predendFreq
            predinitialFreq = initialFreq;
            predendFreq = endFreq;
            
            global matchflg predDS
            matchflg = 1;
            
            % Since matchflg is used, we use the same number of data sets as estimation.
            %
            TolNumDataSet_Prediction = TolNumDataSet_Estimation;
            predDS = DS;
            
            twindow = [-TimeLag TimeLag];
            disp('Now doing prediction.')
            save(fullfile(outputPath, 'predVariables.mat'), 'predDS',...
                'predinitialFreq', 'predendFreq');
            
            %Call cal_PredStrf to do prediciton
            for ntol=1:length(Tol_val)
                errFlg = cal_PredStrf(0,predDS, stim_avg, fullfile(outputPath,...
                    ['strfResult_Tol',num2str(ntol),'.mat']),...
                    ntol, twindow, NBAND);
            end
            
            %Done Prediction Stage.
            disp('Done Prediction Stage.')
            
            
            % ========================================================
            %    4. VALIDATION STEP
            % ========================================================
            
            % Validation Variables
            binWindow = 128;
            
            disp('Now doing validation.');
            % Loop through all the results for all the tol values
            % JXZ 7/13/2005
            %   Add extra loop for Std_val list
            
            % Load part I of PSTH of predict response files 
            spredresult = fullfile(outputPath, 'predResult_avgSpike1.mat');
            avgSpike1 = Check_And_Load(spredresult);
            
            % Load part II of PSTH of predict response files 
            spredresult = fullfile(outputPath, 'predResult_avgSpike2.mat');
            avgSpike2 = Check_And_Load(spredresult);
            
            spike_psth = Check_And_Load(fullfile(outputPath,'spike_psth_count.mat'));
            
            ntrials_proper=spike_psth(end);
            spike_psth=spike_psth(1:end-1);
            ntrials_index=find(spike_psth>=ntrials_proper);
            ntrials = spike_psth;
            clear spike_psth;
            nrec = length(ntrials_index);
            % Load all predicted results and form them to a big matrix
            spike_est1 = [];
            spike_est2 = [];
            for k=1:nrec 
                
                %irec = ntrials_index(k);
                % Load avg_spike1.mat
                spike_est1temp = avgSpike1{k};
                
                % Load avg_spike2
                spike_est2temp = avgSpike2{k}; 
                
                spike_est1 = [spike_est1; spike_est1temp];
                spike_est2 = [spike_est2; spike_est2temp];
            end
            
            for istd = 1:length(Std_val)
                for ntol = 1:length(Tol_val)
                    % calculate Info and CC
                    [infopre{ntol}{istd}, infouppre{ntol}{istd}, infodownpre{ntol}{istd}, info{ntol}{istd},...
                            infoup{ntol}{istd}, infodown{ntol}{istd},...
                            cc_spike_pre_best{ntol}{istd}, cc_spike_pre_constant{ntol}{istd},...
                            cc_two_halves_tmax{ntol}{istd},...
                            cc_two_halves_constant{ntol}{istd}, tmax_pre{ntol}{istd},...
                            cc_ratio_max{ntol}{istd}, cc_ratio_const{ntol}{istd},tmax_ratio{ntol}{istd},...
                            fpxypre{ntol}{istd}, cxy{ntol}{istd}, cxypre{ntol}{istd}, cxyup{ntol}{istd}, cxydown{ntol}{istd},...
                            cxyuppre{ntol}{istd}, cxydownpre{ntol}{istd},...
                            cc_spike_pre_corrval{ntol}{istd}, cc_spike_pre_corrval_low{ntol}{istd}, cc_spike_pre_corrval_high{ntol}{istd},...
                            cc_two_halves_corrval{ntol}{istd},cc_two_halves_corrval_low{ntol}{istd}, cc_two_halves_corrval_high{ntol}{istd},...
                            raw_r{ntol}{istd}]...
                        = do_cached_calc('cal_Validate',spike_est1, spike_est2,ntol, istd, ampsamprate, binWindow);
                    %cal_Validate(spike_est1, spike_est2,ntol, istd, ampsamprate, binWindow); %ampsamprate);
                    
                end % END of ntol
            end  % END of istd
            
            
            % For displaying, save them to display files
            save(fullfile(outputPath,'info_r_result.mat'), 'Tol_val', 'infopre', 'infouppre','smoothVect',...
                'infodownpre', 'info', 'infoup', 'infodown',...
                'cc_spike_pre_best', 'cc_spike_pre_constant','cc_ratio_const',...
                'cc_two_halves_tmax', 'cc_two_halves_constant', 'tmax_pre',...
                'cc_ratio_max', 'tmax_ratio');
            
            save(fullfile(outputPath,'display_CC_result.mat'),'cc_spike_pre_corrval','cc_spike_pre_corrval_low','cc_spike_pre_corrval_high','raw_r',...
                'cc_two_halves_corrval','cc_two_halves_corrval_low','cc_two_halves_corrval_high','cc_ratio_max', 'tmax_ratio',...
                'cc_ratio_const');
            
            save(fullfile(outputPath,'display_INFO_result.mat'),'fpxypre','cxy','cxyup','cxydown',...
                'cxypre', 'cxyuppre', 'cxydownpre', 'infopre', 'infouppre',...
                'infodownpre', 'info', 'infoup', 'infodown', 'raw_r');
    
            disp('Done Validation Stage.')
            
            %You can call any related display functions to display them.
            %Or run "strfFirstGUI" and click "Display xxx" Button.
            
            disp(['Done analyzing of ',current_Path])
            disp(['All results saved to ',outputPath])
        end % end of else section
    end % end of scaleoption
end  % end of isubdir 

% END of _script
