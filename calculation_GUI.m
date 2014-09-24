function varargout = calculation_GUI(varargin)
% CALCULATION_GUI Application M-file for songwave_specgram.fig
%    FIG = CALCULATION_GUI launch songwave_specgram GUI.
%    CALCULATION_GUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 26-Apr-2006 14:37:34
if nargin == 0  % LAUNCH GUI

    fig = openfig(mfilename);

    % for resize property
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
        hUIControls],'units','normalized','fontunits','normalized');

    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);

    % display if TimeLag and Tol_val are already set
    global TimeLag TimeLagUnit
    if ~isempty(TimeLag)
        set(handles.twindow, 'String', TimeLag);
    end

    if ~isempty(TimeLagUnit)
        if strcmp(TimeLagUnit, 'msec')
            set(handles.timelagUnit, 'value', 2);
        elseif strcmp(TimeLagUnit, 'frame')
            set(handles.timelagUnit, 'value',3);
        end
    end
    global Tol_val
    if ~isempty(Tol_val)
        tolval_temp = Tol_val*100;
        set(handles.tol_startvalue, 'String', num2str(tolval_temp, 3));
    end

    % JXZ: 7/12/2005
    % Redisplay if Std_val already assigned.
    global Std_val
    if ~isempty(Std_val)
        set(handles.stdval_list, 'String', num2str(Std_val, 3));
    end

    global timevary_PSTH
    if ~isempty(timevary_PSTH)
        if timevary_PSTH == 1
            set(handles.yesRadio, 'Value', 1);
        else
            set(handles.noRadio, 'Value', 1);
        end
    end

    global setSep
    if ~isempty(setSep)
        if setSep == 1
            set(handles.separable, 'Value', 1);
        else
            set(handles.nonseparable, 'Value', 1);
        end
    end

    % Add smoothing window size for time-varying mean rate
    global smooth_rt
    if ~isempty(smooth_rt)
        set(handles.smoothmeanrate, 'String', smooth_rt);
    end

    handles.tolindex = 1;
    handles.stdindex = 1;
    handles.imageIndex = 0;
    handles.numList = 20;  % subaxes for 2-D strf
    handles.outSpectrum = {};

    % total numbers of image on the screen
    for ii = 1:handles.numList
        axesStr = strcat('axes', num2str(ii));
        subhandle = getfield(handles, axesStr);
        set(subhandle, 'Visible', 'off');
    end
    set(handles.strf, 'visible', 'off');
    %set(handles.tolstd, 'visible', 'off');

    guidata(fig, handles);

    if nargout > 0
        varargout{1} = fig;
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    %try
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    %catch
    %   disp(lasterr);
    %end

end

% --------------------------------------------------------------------
function varargout = ampsamprate_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global ampsamprate fwidthHz
if ~isempty(fwidthHz)
    fwhz = fwidthHz;
else
    fwhz = 250;
end
newStr = get(h, 'String');
NewVal = str2double(newStr);

% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid frequency bandwith (positive number only).',...
        'Variable Error', 'modal')
    return;
elseif NewVal < 2*pi*sqrt(2)*fwhz
    ttt=warndlg('You are undersampling your data.', 'amp_samp_rate warning','modal');
    uiwait(ttt);
end

ampsamprate = NewVal;

% --------------------------------------------------------------------
function varargout = respsamprate_Callback(h, eventdata, handles, varargin)

global respsamprate

newStr = get(h, 'String');
NewVal = str2double(newStr);

% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid frequency bandwith (positive number only).',...
        'Variable Error', 'modal')
    return;

end

respsamprate = NewVal;


% --------------------------------------------------------------------
function varargout = smoothconst_Callback(h, eventdata, handles, varargin)
global psth_smooth

newStr = get(h, 'String');
NewVal = str2double(newStr);

% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter postive number.',...
        'Variable Error', 'modal')
    return;

end

psth_smooth = NewVal;
plot_strfonly(handles);

% --------------------------------------------------------------------
function varargout = compsave_Callback(h, eventdata, handles, varargin)
%  Ask where you want to put your intermediate results
%
global DS outputPath OnlyOne
if isempty(DS)
    errordlg('Please select input data first.',...
        'Data Input Error', 'modal')
    return
else
    if length(DS) ==1  % Whether want to do cross validation
        if isempty(outputPath)

            currentPath = pwd;
            prompt={['Please Enter where you want to put your intermediate results:']};
            def = {fullfile(currentPath, 'Output')};
            dlgTitle='Path for intermediate results';
            lineNo=1;

            % picture feature
            AddOpts.Resize='on';
            AddOpts.WindowStyle='normal';
            AddOpts.Interpreter='tex';
            datadir =inputdlg2(prompt,dlgTitle,lineNo,def,AddOpts);

            % Check if user input valid directory
            if isempty(datadir)
                errordlg('Please enter valid directory.','Input Error','modal')
                return
            end

            outdatafiledir = datadir{1};

            if not(exist(outdatafiledir,'dir'))
                disp('Directory not found. Creating new directory.');
                [p, n, e] = fileparts(outdatafiledir);
                if not(exist(p, 'dir'))
                    errordlg('The upper directory is not found and so exists. ');
                    return
                end
                cd (p)
                mkdir(n)
                % Junli: 11/11/2005
            end

            outputPath = outdatafiledir;
        end

        if isempty(OnlyOne)
            answ = questdlg({'NOTE: This message shows only if you choose ONE data set.',...
                'There are two options for the data division:             ',...
                '    Use the whole data set:                              ',...
                '        The whole data set is used for estimation only.  ',...
                '        Later the new data sets are needed for validation',...
                '                                                         ',...
                '    Split 90%/10% data division:                         ',...
                '        The 90% of this data is used for estimation      ',...
                '        and the 10% of the data is used for validation.  ',...
                '        The splitted data are saved to the output        ',...
                '        directory:', outputPath,...
                '                                                         ',...
                'What is your choice?                                     ',...
                '                                                      '},...
                'Warning Message', 'Use the whole data set',...
                'Split 90%/10% data set', 'Use the whole data set');
            switch answ
                case 'Use the whole data set'
                    clear global predDS
                    OnlyOne = 1;

                case 'Split 90%/10% data set'

                    % OnlyOne flag used for avoiding overfitting
                    clear global predDS
                    global predDS
                    OnlyOne = 1;

                    % Divide the orignal data files into two parts

                    stim_env = Check_And_Load(DS{1}.stimfiles);
                    psth_rec = Check_And_Load(DS{1}.respfiles);
                    oldLength = min(length(stim_env), length(psth_rec));

                    partI = stim_env(:, 1:floor(0.9*oldLength));
                    partII = stim_env(:, floor(0.9*oldLength)+1:oldLength);
                    [p, sn, se] = fileparts(DS{1}.stimfiles);
                    save(fullfile(outputPath,['partI_' sn,'.mat']), 'partI');
                    save(fullfile(outputPath,['partII_' sn, '.mat']), 'partII');

                    DS{1}.stimfiles = fullfile(p, ['partI_', sn, '.mat']);
                    predDS{1}.stimfiles = fullfile(p, ['partII_', sn, '.mat']);

                    DS{1}.nlen = floor(0.9 * oldLength);
                    predDS{1}.nlen = oldLength - DS{1}.nlen;
                    predDS{1}.ntrials = DS{1}.ntrials;

                    % Divide orignal response files

                    [p, rn, re] = fileparts(DS{1}.respfiles);
                    partI = psth_rec(:, 1:floor(0.9*oldLength));
                    partII = psth_rec(:, floor(0.9*oldLength)+1:oldLength);

                    save(fullfile(outputPath,['partI_', rn, '.mat']), 'partI');
                    save(fullfile(outputPath,['partII_', rn, '.mat']), 'partII');

                    DS{1}.respfiles = fullfile(p, ['partI_', rn, '.mat']);
                    predDS{1}.respfiles = fullfile(p, ['partII_', rn, '.mat']);
                otherwise
                    errordlg('Not valid selection', 'modal')
                    return;

            end
        end
    end
    %clear global OnlyOne
end

%%%$$$ begin calculate

% Begin to calculate
global TimeLag Tol_val setSep TimeLagUnit timevary_PSTH smooth_rt
global Std_val preprocessOption
if isempty(TimeLagUnit)
    tt = errordlg('Please specify Time-lag unit (msec or frames).',...
        'Time-lag Unit missing','modal');
    uiwait(tt);
    return;
end
if isempty(TimeLag) | isempty(Tol_val) ...
        | isempty(setSep) | isempty(timevary_PSTH)...
        | isempty(smooth_rt)
    tt= errordlg('Please set all calculation parameter values.',...
        'Parameter Value Missing', 'modal');
    uiwait(tt);
    return;
end
global running_in_script_mode
if ~strcmp(running_in_script_mode,'yes')

    % show the current path
    set(handles.outdirshow, 'string', outputPath);
    set(handles.outdirshow, 'visible', 'on');
    set(handles.outdirBrowser, 'visible', 'on');


    set(handles.figure1, 'Pointer', 'Watch');
end
% =========================================================
% calculate avg. of stimulus and response that used later
% =========================================================

[stim_avg, avg_psth, psth, errFlg] = cal_AVG(DS);

% Check if cal_Avg ends normally
if errFlg == 1
    set(handles.figure1, 'Pointer', 'Arrow');
    tt = errordlg('Cal_AVG has error', 'cal_AVG Error', 'modal');
    uiwait(tt);
    return
end

% =========================================================
% Now calcualting stimulus AutoCorr.
% =========================================================

global NBAND sDim
% Nov 22, 2005
% TimeLag need has ms unit, so I define twindow.
global ampsamprate
if isempty(ampsamprate)
    ampsamprate = 1000;
end
if strcmp(TimeLagUnit, 'msec')

    twindow = [-round(TimeLag*ampsamprate/1000) round(TimeLag*ampsamprate/1000)];
elseif strcmp(TimeLagUnit, 'frame')
    twindow = [-TimeLag TimeLag];
end

if setSep == 0  % Nonseparable space-time algorithm
    disp('Now calculating stim autocorrelation');
    do_long_way = 1;
    [cached_dir,maxsize] = dir_of_caches;
    autocorr_start_time = cputime;
    hashes_of_stims = create_stim_cache_file(outputPath,DS);


    if ~strcmp(cached_dir,'null')
        do_long_way = 0;
        [loaded,order] = sort(hashes_of_stims);  %  Sort to make the checksum invarient to dataset shuffling
        n_trial_array = get_ntrials(DS);
        checksum_for_autocorr_calc = checksum(load_function_text('cal_AutoCorr'),loaded,n_trial_array(order),stim_avg,twindow,NBAND);  % Sort the ntrial array the same way as the stimuli
        [CS,errFlg] = do_cached_calc_checksum_known('cal_AutoCorr',checksum_for_autocorr_calc,1, DS, stim_avg, twindow, NBAND);
    else
        [CS, errFlg] = cal_AutoCorr(1, DS, stim_avg, twindow, NBAND);
    end

    %[DS_data,the_checksum] = load_DS_data(DS,stim_avg,twindow, NBAND);
    %autocorr_start_time = cputime;
    %[CS, errFlg] = do_cached_calc_checksum_known('cal_AutoCorr_for_cache',the_checksum, DS_data, stim_avg,twindow, NBAND);
    autocorr_end_time = cputime;
    disp(['The autocorrelation took ' num2str(autocorr_end_time - autocorr_start_time) ' seconds.']);
    currentPath = pwd;
    global outputPath
    if ~isempty(outputPath)
        cd (outputPath);
    else
        disp('Saving output to Output Dir.');
        stat = mkdir('Output');
        cd('Output');
        outputPath = pwd;
    end

    save('Stim_autocorr.mat', 'CS');
    cd(currentPath);


    % Check if cal_AutoCorr ends normally
    if errFlg == 1
        set(handles.figure1, 'Pointer', 'Arrow');
        return
    end

    % Done calcualtion of stimulus AutoCorr
    % =========================================================
    % Now calcualting stimulus spike Cross Corr.
    % =========================================================

    %  Let's assume that if the user has caching on and is using the GUI
    %  that they might want to evaluate cal_crossCorr more than once; so:
    global running_in_script_mode
    cache_crosscorr = ~strcmp(cached_dir,'null');  % This is the flag to say if we'll cache results specific to the current spike train.

    hashes_of_stims = create_stim_cache_file(outputPath,DS);
    hashes_of_spikes = create_spike_cache_file(outputPath,DS);
    global smooth_rt
    if isempty(smooth_rt)
        smooth_rt = 41;
    end
    if ~exist('psth_option')
        global timevary_PSTH
        if timevary_PSTH == 0
            psth_option = 0;
        else
            psth_option = 1;
        end
    end
    checksum_CrossCorr = checksum(load_function_text('cal_CrossCorr'),hashes_of_spikes,hashes_of_stims,twindow,smooth_rt,psth_option);

    %if cache_crosscorr & ~strcmp(running_in_script_mode,'yes')  % As of version 5.2, use locally-cached folder always
    if ~strcmp(cached_dir,'null')
        %[CSR, CSR_JN, errFlg]= do_cached_calc_checksum_known('cal_CrossCorr',checksum_CrossCorr,DS,stim_avg,avg_psth,psth,twindow,NBAND);
        [CSR, CSR_JN, errFlg]= do_locally_cached_calc_checksum_known(get_local_cache_dir,'cal_CrossCorr',checksum_CrossCorr,DS,stim_avg,avg_psth,psth,twindow,NBAND);
        save(fullfile(outputPath,'StimResp_crosscorr.mat'), 'CSR');
        save(fullfile(outputPath,'SR_crosscorrJN.mat'), 'CSR_JN');


    else

        [CSR, CSR_JN, errFlg]= cal_CrossCorr(DS,stim_avg,avg_psth,...
            psth,twindow,NBAND);
    end

    % Check if cal_CrossCorr ends normally
    if errFlg == 1
        set(handles.figure1, 'Pointer', 'Arrow');
        return
    end
    % =========================================================
    %Done calcualtion of stimulus-spike CrossCorr in GUI window
    % =========================================================
    disp('Calculating strfs for each tol value.');
    %checksum_CrossCorr
    calStrf_script;
    calculation_endtime = cputime;
    disp(['The STRF calculation took ' num2str(calculation_endtime - autocorr_start_time) ' seconds.']);

else % Separable space-time algorithm

    % Provide Space-time separability algorithm to estimate STRF

    %         [CSspace, CStime, errFlg] = cal_AutoCorrSep(DS, stim_avg,...
    %         twindow, NBAND, 1);


    disp('Now calculating stim autocorrelation');
    do_long_way = 1;
    [cached_dir,maxsize] = dir_of_caches;
    autocorr_start_time = cputime;
    hashes_of_stims = create_stim_cache_file(outputPath,DS);


    if ~strcmp(cached_dir,'null')
        do_long_way = 0;
        [loaded,order] = sort(hashes_of_stims);  %  Sort to make the checksum invarient to dataset shuffling
        n_trial_array = get_ntrials(DS);
        checksum_for_autocorr_calc = checksum(load_function_text('cal_AutoCorrSep'),'cal_AutoCorrSep',loaded,n_trial_array(order),stim_avg,twindow,NBAND);  % Sort the ntrial array the same way as the stimuli
        [CSspace, CStime, errFlg] = do_cached_calc_checksum_known('cal_AutoCorrSep',checksum_for_autocorr_calc, DS, stim_avg, twindow, NBAND,1);
    else
        [CSspace, CStime, errFlg] = cal_AutoCorrSep(DS, stim_avg, twindow, NBAND,1);
    end

    %[DS_data,the_checksum] = load_DS_data(DS,stim_avg,twindow, NBAND);
    %autocorr_start_time = cputime;
    %[CS, errFlg] = do_cached_calc_checksum_known('cal_AutoCorr_for_cache',the_checksum, DS_data, stim_avg,twindow, NBAND);
    autocorr_end_time = cputime;
    disp(['The autocorrelation took ' num2str(autocorr_end_time - autocorr_start_time) ' seconds.']);

    % Check if cal_AutoCorrSep ends normally
    if errFlg == 1
        set(handles.figure1, 'Pointer', 'Arrow');
        return
    end

    % Calculate cross-correlation between stimuli and spike


    %  Let's assume that if the user has caching on and is using the GUI
    %  that they might want to evaluate cal_crossCorr more than once; so:

    cache_crosscorr = ~strcmp(cached_dir,'null');  % This is the flag to say if we'll cache results specific to the current spike train.
    if ~exist('psth_option')
        global timevary_PSTH
        if timevary_PSTH == 0
            psth_option = 0;
        else
            psth_option = 1;
        end
    end

    hashes_of_stims = create_stim_cache_file(outputPath,DS);
    hashes_of_spikes = create_spike_cache_file(outputPath,DS);
    checksum_CrossCorr = checksum('sep',load_function_text('cal_CrossCorr'),hashes_of_spikes,hashes_of_stims,twindow,smooth_rt,psth_option);
    global smooth_rt
    if isempty(smooth_rt)
        smooth_rt = 41;
    end
    global running_in_script_mode
    if cache_crosscorr 


        [CSR, CSR_JN, errFlg]= do_locally_cached_calc_checksum_known(get_local_cache_dir,'cal_CrossCorr',checksum_CrossCorr,DS,stim_avg,avg_psth,psth,twindow,NBAND);
        save(fullfile(outputPath,'StimResp_crosscorr.mat'), 'CSR');
        save(fullfile(outputPath,'SR_crosscorrJN.mat'), 'CSR_JN');


    else

        [CSR, CSR_JN, errFlg]= cal_CrossCorr(DS,stim_avg,avg_psth,...
            psth,twindow,NBAND);
    end
    %     [CSR, CSR_JN, errFlg]= cal_CrossCorr(DS,stim_avg,avg_psth,...
    %         psth,twindow,NBAND);

    % Check if cal_CrossCorr ends normally
    if errFlg == 1
        set(handles.figure1, 'Pointer', 'Arrow');
        return
    end

    % Now call calStrfSep_script to calculate STRF, STRF_JN
    % STRF_JNstd for each tol value
    calStrfSep_script;

end

% Now saving the results
% Save global variable for graphical version of STRFPAK
global rawData rawDS ampsamprate stimsamprate respsamprate
global initialFreq endFreq OnlyOne
save(fullfile(outputPath,'Global_Variables.mat'),'DS','NBAND',...
    'sDim','setSep','Std_val','rawData', 'rawDS','ampsamprate','stimsamprate',...
    'TimeLag', 'Tol_val','preprocessOption','respsamprate',...
    'smooth_rt','timevary_PSTH', 'initialFreq', 'endFreq',...
    'OnlyOne')
%%%$$$ end calculate

% major calculation done, then the screen becomes normal.
set(handles.figure1, 'Pointer', 'Arrow');
guidata(handles.figure1, handles);

% =========================================================
% Give the sign for finishing this stage
% =========================================================
disp('Done calculation Stage.');
display_Callback(h, eventdata, handles, varargin);
msgbox('Done STRF Calculation.', 'modal');
% --------------------------------------------------------------------
function varargout = display_Callback(h, eventdata, handles, varargin)

global Std_val Tol_val
tollen = length(Tol_val);
stdlen = length(Std_val);

if tollen ==1
    set(handles.tolvalSlider, 'min', 0, 'max', tollen, 'value',0,'sliderstep', [1 1]);
elseif tollen > 1
    set(handles.tolvalSlider, 'Min', 1, 'Max', tollen,'value',1,...
        'sliderstep', [1/(tollen-1) 1/(tollen-1)]);
end
if stdlen ==1
    set(handles.stdvalSlider, 'min', 0, 'max', stdlen, 'value',0,'sliderstep', [1 1]);
elseif stdlen > 1
    set(handles.stdvalSlider, 'Min', 1, 'Max', stdlen,'value',1,...
        'sliderstep', [1/(stdlen-1) 1/(stdlen-1)]);
end
guidata(h, handles);

plot_strfonly(handles);

% --------------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
global DS outputPath OnlyOne TimeLag Tol_val setSep TimeLagUnit timevary_PSTH smooth_rt
global Std_val preprocessOption NBAND sDim initialFreq endFreq
% Nov 22, 2005
% TimeLag need has ms unit, so I define twindow.
global ampsamprate
if isempty(ampsamprate)
    ampsamprate = 1000;
end
save(fullfile(outputPath,'STRFPAK_script_parameters.mat'), 'OnlyOne' , 'TimeLag' , 'Tol_val' , 'setSep' , 'TimeLagUnit' , 'timevary_PSTH' , 'smooth_rt', ...
    'Std_val' , 'preprocessOption' , 'NBAND' , 'sDim' , 'initialFreq' , 'endFreq' , '-APPEND');
add_to_STRFPAK_script('calculation_GUI.m','calculate');

delete(handles.figure1);

% % --------------------------------------------------------------------
% function varargout = filteroption_Callback(h, eventdata, handles, varargin)
%     v = get(handles.filteroption, 'value');
%     option = get(handles.filteroption, 'String');
%     option = deblank(option(v,:));
%
%     global filteroption
%
%     if strcmp(option, 'Linear')
%         filteroption = 1;
%     elseif strcmp(option, 'Logarithmic')
%         filteroption = 2;
%     end
%

% --------------------------------------------------------------------
function varargout = reset_Callback(h, eventdata, handles, varargin)

% 5/24/2005: Also delete the saved files
current_dir = pwd;
global outputPath
cd(outputPath)
delete *_Spike_time*;
delete *_Stim_*;
cd(current_dir);

%clear DS initialFreq endFreq StimSampRate sDim
clear global rawDS fstep fwidthHz filteroption outputPath
clear global respsamprate samprate psth_smooth
clear global DS NBAND num_trials endFreq initialFreq ampsamprate sDim

set(handles.fwidthHz, 'String', ' ');
set(handles.fwidthms, 'String', ' ');
set(handles.ampsamprate, 'String', ' ');
set(handles.respsamprate, 'String', ' ');
set(handles.samprate, 'String', ' ');
set(handles.smoothconst, 'String', ' ');
set(handles.upfreq, 'String', ' ');
set(handles.lowfreq, 'String',' ');
set(handles.fstep, 'String', ' ');
set(handles.filename, 'String', ' ');
set(handles.rfilename, 'String', ' ');

%v = get(handles.dimension, 'String');
set(handles.filteroption, 'value',1);
set(handles.strf, 'visible', 'off');
child = get(handles.strf, 'Children');
set(child, 'Visible', 'off')
set(handles.psth, 'visible', 'off');
child = get(handles.psth, 'Children');
set(child, 'Visible', 'off')

set(handles.prev, 'visible', 'off');
set(handles.next, 'visible', 'off');
set(handles.smoothpsth,'visible','off');
set(handles.smoothconst,'visible','off');

% Save all handles to figures
guidata(h,handles);
% --------------------------------------------------------------------
function varargout = samprate_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function plot_strfonly(handles)
% --------------------------------------------------------------------


% Get calc parameter value from global variables
global sDim initialFreq endFreq ampsamprate preprocessOption
global outputPath Std_val Tol_val
set(handles.tolvaltext, 'visible', 'on');
set(handles.stdvaltext, 'visible', 'on');
set(handles.tolvalshow, 'visible','on');
set(handles.stdvalshow, 'visible', 'on');
set(handles.tolvalSlider, 'visible', 'on');
set(handles.stdvalSlider, 'visible', 'on');
handles.tolindex = min(handles.tolindex,length(Tol_val));
handles.stdindex = min(handles.stdindex,length(Std_val));
set(handles.tolvalshow, 'string', Tol_val(handles.tolindex));
set(handles.tolvalSlider,'value', handles.tolindex);
set(handles.stdvalshow, 'string', Std_val(handles.stdindex));
set(handles.stdvalSlider,'value', handles.stdindex);

% show the current path
set(handles.outdirshow, 'string', outputPath);
set(handles.outdirshow, 'visible', 'on');
set(handles.outdirBrowser, 'visible', 'on');

strf_result = load(fullfile(outputPath,...
    ['strfResult_Tol',num2str(handles.tolindex),'.mat']));
forwardJN = strf_result.STRFJN_Cell;
forwardJNstd = strf_result.STRFJNstd_Cell;
[nx, nt, nJN]= size(forwardJN);

% Filter the filters
forwardJN_s = zeros(nx, nt, nJN);
stdfilt = Std_val(handles.stdindex);
if nJN == 1
    forwardJN_s = strf_result.STRF_Cell;
else
    for iJN =1:nJN
        %  find the filtered JN STRF
        forwardJN_s(:,:,iJN) = fast_filter_filter(mean(forwardJN(:,:,[1:(iJN-1) (iJN+1):nJN]),3),...
            forwardJNstd(:,:,iJN), stdfilt);
    end
end
forward = squeeze(mean(forwardJN_s, 3));
maxforward = max(max(forward));
minforward = min(min(forward));
absforward = max(abs(minforward),abs(maxforward));

handles.nt = nt;
guidata(handles.figure1,handles);

t = -(nt-1)/2:(nt-1)/2;

% total numbers of small image on the screen
for ii = 1:handles.numList
    axesStr = strcat('axes', num2str(ii));
    subhandle = getfield(handles, axesStr);
    set(subhandle, 'Visible', 'off');
    child = get(subhandle, 'Children');
    set(child, 'Visible', 'off')

end

% Display estimated STRF based on sptial dimension
if strcmp(sDim, '1-D')
    global initialFreq endFreq

    if isempty(initialFreq) | isempty(endFreq)
        f = 1:nx;
        flabel = logspace(log10(1), log10(nx), nx);
        xtext = 'Frequency Band';
    else
        f = (linspace(initialFreq, endFreq, nx))/1000;
        flabel = logspace(log10(initialFreq), log10(endFreq), nx)/1000;
        xtext = 'Frequency (kHz)';
    end

    % 1-D display
    set(handles.strf, 'visible', 'on')

    set(handles.uppertext, 'visible', 'on')
    set(handles.lowertext, 'visible', 'on');
    set(handles.upperfeq, 'visible', 'on');
    set(handles.lowerfeq, 'visible', 'on');
    set(handles.upperfeq, 'string', endFreq);
    set(handles.lowerfeq, 'string', initialFreq);

    axes(handles.strf);

    if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')

        pcolor(ceil(t*1000/ampsamprate), flabel, forward); shading interp;
        caxis([-absforward absforward]);
        axis([ceil(t(1)*1000/ampsamprate) ceil(t(end)*1000/ampsamprate) flabel(1) flabel(end)])
    else
        imagesc(ceil(t*1000/ampsamprate),f,forward,[-absforward absforward]);
        axis xy;
        axis([ceil(t(1)*1000/ampsamprate) ceil(t(end)*1000/ampsamprate) f(1) f(end)])
    end
    %xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
    xlabel('Time (ms)')
    ylabel(xtext);

    title('STRF');
elseif strcmp(sDim, '2-D')

    % Set other axes to invisible
    set(handles.strf, 'visible', 'off');
    global preprocessOption
    if strcmp(preprocessOption, 'Fourier Power transform')
        phasecount=1;
        tbincount = nt;
        kcount = 1;
        spacebincount = nx;
        chancount=spacebincount/phasecount;
        Xmax=sqrt(chancount*2);
        [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);

        tsf=zeros(Xmax*Xmax,tbincount,kcount);
        tsf(cfilt,:,:)=squeeze(sum(reshape(forward,chancount,phasecount,...
            tbincount,kcount),2));

        tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
        STA=reshape(tsf,Xmax,Xmax,tbincount);
    else
        % 2-D display
        splitX = floor(sqrt(nx));
        STA = reshape(forward, splitX, splitX, nt);

    end
    amax = max(abs(STA(:)));
    amin = -amax;
    numList = nt;
    global stimsamprate
    if isempty(stimsamprate)
        binsize = 1;
        titlestr = 'Latency(frame)';
    else
        binsize = ceil((1/stimsamprate)*1000);
        titlestr = 'Latency(msec)';
    end

    % display total numList images on one screen
    for fr = handles.imageIndex * numList+1:numList+handles.imageIndex * numList-(nt-1)/2
        if fr > nt | fr -handles.imageIndex*numList > handles.numList
            break
        end

        % display them in appropriate axes
        axesStr = strcat('axes', num2str(num2str(fr-handles.imageIndex * numList)));
        axes(getfield(handles, axesStr));

        % Display them in specific format

        if amin~=amax,
            ttH = imagesc(STA(:,:,fr+(nt-1)/2),[amin,amax]);

        else
            %ttH = imagesc(zeros(size(STA,1),size(STA,2)));
            ttH = imagesc(STA(:,:,fr+(nt-1)/2));
        end

        % Display them without label
        if size(STA,1)>1 & size(STA,2)>1,
            axis image
            set(get(ttH,'Parent'),'YTickLabel',[]);
            %set(get(ttH,'Parent'),'YTick',[]);
            set(get(ttH,'Parent'),'XTickLabel',[]);
            %set(get(ttH,'Parent'),'XTick',[]);
            %axis off
        end

        %set(h2,'YTickLabel',[]);
        %set(h2,'XTickLabel',[]);

        if (fr -1 ) ~= 0
            %curTitle = num2str((fr-1-(nt-1)/2)*binsize);
            curTitle = num2str((fr-1)*binsize);
        else
            curTitle = titlestr;
        end
        title(curTitle);
        colormap(redblue);
    end

elseif strcmp(sDim, '0-D')
    axes(handles.strf);
    plot(ceil(t*1000/ampsamprate),forward);
    xlabel('Time (ms)');
    axis tight;

end

% Display appropriate tol-value on the window
guidata(handles.figure1, handles);

% --- Executes on button press in smoothpsth.
% --------------------------------------------------------------------
function smoothpsth_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------

helpdlg({'Smoothing PSTH: ',...
    ' ',...
    ' This flag is used for smoothing psth when displaying psth.',...
    ' It is in ms. You can type new number in the text box to modify it.',...
    ' '},...
    'Smooth PSTH Help');


% --------------------------------------------------------------------
function varargout = help_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
ttlStr = get(handles.figure1, 'Name');
hlpStr = [...
    '   STRFPAK: Calculation Window (HELP WINDOW)                           '
    '                                                                       '
    'The right panel is calculation parameters panel. All parameter buttons '
    ' can be clicked to show a brief explaination                           '
    '                                                                       '
    'These parameters are:                                                  '
    '                                                                       '
    '     Time-lag -- is the time lag used for calculating the spike        '
    'triggered average and the stimulus auto-correlation. Time-lag should   '
    ' be long enough to capture the memory of the neuron and the correlation'
    '  time of the stimulus. The unit of Time-lag is ms.                    '
    '                                                                       '
    ' Tol Value: a regularization hyperparameter used to estimate the       '
    ' STRF. Larger tol values yield smooth STRFs. Lower tol values remove   '
    'more of the stimulus correlation but can amplify noise. Use a range of '
    'them to find the best STRF for your specific stimulus and neuron.      '
    'You can separate them by ",", "space" or "tab". For example:           '
    'Tol_Val = [0.1 0.05 0.001 0.0005].                                     '
    '                                                                       '
    ' Sparse paramter: a second hyperparameter in the STRF estimation that  '
    'enforces sparseness. High Std values limits the number of significant  '
    'pixels in the STRF. Use a range of them to find the best STRF for your '
    'specific stimulus and neuron. You can separate them by ",", "space" or '
    ' "tab". For example: Std_Val = [0 0.1 2 4].                            '
    '                                                                       '
    '   Model Selection: Model #1, r+STRF, means that overall mean firing   '
    'rate is removed from neuronal response when estimating STRF.           '
    ' Model #2, r(t)+STRF, means that time-varying mean firing rate is      '
    'removed from neuron response. When Model #2 is selected the smoothing  '
    'window for the time-varying mean rate needs to be specified. The       '
    'smoothing window is in ms.  This model is useful when short stimuli are'
    'used and the neural response to all stimuli has a stong an onset       '
    'componnent. The predictions of the models to new stimuli will include  '
    'both the mean rate (constant or time-varying) and the STRF prediction. '
    '                                                                       '
    '                                                                       '
    '  Compute and Save: Compute the spectrogram of the signal and save     '
    '  the results into the directory where you will be asked to input.     '
    ' The computing status bar also shows up so that you can know progress. '
    '  Display: Graphically display the spectrogram of the stimulus and     '
    '  the smoothed psth with  Smooth_PSTH window size. The smoothing width '
    ' for psth is shown here. You can modify it by typing different number. '
    '  Reset: Reset all the parameters and the data sets chosen.            '
    '  Close: Close window and save all the parameters and all the results. '
    '                                                                       '
    ' FOR DETAILED INFORMATION, PLEASE REFER TO STRFPAK DOCUMENTATION.      '
    '                                                                       '
    ' Updated by Junli, May 2006.                                           '
    ];
myFig = handles.figure1;
helpwin(hlpStr, ttlStr);


% --- Executes during object creation, after setting all properties.
function upfreq_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');



function upfreq_Callback(hObject, eventdata, handles)
newStr_respsamprate = get(hObject, 'String');
NewVal = str2double(newStr_respsamprate);

% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid high freq (positive number only).',...
        'Variable Error', 'modal')
    return;
end

% Assign global variable to ampsamprate
global endFreq;
if isempty(endFreq)
    endFreq = NewVal;
else
    endFreq = NewVal;
    compsave_Callback(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function lowfreq_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');




function lowfreq_Callback(hObject, eventdata, handles)

newStr_respsamprate = get(hObject, 'String');
NewVal = str2double(newStr_respsamprate);

% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid low frequency (positive number only).',...
        'Variable Error', 'modal')
    return;
end

% Assign global variable to ampsamprate
global initialFreq;
if isempty(initialFreq)
    initialFreq = NewVal;
else
    initialFreq = NewVal;
    compsave_Callback(hObject, eventdata, handles);
end

function twindow_Callback(hObject, eventdata, handles)

% Get twindow values from fig
newStrValue = get(hObject, 'String');
NewVal = str2double(newStrValue);

% Check if the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    errordlg('Please enter valid Time-lag value.(positive number only)',...
        'Global Variable Error', 'modal');
    return;
end

handles.settwindow = 1;
guidata(hObject, handles);

% Assign the global variable TimeLag
global TimeLag
TimeLag = NewVal;

% --- Executes during object creation, after setting all properties.
function twindow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in twindowhelp.
function twindowhelp_Callback(hObject, eventdata, handles)
helpdlg({'Time-lag: ',...
    ' ',...
    '-- is the time lag used for calculating the spike triggered',...
    ' average and the stimulus auto-correlation.  Time-lag should',...
    ' be long enough to capture the memory of the neuron and the', ...
    ' correlation time of the stimulus. The unit of Time-lag is ms.',...
    ' '},...
    'Time-lag Help');

% --- Executes on button press in yesRadio.
function yesRadio_Callback(hObject, eventdata, handles)
global timevary_PSTH
newVal = get(hObject, 'Value');
if newVal == 1
    timevary_PSTH = 1;
    set(handles.noRadio, 'Value', 0);
end


function tol_startvalue_Callback(hObject, eventdata, handles)
InputStr = get(hObject, 'String');

if isempty(InputStr)
    errordlg('Please enter nonempty tol value list', 'Input Error', 'modal')
    return;
end

global Tol_val
% parse the input string to a set of data files
% delimited by ' '
[stoken, rem] = strtok(InputStr, ' \t,');
if isempty(stoken)
    errordlg('Please enter valid Tol_Val.', 'Input Error', 'modal')
    return;
end
tolval = str2double(stoken);
if isnan(tolval) | tolval < 0
    errordlg('Please enter valid Tol_Val.', 'Input Error', 'modal')
    return;
end

Tol_val = [tolval/100];

while  ~isempty(rem)
    [stoken, rem] = strtok(rem, ' \t,');
    if isempty(stoken)
        break
    end
    tolval = str2double(stoken);
    if isnan(tolval) | tolval < 0
        errordlg('Please enter valid Tol_Val.', 'Input Error', 'modal')
        return;
    end
    Tol_val = [Tol_val tolval/100];
end
Tol_val = unique(Tol_val);
Tol_val = -sort(-Tol_val);

handles.setstart = 1;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tol_startvalue_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tolvalhelp.
function tolvalhelp_Callback(hObject, eventdata, handles)
helpdlg({'Tol. Values: ',...
    ' ',...
    ' -- a regularization hyperparameter used to estimate the STRF.',...
    ' Larger tol values yield smooth STRFs. Lower tol values remove',...
    ' more of the stimulus correlation but can amplify noise.', ...
    ' Use a range of them to find the best STRF for your specific',...
    ' stimulus and neuron. You can separate them by ",", "space"',...
    ' or "tab". For example: Tol_Val = [0.1 0.05 0.001 0.0005].',...
    ' '}, 'Tol Value Help');

% --- Executes on button press in noRadio.
function noRadio_Callback(hObject, eventdata, handles)
global timevary_PSTH
newVal = get(hObject, 'Value');
if newVal == 1
    timevary_PSTH = 0;
    set(handles.yesRadio, 'Value', 0);
end

% --- Executes on button press in stdvalhelp.
function stdvalhelp_Callback(hObject, eventdata, handles)
helpdlg({'Std Values: ',...
    ' ',...
    ' -- a second hyperparameter in the STRF estimation that enforces',...
    ' sparseness. High Std values limits the number of significant',...
    ' pixels in the STRF. Use a range of them to find the best STRF',...
    ' for your specific stimulus and neuron. You can separate them',...
    ' by ",", "space" or "tab". For example: ',...
    ' Std_Val = [0 0.1 2 4].',...
    ' '}, 'Tol Value Help');

function stdval_list_Callback(hObject, eventdata, handles)
InputStr = get(hObject, 'String');

if isempty(InputStr)
    errordlg('Please enter nonempty std value list', 'Input Error', 'modal')
    return;
end

global Std_val DS
% parse the input string to a set of data files
% delimited by ' '
[stoken, rem] = strtok(InputStr, ' \t,');
if isempty(stoken)
    errordlg('Please enter valid Std_Val.', 'Input Error', 'modal')
    return;
end
tolval = str2double(stoken);
if isnan(tolval) | tolval < 0
    errordlg('Please enter valid Std_Val.', 'Input Error', 'modal')
    return;
end

Std_val = [tolval];

while  ~isempty(rem)
    [stoken, rem] = strtok(rem, ' \t,');
    if isempty(stoken)
        break
    end
    tolval = str2double(stoken);
    if isnan(tolval) | tolval < 0
        errordlg('Please enter valid Std_Val.', 'Input Error', 'modal')
        return;
    end
    Std_val = [Std_val tolval];
end
Std_val = unique(Std_val);
Std_val = sort(Std_val);
if length(DS) == 1
    Std_val = [0];
    set(handles.stdval_list, 'String', num2str(Std_val, 3));
end
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function stdval_list_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in modelhelp.
function modelhelp_Callback(hObject, eventdata, handles)
helpdlg({'Model Selection: ',...
    ' ',...
    ' Model #1, r+STRF, means that overall mean firing rate is removed',...
    ' from neuronal response when estimating STRF.',...
    ' Model #2, r(t)+STRF, means that time-varying mean firing rate ',...
    ' is removed from neuron response. When Model #2 is selected',...
    ' the smoothing window for the time-varying mean rate needs to',...
    ' be specified. The smoothing window is in ms.  This model is useful',...
    ' when short stimuli are used and the neural response to all stimuli',...
    ' has a stong an onset componnent.',...
    ' The predictions of the models to new stimuli will include both the mean rate (constant or time-varying)',...
    ' and the STRF prediction',...
    ' '},...
    'Model Selection Help');



function smoothmeanrate_Callback(hObject, eventdata, handles)

newStrValue = get(hObject, 'String');
NewVal = str2double(newStrValue);

% Check if the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    errordlg('Please enter valid smoothing window.(positive number only)',...
        'Global Variable Error', 'modal');
    return;
end

% Assign the global variable TimeLag
global smooth_rt
smooth_rt = NewVal;

% --- Executes during object creation, after setting all properties.
function smoothmeanrate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothmeanrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calcalghelp.
function calcalghelp_Callback(hObject, eventdata, handles)
helpdlg({'Calculation algorithms: ',...
    ' ',...
    'When estimating STRFs of sensory neurons, the higher order statistics of ',...
    'stimulus is needed to be calculated, eg. auto correlation matrix, etc. ',...
    'STRFpak provides two methods: ',...
    '  1. space-time separable algorithm (fast)',...
    '  2. space-time nonseparable algorithm (more general).     ',...
    '   ',...
    '                                                          '},...
    'Calculation Algorithm HELP');

% --- Executes on button press in separable.
function separable_Callback(hObject, eventdata, handles)
global setSep
newVal = get(hObject, 'Value');
if newVal == 1
    setSep = 1;
    set(handles.nonseparable, 'Value', 0);
end

% --- Executes on button press in nonseparable.
function nonseparable_Callback(hObject, eventdata, handles)
global setSep
newVal = get(hObject, 'Value');
if newVal == 1
    setSep = 0;
    set(handles.separable, 'Value', 0);
end


% --- Executes on button press in outdirBrowser.
function outdirBrowser_Callback(hObject, eventdata, handles)
helpdlg({'Output Directory: ',...
    ' ',...
    ' The directory is used for saving the results and storing the ',...
    ' tempory files. It is provided from preprocessing stage.',...
    ' '},...
    'Output Directory Help');

function outdirshow_Callback(hObject, eventdata, handles)
global outputPath
set(handles.outdirshow, 'String', outputPath);

% --- Executes during object creation, after setting all properties.
function outdirshow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outdirshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in timelagUnit.
function timelagUnit_Callback(hObject, eventdata, handles)
v = get(handles.timelagUnit, 'value');
global TimeLagUnit;
if v == 1
    errordlg('Please choose right unit', 'TimeLag Unit error', 'modal');
    return;
elseif v == 2
    TimeLagUnit = 'msec';
elseif v == 3
    TimeLagUnit = 'frame';
end


% --- Executes during object creation, after setting all properties.
function timelagUnit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timelagUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function tolvalSlider_Callback(hObject, eventdata, handles)
global Tol_val
newSlider = round(get(hObject, 'Value'));
if newSlider <= 0
    set(handles.tolvalshow, 'String', Tol_val(handles.tolindex));
else
    newSlider = round(newSlider);
    set(handles.tolvalshow, 'String', Tol_val(newSlider));
    handles.tolindex = newSlider;
end
guidata(handles.figure1, handles);
plot_strfonly(handles)

% --- Executes during object creation, after setting all properties.
function tolvalSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolvalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function tolvalshow_Callback(hObject, eventdata, handles)
global Tol_val
maxlength = length(Tol_val);
newIndex = str2double(get(hObject,'String'));
if newIndex <=0 || newIndex > maxlength
    errordlg(['Tol_val_index is out of range. It need to be between 1 and ', num2str(maxlength)]);
    return;
else
    newIndex = round(newIndex);
    handles.tolindex = newIndex;
    %set(hObject, 'String', num2str(handles.Index));
    % Update slider value based on dataSet value
    set(handles.tolvalSlider, 'String', newIndex);
end
guidata(handles.figure1, handles);
plot_strfonly(handles);

% --- Executes during object creation, after setting all properties.
function tolvalshow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolvalshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function stdvalSlider_Callback(hObject, eventdata, handles)
global Std_val
newSlider = round(get(hObject, 'Value'));
if newSlider <= 0
    set(handles.stdvalshow, 'String', Std_val(handles.stdindex));
else
    newSlider = round(newSlider);
    set(handles.stdvalshow, 'String', Std_val(newSlider));
    handles.stdindex = newSlider;
end
guidata(handles.figure1, handles);
plot_strfonly(handles)


% --- Executes during object creation, after setting all properties.
function stdvalSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stdvalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function stdvalshow_Callback(hObject, eventdata, handles)
% hObject    handle to stdvalshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stdvalshow as text
%        str2double(get(hObject,'String')) returns contents of stdvalshow as a double
global Std_val
maxlength = length(Std_val);
newIndex = str2double(get(hObject,'String'));
if newIndex <=0 || newIndex > maxlength
    errordlg(['Std_val_index is out of range. It need to be between 1 and ', num2str(maxlength)]);
    return;
else
    newIndex = round(newIndex);
    handles.stdindex = newIndex;
    %set(hObject, 'String', num2str(handles.Index));
    % Update slider value based on dataSet value
    set(handles.stdvalSlider, 'String', newIndex);
end
guidata(handles.figure1, handles);
plot_strfonly(handles);

% --- Executes during object creation, after setting all properties.
function stdvalshow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stdvalshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upperfeq_Callback(hObject, eventdata, handles)
% hObject    handle to upperfeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperfeq as text
%        str2double(get(hObject,'String')) returns contents of upperfeq as a double

newStr_respsamprate = get(hObject, 'String');
NewVal = str2double(newStr_respsamprate);

% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid high freq (positive number only).',...
        'Variable Error', 'modal')
    return;
end

% Assign global variable to ampsamprate
global endFreq initialFreq

endFreq = NewVal;
if ~isempty(initialFreq)
    if initialFreq < endFreq
        plot_strfonly(handles);
        %     else
        %         errordlg('The lower freq should be smaller than upper freq.', 'Freq error','modal');
        %         return;
    end
end
% --- Executes during object creation, after setting all properties.
function upperfeq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lowerfeq_Callback(hObject, eventdata, handles)
newStr = get(hObject, 'String');
NewVal = str2double(newStr);

% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid low frequency (positive number only).',...
        'Variable Error', 'modal')
    return;
end

% Assign global variable to ampsamprate
global initialFreq endFreq
initialFreq = NewVal;
if ~isempty(endFreq)
    if initialFreq < endFreq
        plot_strfonly(handles);
        %     else
        %         errordlg('The lower freq should be smaller than upper freq.', 'Freq error','modal');
        %         return;
    end
end


% --- Executes during object creation, after setting all properties.
function lowerfeq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowerfeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
