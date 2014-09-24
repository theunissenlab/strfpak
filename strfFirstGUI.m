function varargout = strfFirstGUI(varargin)
% STRFFIRSTGUI Application M-file for strfFirstGUI.fig
%    FIG = STRFFIRSTGUI launch strfFirstGUI GUI.
%    STRFFIRSTGUI('callback_name', ...) invoke the named callback.
%
% Last Modified by GUIDE v2.5 19-Jan-2007 14:42:13
% STRFFIRSTGUI - the first window when the user run strfpak software.
%      FOUR stages:
%          INPUT
%              -- getfiles window
%              -- set calc. parameter window
%              -- display input window
%          CALCULATION
%              -- calculate
%              -- display stimulus statitics
%              -- display strf option
%          PREDICTION
%              -- get predicted data file
%              -- predict
%              -- display predicted result
%          VALIDATE
%              -- validate
%              -- display info value and cc
%              -- display bestfit among different tolvalue
%             STRFPAK: STRF Estimation Software
% Copyright ÃÂÃÂ¯ÃÂÃÂ¿ÃÂÃÂ½?003. The Regents of the University of California (Regents).
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

% Created by JXZ, Sept, 2002.
% Updated by JXZ, Sept, 2005.

%  Loads a preferences file, if available.
global loaded_preferences
pref_dir = find_prefs_dir;
filename = fullfile(pref_dir,'STRFPAK_script_parameters.mat');
if exist(filename,'file') & isempty(loaded_preferences)
    load(filename);
    disp(['loading preferences file ' filename '...']);
    loaded_preferences = 'Done';
end



fig = openfig(mfilename,'reuse');
% for resize property
hAxes = findall(fig,'type','axes');
hText  = findall(hAxes,'type','text');
hUIControls = findall(fig,'type','uicontrol');
set([hAxes; hText;...
    hUIControls],'fontname', 'times','units','normalized','fontunits','normalized');

% Generate a structure of handles to pass to callbacks, and store it.
handles = guihandles(fig);
guidata(fig, handles);

if nargout > 0
    varargout{1} = fig;
end

guidata(fig, handles);
if nargin == 0  % LAUNCH GUI



elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    % try
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    %     catch
    %         disp(lasterr);
    %     end

end
global ending
if strcmp(ending,'true')
    clear
    clear global
else
    update_guibuttons;
end


% --------------------------------------------------------------------
function varargout = preprocess_new_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Update mini command window's information
global preprocessOption
global rawDS DS predDS finalDS outputPath originalDS
DS = {};
predDS = {};
finalDS = {};
if exist(fullfile(outputPath,'finalResult_avgSpike1.mat'),'file') | exist(fullfile(outputPath,'finalResult_avgSpike2.mat'),'file')
    delete(fullfile(outputPath,'finalResult_avgSpike1.mat'));  % puts avgspike1 into workspace
    delete(fullfile(outputPath,'finalResult_avgSpike2.mat'));  % puts avgspike2 into workspace
end
if exist(fullfile(outputPath,'predResult_avgSpike1.mat'),'file') | exist(fullfile(outputPath,'predResult_avgSpike2.mat'),'file')
    delete(fullfile(outputPath,'predResult_avgSpike1.mat'));  % puts avgspike1 into workspace
    delete(fullfile(outputPath,'predResult_avgSpike2.mat'));  % puts avgspike2 into workspace
end
if exist(fullfile(outputPath,'best_strf.mat'),'file')
    delete(fullfile(outputPath,'best_strf.mat'));
    disp('Deleting best_strf.mat');
end
if exist(fullfile(outputPath,'spike_hashes'),'dir')
    rmdir(fullfile(outputPath,'spike_hashes'),'s');
end
if exist(fullfile(outputPath,'stim_hashes'),'dir')
    rmdir(fullfile(outputPath,'stim_hashes'),'s');
end


if ~isempty(rawDS)
    oldStrVal = get(handles.strf_MCW, 'String');
    newStrVal = {[' % PreprocessOption: ',preprocessOption],...
        [' % Now you are in PREPROCESSING stage.'],...
        [' % There are more options to do preprocessing. STRFPAK'],...
        [' % preprocesses your input data and save them to the files. '],...
        [' % After you have done that, you can choose some or all of inputs for estiamtion.']};
    handles.orgMCW = oldStrVal;
    guidata(h, handles);
    set(handles.strf_MCW, 'String', newStrVal);
    ttt = preprocessmenu;
    uiwait(ttt);
    originalDS = DS;
    delete(fullfile(outputPath,'*_mean_removed.mat'));

else
    errordlg('You need to select data files first.',...
        'Selection Error', 'modal')
    return;
end
handles.validationOK = 0;
handles.predictionOK = 0;
handles.getfileOK = 1;
guidata(h, handles);
strfFirstGUI;

% --------------------------------------------------------------------
function varargout = get_filesGUI_new_Callback(h, eventdata, handles, varargin)
% Update mini command window's information
global preprocessOption
if isempty(preprocessOption)
    preprocessOption = 'Already Preprocessed';
end
oldStrVal = get(handles.strf_MCW, 'String');
newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    [' % Now you are in GET_PROPRESSED_FILES stage.'],...
    [' % After it returns, next you need set parameters'],...
    [' % to calculate STRF. But now you can diplay input by'],...
    [' % clicking the Display Input button.']};
set(handles.strf_MCW, 'String', newStrVal);
handles.orgMCW = oldStrVal;
guidata(h, handles);

% Now call get_files GUI
%clear global DS predDS sDim
global outputPath
if length(outputPath) == 0
    get_filesGUI;
else
    answ2 = questdlg({['You have already loaded data.  Would you like' ...
        char(10)  'to unload current data and load new data?' ]},'Choose new data?','Load new data','Do not load anything','Load new data');
    switch answ2
        case 'Load new data'
            clear global outputPath
            clear global DS
            tt = get_filesGUI;
            uiwait(tt);
            handles.validationOK = 0;
            handles.predictionOK = 0;
            handles.getfileOK = 1;
            guidata(h, handles);
            strfFirstGUI;
    end
end

% --------------------------------------------------------------------
function varargout = strf_parameters_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Update mini command window's information
global preprocessOption
newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    ['% Now you are in SET PARAMETER stage.'],...
    ['% After you have set all the calc. parameters,'],...
    ['% you can go to Calcuation Stage. But now you'],...
    ['% can display input if you already get files.'],...
    ['% Otherwise, click GET_PREPROCESSED_Files button.']};
set(handles.strf_MCW, 'String', newStrVal);

% Set default calculation parameters
global Tol_val TimeLag smooth_rt Std_val timevary_PSTH
if isempty(Tol_val)
    Tol_val = [0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005];
end
if isempty(smooth_rt)
    smooth_rt = 41;
end

global sDim DS
if isempty(Std_val)
    if ~isempty(DS) & length(DS) ==1
        Std_val = [0];
    else
        Std_val = [0 0.25 0.5 0.75 1 1.5 2 4 6];
    end
end
if isempty(TimeLag)
    if ~isempty(sDim)
        if strcmp(sDim, '1-D')
            TimeLag = 300;
        else
            TimeLag = 200;
        end
    end
end

% Call calc_parameters to update calculation parameters
calc_parameter2GUI;

handles.setparameterOK = 1;
guidata(h,handles);

% --------------------------------------------------------------------
function varargout = displayrawstim_new_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------

global rawDS rawData
global preprocessOption

if rawData ==0
    preprocessOption = 'Already Preprocessed Data';
    errordlg('This option is for raw stimulus only.',...
        'Selection Error', 'modal')
    return;
elseif rawData == 1
    preprocessOption = 'Raw Data';
end

if ~isempty(rawDS)
    % Update mini command window's information
    newStrVal = {[' % PreprocessOption: ',preprocessOption],...
        ['% Now you are in DISPLAYINPUT stage.'],...
        ['% You can set the calc. parameters and'],...
        ['% reselect data files or you can go to '],...
        ['% Calcuation Stage.']};
    set(handles.strf_MCW, 'String', newStrVal);
    %displayinput_GUI;
    displayrawstim_GUI;
else
    errordlg('You need to select data files first.',...
        'Selection Error', 'modal')
    return;
end



% --------------------------------------------------------------------
function varargout = strf_calculate_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global outputPath DS originalDS predDS finalDS
if ~exist(fullfile(outputPath,'Global_Variables.mat'))
    save(fullfile(outputPath,'Global_Variables.mat'),'originalDS');
else
    save(fullfile(outputPath,'Global_Variables.mat'),'originalDS','-APPEND');
end
save(fullfile(outputPath,'Global_Variables.mat'),'DS','predDS','finalDS','-APPEND');

if (length(DS) + length(predDS) + length(finalDS)) ~= length(originalDS)
    if length(originalDS) > 0
        DS = originalDS;
    else
        if ~isempty(predDS)
            DS = {DS{:} predDS{:}};
        end
        if ~isempty(finalDS)
            DS = {DS{:} finalDS{:}};
        end
    end
    predDS = {};
    finalDS = {};
end

if exist(fullfile(outputPath,'finalResult_avgSpike1.mat'),'file') | exist(fullfile(outputPath,'finalResult_avgSpike2.mat'),'file')
    delete(fullfile(outputPath,'finalResult_avgSpike1.mat'));  % puts avgspike1 into workspace
    delete(fullfile(outputPath,'finalResult_avgSpike2.mat'));  % puts avgspike2 into workspace
    disp('Removing "finalResult_avgSpike" files, which get in the way.');
end
if exist(fullfile(outputPath,'predResult_avgSpike1.mat'),'file') | exist(fullfile(outputPath,'predResult_avgSpike2.mat'),'file')
    delete(fullfile(outputPath,'predResult_avgSpike1.mat'));  % puts avgspike1 into workspace
    delete(fullfile(outputPath,'predResult_avgSpike2.mat'));  % puts avgspike2 into workspace
    disp('Removing "predResult_avgSpike" files, which get in the way.');

end
if exist(fullfile(outputPath,'best_strf.mat'),'file')
    delete(fullfile(outputPath,'best_strf.mat'));
    disp('Deleting best_strf.mat');
end
%  First, let's figure out what kind of cross-validation to do.
message = ['Before calculation begins, please indicate' char(10) ...
    'if you would like to calculate STRFs using' char(10) ...
    'a double-jackknife (the default validation)' char(10) ...
    'or by manually selecting data to be left out' char(10) ...
    'for validation (requires more work; probably' char(10) ...
    'not a great idea unless there''s a specific' char(10) ...
    'reason for it).'];
answer3 = questdlg(message,'Decide how to validate','Jackknife','Manual','Jackknife');


if strcmp(answer3,'Jackknife')
    %%%$$$ begin select_validation
    global matchflg;

    global predDS predinitialFreq predendFreq predampsamprate sDim finalDS;
    global DS initialFreq endFreq ampsamprate;
    if (~isempty(predDS) & matchflg == 0)
        DS = {DS{:} predDS{:}};
        predDS = {};
    end
    if (~isempty(finalDS) & matchflg == 0)
        DS = {DS{:} finalDS{:}};
        finalDS = {};
    end
    matchflg = 1;

    if isempty(DS)
        errordlg('Please assign global variable DS first.',...
            'Global Variable Error', 'modal')
        return;
    end
    %predDS = DS;
    predinitialFreq = initialFreq;
    predendFreq = endFreq;
    predampsamprate = ampsamprate;
    %%%$$$ end select_validation
    global outputPath
    save(fullfile(outputPath,'STRFPAK_script_parameters.mat'), 'predinitialFreq' , 'predendFreq' , 'matchflg' , 'predampsamprate' , 'sDim'  , '-APPEND');
    add_to_STRFPAK_script('strfFirstGUI.m','select_validation');
    %  Note that the STRFPAK script here will generate its own predDS

elseif strcmp(answer3,'Manual')
    global matchflg outputPath predinitialFreq predendFreq predampsamprate sDim predDS finalDS
    if matchflg == 1
        predDS = {};
        finalDS = {};
    end
    matchflg = 0;

    ttt=select_validation_files;  %%% We still need to add STRFPAK_script support here...
    uiwait(ttt);
    save(fullfile(outputPath,'STRFPAK_script_parameters.mat'), 'predinitialFreq' , 'predendFreq' , 'matchflg' , 'predampsamprate' , 'sDim'  , '-APPEND');

end

% Check if preprocess Input or get_files have been chosen
% if not, can not go on to this step

% Set default calculation parameters
global Tol_val TimeLag smooth_rt Std_val timevary_PSTH
%disp(['Before, tol is ' num2str(Tol_val)]);
if isempty(Tol_val)
    Tol_val = [0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005];
end
if isempty(smooth_rt)
    smooth_rt = 41;
end

global sDim DS
if isempty(Std_val)
    if ~isempty(DS) & length(DS) ==1
        Std_val = [0];
    else
        Std_val = [0 1 2 6];
    end
end
if isempty(TimeLag)
    if ~isempty(sDim)
        if strcmp(sDim, '1-D')
            TimeLag = 300;
        else
            TimeLag = 200;
        end
    end
end
%calculation_GUI;
tt = calculation_GUI;
uiwait(tt);
strf_predict_Callback(h, eventdata, handles, varargin);

% disp(['After, tol is ' num2str(Tol_val)]);
handles.validationOK = 0;
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = strf_displaystimstat_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------

% Update mini command window's information
global preprocessOption
newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    ['% Now you are in DISPLAYSTIMSTAT stage.'],...
    ['% You can go to DISPLAYSTRF stage or go to PREDICT stage.']};
set(handles.strf_MCW, 'String', newStrVal);

% Check sptio-domain to display appropriately
displaystimstat2d_GUI;

% --------------------------------------------------------------------
function varargout = strf_displaystrf_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Update mini command window's information
global preprocessOption
newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    ['% Now you are in DISPLAYSTRF stage.'],...
    ['% You can go to DISPLAYSTIMSTAT stage or go to PREDICT stage.']};
set(handles.strf_MCW, 'String', newStrVal);

displaystrf_GUI;




% --------------------------------------------------------------------
function varargout = strf_predict_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Update mini command window's information
global preprocessOption
newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    ['% Now you are in Predict stage.'],...
    ['% You can go to displaypredstrf stage or go to Validate stage.']};
set(handles.strf_MCW, 'String', newStrVal);

guidata(h, handles);

%%%$$$ begin validation

% Get calc parameters
global TimeLag Tol_val sDim ampsamprate
global NBAND
global outputPath matchflg DS predDS

if matchflg
    predDS = DS;
end
% Could be new set of pred data
global predinitialFreq predendFreq

if isempty(predDS)
    errordlg('Please assign prediction data sets first.',...
        'Global Variable Error', 'modal');
    return
end

if isempty(Tol_val) | isempty(TimeLag)
    errordlg('Please assign Tol_val and TimeLag  first.',...
        'Global Variable Error', 'modal');
    return
end

global TimeLagUnit
if strcmp(TimeLagUnit, 'msec')
    global ampsamprate
    if isempty(ampsamprate)
        ampsamprate = 1000;
    end
    twindow = [-round(TimeLag*ampsamprate/1000) round(TimeLag*ampsamprate/1000)];
else
    twindow = [-TimeLag TimeLag];
end

% Calculate stim_avg if not calculated before
avg_done_check = fullfile(outputPath,'stim_avg.mat');
if ~exist(avg_done_check)
    stim_avg = cal_AVG(predDS, NBAND,1);    % this seems innefficient - how about having this somewhere
else
    loadpsth = load(fullfile(outputPath, 'stim_avg.mat'), 'stim_avg');
    stim_avg = loadpsth.stim_avg;

end

disp('Now doing prediction.');
if ~exist('running_in_script_mode','var')

    % show the pointer to see how prediction process goes.
    set(handles.figure1,'Pointer', 'watch');
end
running_flag = 1;
if running_flag == 1
    tempWait = waitbar(0,...
        sprintf('Validating neuron response for %d Tol values', length(Tol_val)));
end

fname = fullfile(outputPath,'predResult_avgSpike1.mat');
if exist(fname,'file')
    delete(fname);
end

fname = fullfile(outputPath,'predResult_avgSpike2.mat');
if exist(fname,'file')
    delete(fname);
end

global now_do_untouched
now_do_untouched = 'no';

for ntol=1:length(Tol_val)
    if running_flag == 1
        waitbar(ntol/length(Tol_val), tempWait);
    end
    if ~exist('running_in_script_mode','var')
        % Call cal_PredStrf to do prediction
        errFlg = cal_PredStrf(1,predDS, stim_avg,...
            fullfile(outputPath,['strfResult_Tol',num2str(ntol),'.mat']),...
            ntol, twindow, NBAND);
    else
        errFlg = cal_PredStrf(1,predDS, stim_avg,...
            fullfile(outputPath,['strfResult_Tol',num2str(ntol),'.mat']),...
            ntol, twindow, NBAND,running_in_script_mode);
    end



    % Check if cal_PredStrf normally ends
    if errFlg == 1
        set(handles.figure1, 'Pointer', 'Arrow');
        return
    end
end
if running_flag == 1
    close(tempWait)
end
%%%$$$ end validation
global outputPath
add_to_STRFPAK_script('strfFirstGUI','validation');

% Could be new set of pred data
global predDS predinitialFreq predendFreq


% Set pointer to Arrow to show it ends normally
set(handles.figure1, 'Pointer', 'Arrow');

save(fullfile(outputPath, 'pred_Variables.mat'), 'predDS',...
    'predinitialFreq', 'predendFreq');

disp('Done Validation.')
handles.predictionOK = 1;
tt = msgbox('Done Validation Stage.', 'modal');
uiwait(tt);
strfFirstGUI;
if nargout > 0
    varargout{1} = 'Done prediction.';
end

% --------------------------------------------------------------------
function varargout = strf_displaypredstrf_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Update mini command window's information
% if ~isfield(handles, 'predictionOK')
%     errordlg(['Please do validation/cross validation first. Thank you!'],...
%         '!!DO Validation First!!','modal');
%     return;
% end
% Update mini command window's information
global preprocessOption

newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    ['% Now you are in displayPredStrf stage.'],...
    ['% You can go to the Validate stage.']};
set(handles.strf_MCW, 'String', newStrVal);

displaypredstrf_GUI;

% --------------------------------------------------------------------
function varargout = strf_validate_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Update mini command window's information
global preprocessOption
newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    ['% Now you are in VALIDATE stage.'],...
    ['% After you have done, you can go to the DISPLAY CC/INFO stage'],...
    ['% or the DISPLAY BESTFIT stages.']};
set(handles.strf_MCW, 'String', newStrVal);

global ampsamprate sDim outputPath
global Tol_val Std_val
if isempty(Tol_val) | isempty(ampsamprate) | isempty(sDim) | isempty(Std_val)
    errordlg('There is no value for global Tol_val and ampsamprate.',...
        'Global Variable Error', 'modal');
    return
end
set(handles.figure1,'Pointer', 'watch');

% Jxz, 9/15/2005
% First check whether you have proper trial numbers (must > 1) for validation.
spike_psth = Check_And_Load(fullfile(outputPath,'spike_psth_count.mat'));
ntrials_proper=spike_psth(end);


% If trial number is less than one, we can not go ahead for validation.
if ntrials_proper <= 1
    errordlg(['Since the most common trial number of prediction data is only one,',...
        'the current validation methods are not allowed. Sorry!']);
    handles.validationOK = 0;
    guidata(h, handles);
    return;
end

spike_psth=spike_psth(1:end-1);
ntrials_index=find(spike_psth>=ntrials_proper);
ntrials = spike_psth;

% Variable binWindow used for validation
% The binWindow really depends on  width of autocorrelation of psth
if ~isempty(sDim)
    if strcmp(sDim, '1-D')
        binWindow = 128;
    else
        binWindow = 20;
    end
end


disp('Now doing validation.');


% JXZ, 8/23/2005
%  Ask the user to input smoothing range for cc
prompt={'Smallest smoothing window width (in ms):', 'Smoothing window width step (in ms):',...
    'Largest smoothing window width (in ms):','Fixed smoothing window width (in ms):'};
def = {'5','8','29','21'};
dlgTitle='Hanning window range for correlation coefficient';
lineNo=1;

% picture feature
AddOpts.Resize='on';
AddOpts.WindowStyle='normal';
AddOpts.Interpreter='tex';
smooth_V =inputdlg(prompt,dlgTitle,lineNo,def,AddOpts);
% show the waitbar to see how calculate process goes.
pause(.3);

global smoothVect psth_smoothconst
smoothVect = [str2double(smooth_V{1}) str2double(smooth_V{2}) str2double(smooth_V{3})];
psth_smoothconst = str2double(smooth_V{4});


% Loop through all the results for all the tol values
% JXZ 7/13/2005
%   Add extra loop for Std_val list

%%%$$$ begin goodness_of_fit
global now_do_untouched
now_do_untouched = 'no';
% Load part I of PSTH of predict response files
spredresult = fullfile(outputPath, 'predResult_avgSpike1.mat');
avgSpike1 = Check_And_Load(spredresult);

% Load part II of PSTH of predict response files
spredresult = fullfile(outputPath, 'predResult_avgSpike2.mat');
avgSpike2 = Check_And_Load(spredresult);

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
            = do_locally_cached_calc(get_local_cache_dir,'cal_Validate',spike_est1, spike_est2,ntol, istd, ampsamprate, binWindow); %ampsamprate);

    end % END of ntol
end  % END of istd


% save to the file for each data pair
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

%%%$$$ end goodness_of_fit
add_to_STRFPAK_script('strfFirstGUI.m','goodness_of_fit');

set(handles.figure1, 'Pointer', 'Arrow');
disp('Done validation.')
msgbox('Done Validation Stage.', 'modal')

% --------------------------------------------------------------------
function varargout = goodnessfitting_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Update mini command window's information
global preprocessOption
newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    ['% Now you are in DISPLAYINFO/CC stage.'],...
    ['% After you have done, you can go anywhere from here.']};
set(handles.strf_MCW, 'String', newStrVal);

global outputPath;
% First check whether you have proper trial numbers (must > 1) for validation.
spikefile = fullfile(outputPath,'spike_psth_count.mat');
if not(exist(spikefile))
    errordlg('No validation results yet, please do validation first', 'Display InfoCC error','modal');
    return;
end
spike_psth = Check_And_Load(spikefile);
ntrials_proper=spike_psth(end);

if ntrials_proper <= 1
    global novalidation
    novalidation = 1;
    errordlg(['Since the most common trial number of prediction data is only one,',...
        'the current validation methods are not allowed. But you can go ahead to do prediction!']);
    return;
end
% Release Version 1.1
%displayinfocc_GUI_before;
tt = displayinfocc_GUI;
uiwait(tt);
strfFirstGUI;
handles.validationOK = 1;
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = strf_displaybestfit_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Update mini command window's information
% if ~isfield(handles, 'validationOK')
%     errordlg(['Please do validation first!']);
%     return;
% end

global preprocessOption
newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    ['% Now you are in DISPLAYBESTFIT stage.'],...
    ['% You can go anywhere from here.']};
set(handles.strf_MCW, 'String', newStrVal);
global outputPath;
% First check whether you have proper trial numbers (must > 1) for validation.
spike_psth = Check_And_Load(fullfile(outputPath,'spike_psth_count.mat'));
ntrials_proper=spike_psth(end);

if ntrials_proper <= 1
    global novalidation
    novalidation = 1;
    errordlg(['Since the most common trial number of prediction data is only one,',...
        'the current validation methods are not allowed. But you can go ahead to do prediction!']);
    return;
end
displaybestfit_GUI;
disp('Done displaying best STRF.');
%msgbox('Coming soon...', 'modal');


% --------------------------------------------------------------------
function varargout = loadprevresult_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global outputPath
if ~isempty(outputPath)
    defaultdir = outputPath;
else
    defaultdir = pwd;
end

outputdir = uigetdir(defaultdir, 'Pick the output directory');

if not(exist(outputdir,'dir'))
    text = ['Directory not found... exiting.'  char(10) ...
        'On some Linux systems the function "uigetdir" does not work.'];
    errordlg(text);
    return
end
outputPath = outputdir;
if exist(fullfile(outputPath,'finalResult_avgSpike1.mat'),'file') | exist(fullfile(outputPath,'finalResult_avgSpike2.mat'),'file')
    handles.predictionOK = 1;
end
if exist(fullfile(outputPath,'best_strf.mat'),'file')
    handles.validationOK = 1;
end
%guidata(fig, handles);

load(fullfile(outputPath,'STRFPAK_script_dataset.mat'));
load(fullfile(outputPath,'STRFPAK_script_parameters.mat'));

% Load old copy of global variables
tempGV = fullfile(outputPath, 'Global_Variables.mat');
if not(exist(tempGV, 'file'))
    errordlg(['No input files yet. Are you sure you aready',...
        ' ran STRFPAK?.'])
    return
end
load(tempGV)
% Load old copy of global variables
tempGV = fullfile(outputPath, 'pred_Variables.mat');
if not(exist(tempGV, 'file'))
    errordlg('No prediction input yet.')
    return
end
load(tempGV)
tt = msgbox('Directory has been verified. Please go ahead');
uiwait(tt);
strfFirstGUI;
% --------------------------------------------------------------------
function varargout = strf_help_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
ttlStr = get(handles.figure1, 'Name');
hlpStr = [...
    'The STRFPAK software is a matlab toolbox that can be used to estimate the   '
    'stimulus-response transfer function of sensory neurons. This version of     '
    'STRFPAK is designed as a framework in the study of STRF of sensory neurons. '
    'It currently implements modified reverse correlation technique to estimate  '
    'the STRFs and implements Short-time Fourier transform and Wavelet transform '
    'algorithms to preprocess 1-D stimuli. It also implements two validation     '
    'methods: correlation coefficients and information values. It includes demo  '
    'data under the source codes directory and also provides limited help with   '
    'references on each GUI window. Right now it includes four stages(modules):  '
    '                                                                            '
    '        LOAD INPUT                                                          '
    '            - load raw stimulus and spike train or preprocessed data sets   '
    '            1. Load Data: load stimulus and response window.                '
    '            2. Display Stim Only: This window is only for raw stimulus      '
    '               displaying. It is useful when you want to display what new   '
    '               stimulus looks like.                                         '
    '            3. Display Raw Data: It displays raw stimuli and its correspond '
    '               -ed spike train. Here raw stimuli mean songwave files for    '
    '               auditory data and movie frames for vision data.              '
    '        Preprocess                                                          '
    '            - after you load the new data set, you can preprocess your raw  '
    '              data set if the loaded data is raw data.                      '
    '            1. Preprocess Data: preprocessing raw data and converting to    '
    '                proper data format for STRFPAK.                             '
    '            2. DISPLAY Preprocessed_Data: displaying loaded preprocessed    '
    '               data visually.                                               '
    '                                                                            '
    '        ESTIMATE                                                            '
    '            - Estimate the STRFs and display the intermediate results.      '
    '            1. CALCULATE: calculating the STRF.                             '
    '            2. DISPLAY STIMSTAT: displaying stimulus statistics and STA.    '
    '            3. DISPLAY STRF: displaying the STRFs and regularizated STRFS.  '
    '                                                                            '
    '         Validate                                                           '
    '            - Get Valfiles: When this button is clicked, you will be asked  '
    '              the validation options: validation on new data set or do cross'
    '              validation.                                                   '
    '            - Validation/Cross Validation.                                  '
    '            - Display prePSTH: display the validation results.              '
    '            - Goodness of fitting: Display how good it fits quantatively.   '
    '                                                                            '
    '        Predict                                                             '
    '            - Display the Best estimated STRF. The best STRF is chosen based'
    '              on the predicted information values.                          '
    '            - Get New Stim Only: To load totally new stimulus. You can      '
    '              predict the neural response on the new stimulus using the     '
    '              estimated STRFs.                                              '
    '            - Predict and Display Results: Show the STRF, the new Stimulus  '
    '              and predicted neural response in one window. The results are  '
    '              saved under output direcotry. They can be used for the neural '
    '              prediction challenge project.                                 '
    '                                                                            '
    '                                                                            '
    ' MINI COMMAND WINDOW: The left panel of the main window is the mini         '
    ' command/help box. After you click any button from above, it shows what     '
    ' process you just click and where you can go from there.                    '
    '                                                                            '
    ' LOAD PREV RESULT: This option is very useful once you want to display      '
    ' previous interactive results or batch results. After you click it, it  will'
    ' ask you the output path.                                                   '
    '                                                                            '
    ' HELP: This gives you a brief explaination for this main window.  For more  '
    ' detailed manual or documents, please go to http://strfpak.berkeley.edu.    '
    '                                                                            '
    ' CLOSE: This button close all the figures and windows.                      '
    '                                                                            '
    ' Updated by JXZ, May, 2006.                                                 '
    '                                                                            '
    '                                                                            '];
myFig = handles.figure1;
helpwin(hlpStr, ttlStr);
%


% --------------------------------------------------------------------
function varargout = strf_close_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% close this main window
pref_dir = find_prefs_dir;
filename = fullfile(pref_dir,'STRFPAK_script_parameters.mat');
global outputPath
sourcename = fullfile(outputPath,'STRFPAK_script_parameters.mat');

answ2 = questdlg({['STRFPAK can save your parameters in the following directory:' ...
    char(10) pref_dir char(10) ...
    'Would you like to save your parameters as defaults?' ]},'Save Parameters for Later?','Save Parameters','Do Not Save','Save Parameters');
switch answ2
    case 'Save Parameters'
        copyfile(sourcename,filename,'f');
    otherwise
        disp('Preferences not saved.');
end

delete(handles.figure1);


global ending
ending = 'true';


% --- Executes on button press in ppdata_display.
function ppdata_display_Callback(hObject, eventdata, handles)
global preprocessOption
global DS sDim
if ~isempty(DS) | ~isempty(sDim)
    % Update mini command window's information
    newStrVal = {[' % PreprocessOption: ',preprocessOption],...
        ['% Now you are in DISPLAYINPUT stage.'],...
        ['% You can set the calc. parameters and'],...
        ['% reselect data files or you can go to '],...
        ['% Calcuation Stage.']};
    set(handles.strf_MCW, 'String', newStrVal);
    displayppinput_GUI;
else
    errordlg('You need preprocess data first.',...
        'Selection Error', 'modal')
    return;
end


% --- Executes on button press in displayrawstimpsth_new.
function displayrawstimpsth_new_Callback(hObject, eventdata, handles)

global rawDS rawData
global preprocessOption

if rawData ==0
    preprocessOption = 'Already Preprocessed Data';
    errordlg('This option is for raw stimulus only.',...
        'Selection Error', 'modal')
    return;
elseif rawData == 1
    preprocessOption = 'Raw Data';
end

if ~isempty(rawDS)
    % Update mini command window's information
    newStrVal = {[' % PreprocessOption: ',preprocessOption],...
        ['% Now you are in DISPLAYINPUT stage.'],...
        ['% You can set the calc. parameters and'],...
        ['% reselect data files or you can go to '],...
        ['% Calcuation Stage.']};
    set(handles.strf_MCW, 'String', newStrVal);
    displayinput_GUI;
else
    errordlg('You need select data files first.',...
        'Selection Error', 'modal')
    return;
end

% --- Executes on button press in getpredstim.
function getpredstim_Callback(hObject, eventdata, handles)
% Update mini command window's information
global preprocessOption
newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    ['% Now you are in "Get New Stim" stage.'],...
    ['% This option is useful when you participate the NPC contest.'],...
    };
set(handles.strf_MCW, 'String', newStrVal);
ttt =loadstimonly_GUI;
uiwait(ttt);

% --- Executes on button press in predNdisplay.
function predNdisplay_Callback(hObject, eventdata, handles)
% Update mini command window's information
global preprocessOption
newStrVal = {[' % PreprocessOption: ',preprocessOption],...
    ['% Now you are in "Predict and Display Results" stage.'],...
    ['% This option is useful when you participate the NPC contest.'],...
    ['% You can take the result files to submit at http://neuralprediction.berkeley.edu.'],...
    };
set(handles.strf_MCW, 'String', newStrVal);
predndisplay_GUI;

% --- Executes on button press in Special_options.
function Special_options_Callback(hObject, eventdata, handles)
% hObject    handle to Special_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ttt = special_optionsGUI;
uiwait(ttt);


global allow_negative_rates use_alien_space alien_space_file outputPath
if allow_negative_rates | use_alien_space
    set(handles.Special_options_warning,'Visible','on');
else
    set(handles.Special_options_warning,'Visible','off');
end
save(fullfile(outputPath,'STRFPAK_script_parameters.mat'),'allow_negative_rates','use_alien_space','alien_space_file','-APPEND');
% --------------------------------------------------------------------
% END of strfFirstGUI
% --------------------------------------------------------------------


