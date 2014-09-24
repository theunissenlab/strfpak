function varargout = displayppinput_GUI(varargin)
% DISPLAYPPINPUT_GUI Application M-file for songwave_specgram.fig
%    FIG = DISPLAYPPINPUT_GUI launch songwave_specgram GUI.
%    DISPLAYPPINPUT_GUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 10-Apr-2006 15:49:36

if nargin == 0  % LAUNCH GUI

	%fig = openfig(mfilename,'reuse');
    fig = openfig(mfilename);

    % for resize property
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
           hUIControls],'units','normalized','fontunits','normalized');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    handles.index = 1;
	guidata(fig, handles);
  
    % Check whether the input is right for preprocessing
    global DS sDim
    if isempty(DS)
        tt=errordlg('You need preprocess your data first.',...
            'Display Type Error', 'modal');
        uiwait(tt);
        delete(handles.figure1);
        return;
    else
        [path, name, ext] = fileparts(DS{handles.index}.stimfiles); 
        set(handles.filename, 'String',[name ext]);
        [path, name, ext] = fileparts(DS{handles.index}.respfiles); 
        set(handles.rfilename, 'String',[name ext]);
        set(handles.inputdirshow, 'String', path);
    end
    
    if isempty(sDim)
        set(handles.spatialdomain, 'value', 1);
    elseif strcmp(sDim, '1-D')
        set(handles.spatialdomain, 'value', 2);
    elseif strcmp(sDim, '2-D')
        set(handles.spatialdomain, 'value', 3);
    end
    
    global preprocessOption
    if ~isempty(preprocessOption)
        set(handles.preprocessoption_show, 'String', preprocessOption);
    else
        preprocessOption = 'Already preprocessed data';
        set(handles.preprocessoption_show, 'String', preprocessOption);
    end
    
    % Set values for dataset
    if length(DS) ==1
        set(handles.datasetSlider, 'min', 0, 'max', length(DS), 'value',0,'sliderstep', [1 1]);
    elseif length(DS) > 1
        set(handles.datasetSlider, 'Min', 1, 'Max', length(DS),'value',1,...
            'sliderstep', [1/(length(DS)-1) 1/(length(DS)-1)]);
    end
    set(handles.DataSet, 'string', num2str(handles.index));
    
    global outputPath
    if ~isempty(outputPath)
        set(handles.outdirshow, 'visible', 'on', 'String', outputPath);
        set(handles.outdirBrowser, 'visible', 'on');
    end
    
    global stimsamprate respsamprate
    if ~isempty(stimsamprate)
        set(handles.stimsamprateD, 'string', stimsamprate);
    end
    if ~isempty(respsamprate)
        set(handles.respsamprateD, 'string', respsamprate);
    end
    
    displayspecgram(handles);
    
	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%
% --------------------------------------------------------------------
function varargout = startidx_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    global startidx_v stopidx_v
    newStr = get(h, 'String');
    NewVal = str2double(newStr);

    % Check that the entered value falls within the allowable range
    
    if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal) | startidx_v > stopidx_v
        % Set default value for resp sample rate
        errordlg('Please enter valid starting frame (positive number only).',...
            'Variable Error', 'modal')
        return; 
    end
    
    % Assign global variable to StimSampRate
    startidx_v = NewVal;
    
    

% --------------------------------------------------------------------
function varargout = stopidx_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
     global stopidx_v startidx_v
    newStr = get(h, 'String');
    NewVal = str2double(newStr);

    % Check that the entered value falls within the allowable range
    
    if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal) | NewVal < startidx_v
        % Set default value for resp sample rate
        errordlg('Please enter valid starting frame (positive number only).',...
            'Variable Error', 'modal')
        return; 
    end
    
    % Assign global variable to StimSampRate
    stopidx_v = NewVal;

% --------------------------------------------------------------------
function varargout = bwindow_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
     global bwindow_v
     
     newStr = get(h, 'String');
     NewVal = str2double(newStr);

     % Check that the entered value falls within the allowable range
     if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
         % Set default value for resp sample rate
         errordlg('Please enter valid frequency bandwith (positive number only).',...
                'Variable Error', 'modal')
         return; 
     end
     
    bwindow_v = NewVal;
     
% --------------------------------------------------------------------
function varargout = btakepower_Callback(h, eventdata, handles, varargin)

     global btakepower_v
     
     newStr = get(h, 'String');
     NewVal = str2double(newStr);

     % Check that the entered value falls within the allowable range
     if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
         % Set default value for resp sample rate
         errordlg('Please enter valid frequency bandwith (positive number only).',...
                'Variable Error', 'modal')
         return; 
     
     end
     
     btakepower_v = NewVal;
     
     
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
     displayspecgram(handles);
     
% --------------------------------------------------------------------
function varargout = fstep_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    global fstep;
    
    newStr = get(h, 'String');
    NewVal = str2double(newStr);

    % Check that the entered value falls within the allowable range
    
    if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
        % Set default value for resp sample rate
        errordlg('Please enter valid frequency bandwith (positive number only).',...
            'Variable Error', 'modal')
        return; 
    end
    
    % Assign global variable to StimSampRate
    fstep = NewVal;
    
% --------------------------------------------------------------------
function varargout = loadsongfile_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    % Load raw files by calling loadrawfile
    tt=loadrawfile;
    uiwait(tt);
    global DS
    if ~isempty(DS)
        [path, name, ext] = fileparts(DS{handles.index}.stimfiles); 
        set(handles.filename, 'String',[name ext]);
        [path, name, ext] = fileparts(DS{handles.index}.respfiles); 
        set(handles.rfilename, 'String',[name ext]);
    end
    guidata(h, handles);
% --------------------------------------------------------------------
function varargout = loadresp_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
   %msgbox('Coming soon...', 'modal');
   loadrawfile;

% --------------------------------------------------------------------
function varargout = display_Callback(h, eventdata, handles, varargin)
    % 1. Load display result
    handles.index = 1;
    guidata(h, handles);
    displayspecgram(handles);

% --------------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
global ampsamprate respsamprate
if isempty(respsamprate)
    warndlg('Response sampling rate is not set and will be set to its default value: 1000Hz',...
        'Resp Samprate Info','modal');
    respsamprate = 1000;
end
delete(handles.figure1);
    

% --------------------------------------------------------------------
function varargout = filteroption_Callback(h, eventdata, handles, varargin)
    v = get(handles.filteroption, 'value');
    option = get(handles.filteroption, 'String');
    option = deblank(option(v,:));
    
    global filteroption
    
    if strcmp(option, 'Linear')
        filteroption = 1;
    elseif strcmp(option, 'Logarithmic')
        filteroption = 2;
    end
    
% --------------------------------------------------------------------
function varargout = prev_Callback(h, eventdata, handles, varargin)

    global DS;
    index = handles.index;
    % update index
	i = index - 1;
	% If the index is less then one then set it equal to the index of the
	% last element in the Addresses array
	if i < 1
	    i = length(DS);
	end

	handles.index = i;
    guidata(h, handles);
    
    displayspecgram(handles)
    
% --------------------------------------------------------------------
function varargout = next_Callback(h, eventdata, handles, varargin)

    
    index = handles.index;
    % update index
	i = index + 1;
	% If the index is less then one then set it equal to the index of the
	% last element in the Addresses array
	if i > length(DS)
	    i = 1;
	end

	handles.index = i;
    guidata(h, handles);
    
    displayspecgram(handles)

% --------------------------------------------------------------------  
function displayspecgram(handles)
    
    %  Prepare for axes
    %
    global DS psth_smooth ampsamprate
    global outputPath sDim
    
    if isempty(psth_smooth)
        psth_smooth = 21;
    end
    if isempty(ampsamprate)
        ampsamprate = 1000;
    end
    
    forward = Check_And_Load(DS{handles.index}.stimfiles);
    global NBAND
    NBAND = size(forward,1);
    % Reupdate length of data pair
    DS{handles.index}.nlen = min(DS{handles.index}.nlen, size(forward,2));
    
    rawResp = Check_And_Load(DS{handles.index}.respfiles);
    % Validate the MAT-file

    if isempty(DS{handles.index}.ntrials)
        DS{handles.index}.ntrials = length(rawResp);
    end
    if iscell(rawResp)
        spiketrain = zeros(length(rawResp),DS{handles.index}.nlen);
        for trial_ind =1:DS{handles.index}.ntrials

            spiketrain(trial_ind, rawResp{trial_ind}) = ones(1,length(rawResp{trial_ind}));
        end

        newpsth = spiketrain;
%         newpsth = newpsth'; % make sure new response data is trials x T.
%         newpsth(find(newpsth < 0)) = 0;
    else
        newpsth = rawResp;
    end
    [xs, ys] = size(newpsth);

    if xs > ys   % make sure response data is trials X TimeFrame
        newpsth = newpsth';
        DS{handles.index}.nlen = xs;
    else
        DS{handles.index}.nlen = ys;
    end
    
    if xs > 1 & ys >1   % check whether multiple trials
        newpsth = mean(newpsth);
    end
    rgoodidx = find(~isnan(newpsth));
    newpsth = newpsth(rgoodidx);
    global  allow_negative_rates
    if isempty(allow_negative_rates)
        allow_negative_rates = 0;
    end
    if ~allow_negative_rates
        newpsth(find(newpsth < 0)) = 0;
    end

    % Check if raw data is loaded as preprocessed data
    if length(size(forward)) > 2
        errordlg('This window is for preprocessed data only. Your data is identified as raw data.',...
            'modal');
        return;
    end
    [xm nt] = size(forward);
    
    maxforward = max(max(forward));
    minforward = min(min(forward));
    absforward = max(abs(minforward),abs(maxforward));
    
    % Now display input stimulus file
    
    axes(handles.stim);
    global sDim
    if strcmp(sDim, '0-D')
        ttH = plot(forward);
    else
        ttH = imagesc(forward);
        %ttH = imagesc(forward,[-absforward absforward]);
    end
    axis tight;
    axis xy;
    title('Stimulus')
    
    set(handles.smoothconst, 'String', psth_smooth);
    set(handles.stim, 'visible', 'on');
    set(handles.smoothpsth, 'visible','on');
    set(handles.smoothconst, 'visible','on');
    set(handles.mstag, 'visible', 'on');
   
    spikeavg = newpsth;
    nt = min(nt, length(spikeavg));
    % Junli: 9/29/04
    % Check if psth_smoothwindow is even or odd and to set cutsize
    if psth_smooth == 0
        psth_smooth = 1;
    end

    if mod(psth_smooth, 2) == 0
        cutsize = 0;
    else
        cutsize = 1;
    end
    
    global outputPath
    [ss, Avg_psth, pp, constmeanrate] = cal_AVG(DS, NBAND,1);
    %loadpsth = load(fullfile(outputPath, 'stim_avg.mat'), 'Avg_psth', 'constmeanrate');
    vary_rate = Avg_psth;
    const_rate = constmeanrate;
    clear ss pp
    
    %Show smoothed version of psth at window = 15
    xlabelrange=1:nt;
    halfwinsize = floor(psth_smooth/2);
    wind1 = hanning(psth_smooth)/sum(hanning(psth_smooth));
    svagsm=conv(newpsth(1:nt),wind1);
    smoothpsth = svagsm(halfwinsize+cutsize:length(svagsm)-halfwinsize);
    axes(handles.psth);
    if strcmp(sDim, '1-D')
        plot(xlabelrange,smoothpsth(xlabelrange)*ampsamprate,'b'); hold on;
        plot(xlabelrange,(smoothpsth(xlabelrange)-vary_rate(handles.index,xlabelrange))*ampsamprate, 'r'); hold on;
        plot(xlabelrange, (smoothpsth-const_rate)*ampsamprate, 'g');
        axis tight;
    else
        plot(xlabelrange,smoothpsth(xlabelrange)*ampsamprate,'b'); hold on;
        plot(xlabelrange,(smoothpsth(xlabelrange)-vary_rate(handles.index,1:size(smoothpsth(xlabelrange),2)))*ampsamprate, 'r'); hold on;
        plot(xlabelrange,(smoothpsth(xlabelrange)-const_rate)*ampsamprate, 'g');
        
    end
    ylabel('Smoothed PSTH (spikes/s)')
    axis tight;
    %xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
    legend('PSTH', 'Time-varying removed', 'Mean removed');
    hold off;
    
    global stimsamprate sDim
    if ~isempty(stimsamprate)
        if sDim == '2-D'
            timebin = round(1000/stimsamprate);  % in ms
            xlabel(['Time (in ', num2str(timebin), ' ms)']);
            set(handles.stimsamprateD, 'String', stimsamprate);
        else

            xlabel('Time (in ms)');
        end
    else
        xlabel('Frames');
    end
   
    
    if ~isempty(DS)
        [path, name, ext] = fileparts(DS{handles.index}.stimfiles); 
        set(handles.filename, 'String',[name ext]);
        [path, name, ext] = fileparts(DS{handles.index}.respfiles); 
        set(handles.rfilename, 'String',[name ext]);
        set(handles.inputdirshow, 'String', path);
    end
    
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
            '   STRFPAK: Display Preprocessed Input Window                  '
            '                                                               '
            ' Only after the input preprocessed data has been loaded, you   '
            '  can display the stimulus and its response data graphically.  '
            ' On this window, there are two panels. The left panel shows    '
            ' the plots. The right panel shows all the related information  '
            ' about the plots.                                              '
            '                                                               '
            '     1. Plot of stimulus:                                      '
            '             X-axis is TIME (in milliseconds);                 '
            '             Y-axis is Spatial domain (will be frequency in Hz '
            '             if spatio-domain is 1-Dimension and wrapped space '
            '             axes if 2-Dimension).                             '
            '     2. Plot of psth (TIME)                                    '
            '        PSTH is computed as average across all the trials and  '
            '        is shown in blue line. The green line shows the psth   '
            '        removed constant mean rate. The red lines shows the    '
            '        psth removed time-varying mean rate.                   '
            '                                                               '
            '     Spatial Domain: shows the dimension size of stimulus.     '
            '     Stim files: shows the filename of the stimulus.           '
            '     Resp files: shows the filename of the response.           '
            '     Total_trials: shows total number of trials for response.  '
            '     Prev_10Trials:  if total_trials is more than 10 trials,   '
            '                   the left middle panel shows first 10 trials.'
            '     Next_10Trials: show next 10 or less trials.               '
            '     Smooth_psth: the parameter used for smoothing psth window,'
            '          you can type different values for different psth plot'
            '     Display whole psth: (for 2-D only) It plots the psth with '
            '     the whole time duration instead of 10 frames only.        '
            '                                                               '
            '                                                               '
            'You can also see which data set you are displaying from stim   '
            'file text field and response file text field.                  '
            '                                                               '
            'Zoom in/out properties:                                        '
            '   MATLAB has built-in function for zooming the figure. If you '
            '   want to zoom in/out the figure, you go to VIEW menu and     '
            '   then choose FIGURE TOOLBAR. Then you choose ZOOM IN/OUT     '
            '   button and select the part of the figure for zooming.       ' 
            '                                                               '
            ' Updated by Junli, May 2006.                                   '
            '                                                               '];
    myFig = handles.figure1;
    helpwin(hlpStr, ttlStr);


function btempnl_Callback(hObject, eventdata, handles)

newStr_respsamprate = get(hObject, 'String');
NewVal = str2double(newStr_respsamprate);

% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid low frequency (positive number only).',...
        'Variable Error', 'modal')
    return; 
end

% Assign global variable to bwindow
global initialFreq;
if isempty(initialFreq)
    initialFreq = NewVal;
else
    initialFreq = NewVal;
    compsave_Callback(hObject, eventdata, handles);
end

% --- Executes on slider movement.
function datasetSlider_Callback(hObject, eventdata, handles)
global DS
maxlength = length(DS);
newSlider = round(get(hObject, 'Value'));
if newSlider <= 0 | newSlider > maxlength
    set(handles.DataSet, 'String', handles.index);
else
    set(handles.DataSet, 'String', newSlider);
    handles.index = newSlider;
end
guidata(handles.figure1, handles);
displayspecgram(handles)

% --- Executes during object creation, after setting all properties.
function datasetSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datasetSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function DataSet_Callback(hObject, eventdata, handles)
global DS
maxlength = length(DS);
newIndex = str2double(get(hObject,'String'));
if newIndex <=0 || newIndex > maxlength
    warndlg(['DataSet is out of range. It need to be between 1 and ', num2str(maxlength)]);
    set(hObject, 'String', num2str(handles.index));
else
    handles.index = newIndex;
    %set(hObject, 'String', num2str(handles.Index));
    % Update slider value based on dataSet value
    set(handles.datasetSlider, 'String', newIndex);
end
guidata(handles.figure1, handles);
displayspecgram(handles)
     
% --- Executes during object creation, after setting all properties.
function DataSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in outdirBrowser.
function outdirBrowser_Callback(hObject, eventdata, handles)
 helpdlg({'Output Directory: ',...
      ' ',...
      ' The directory is used for saving the results and storing the ',...
      ' tempory files. It is provided from preprocessing stage.',...
      ' '},...
      'Output Directory Help');
% global outputPath
% outputdir = uigetdir(pwd, 'Pick the response data directory');
% if not(exist(outputdir,'dir'))
%      errordlg('There is no selected directory', 'Output Dir Error', 'modal');
%      return;
% end
% set(handles.outdirshow, 'String', outputdir);
% outputPath = outputdir;

function outdirshow_Callback(hObject, eventdata, handles)
global outputPath
set(handles.outdirshow, 'String', outputPath);
 

% --- Executes during object creation, after setting all properties.
function outdirshow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spatialdomain_Callback(hObject, eventdata, handles)
v = get(handles.spatialdomain, 'value');
global sDim;
if v == 1
    sDim = '0-D';
elseif v == 2
    sDim = '1-D';
elseif v == 3
    sDim = '2-D';
end

% --- Executes during object creation, after setting all properties.
function spatialdomain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spatialdomain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1);



function stimsamprateD_Callback(hObject, eventdata, handles)

newStrValue = get(hObject, 'String');
NewVal = str2double(newStrValue);

% Check if the entered value falls within the allowable range

if (NewVal< 0) | ~isnumeric(NewVal)
    % Set default value for step
    errordlg('Not valid value for stim sampling rate.');
    return;
elseif isnan(NewVal)
    clear global stimsamprate;
    set(handles.stimsamprateD, 'String', '');
    displayspecgram(handles);
    return;
end
% Assign the global variable
global stimsamprate
stimsamprate = NewVal;
set(handles.stimsamprateD, 'String', num2str(stimsamprate));
displayspecgram(handles);

% --- Executes during object creation, after setting all properties.
function stimsamprateD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function respsamprateD_Callback(hObject, eventdata, handles)
newStrValue = get(hObject, 'String');
NewVal = str2double(newStrValue);

% Check if the entered value falls within the allowable range
if  (NewVal< 0)
    % Set default value for step
    warndlg('Not valid value for resp sampling rate.');
    NewVal = 0;
end

% Assign the global variable
global respsamprate ampsamprate
respsamprate = NewVal;
if isempty(ampsamprate)
    ampsamprate = NewVal;
end
set(handles.respsamprateD, 'String', num2str(respsamprate));
%checknload(rawDS, handles.Index,handles);
displayspecgram(handles);

% --- Executes during object creation, after setting all properties.
function respsamprateD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function frame3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


