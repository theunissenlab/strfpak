function varargout = calc_parameterGUI(varargin)
% CALC_PARAMETER2GUI Application M-file for calc_parameterGUI.fig
%    FIG = CALC_PARAMETER2GUI launch calc_parameterGUI GUI.
%    CALC_PARAMETER2GUI('callback_name', ...) invoke the named callback.
%
%             STRFPAK: STRF Estimation Software
% Copyright ?003. The Regents of the University of California (Regents).
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

% Created by JXZ, 2003.
% 
% Modified by JXZ, 7/12/2005
%    1. Add global variable Std_val
%    2. Add button Std_val and editable text field stdval_list
%    

if nargin == 0  % LAUNCH GUI

    fig = openfig(mfilename,'reuse');
    %set(fig,'Color',[0.7 0.7 0.7]);

    % for resize property
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
           hUIControls],'fontname', 'Times New Roman','units','normalized','fontunits','normalized');

    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    guidata(fig, handles);

    % display if TimeLag and Tol_val are already set
    global TimeLag
    if ~isempty(TimeLag)
        set(handles.twindow, 'String', TimeLag);
    end
    
    global Tol_val
    if ~isempty(Tol_val)
        set(handles.tol_startvalue, 'String', num2str(Tol_val, 3));
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
    
   
    % save handles to figure 
    guidata(fig, handles);
    
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


% --------------------------------------------------------------------
function varargout = twindow_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    % Get twindow values from fig
    newStrValue = get(h, 'String');
    NewVal = str2double(newStrValue);

    % Check if the entered value falls within the allowable range
    if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
        errordlg('Please enter valid Time-lag value.(positive number only)',...
            'Global Variable Error', 'modal');
        return;
    end

    handles.settwindow = 1;
    guidata(h, handles);

    % Assign the global variable TimeLag
    global TimeLag
    TimeLag = NewVal;


% --------------------------------------------------------------------
function varargout = tol_startvalue_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    InputStr = get(h, 'String');

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
    
    Tol_val = [tolval];

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
	    Tol_val = [Tol_val tolval];
    end

    handles.setstart = 1;
    guidata(h, handles);


% --------------------------------------------------------------------
function varargout = tol_stopvalue_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    newStrValue = get(h, 'String');
    NewVal = str2double(newStrValue);

    % Check if the entered value falls within the allowable range
    if  isempty(NewVal) | (NewVal< 0) 
        % Set default value for stop 
        NewVal = 1;
    end

    handles.setstop = 1;
    guidata(h, handles);

    % Assign the global variable TimeLag
    global Tol_stop
    Tol_stop = NewVal;

% --------------------------------------------------------------------
function varargout = tol_stepvalue_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    newStrValue = get(h, 'String');
    NewVal = str2double(newStrValue);

    % Check if the entered value falls within the allowable range
    if  isempty(NewVal) | (NewVal< 0) 
        % Set default value for step
        NewVal = 1;
    end

    handles.setstep = 1;
    guidata(h, handles);

    % Assign the global variable tol_step
    global Tol_step
    Tol_step = NewVal;

% --------------------------------------------------------------------
function varargout = parameter_help_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    ttlStr = get(handles.figure1, 'Name');
    hlpStr = [...
             '                                                                       '
             '  STRFPAK: Calculatin Parameters                                       '
             '                                                                       '
             'The CALC_PARAMETER Window is needed before doing calculation. It helps '
             'to set all calculation parameters. STRFPAK also provides a brief       '
             'explaination for each parameters after the parameter button is clicked.'
             '                                                                       '
             'These parameters are:                                                  '
             '                                                                       '
             '     TOL_VALUE  -- a regularization hyperparameter used to estimate the'
             ' STRF. Larger tol values yield smooth STRFs. Lower tol values remove   '
             'more of the stimulus correlation but can amplify noise. Use a range of '
             'them to find the best STRF for your specific stimulus and neuron.      '
             'You can separate them by ",", "space" or "tab". For example:           '
             'Tol_Val = [0.1 0.05 0.001 0.0005].                                     '
             '                                                                       '
             '     Time-lag -- is the time lag used for calculating the spike        '
             'triggered average and the stimulus auto-correlation. Time-lag should   '
             ' be long enough to capture the memory of the neuron and the correlation'
             '  time of the stimulus. The unit of Time-lag is ms.                    '
             '                                                                       '
             '    Std. Value -- a second hyperparameter in the STRF estimation that  '
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
             'Updated by JXZ, Sept. 2005.                                            '
             ];
    myFig = handles.figure1;
    helpwin(hlpStr, ttlStr);          

% --------------------------------------------------------------------
function varargout = parameter_close_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    delete(handles.figure1);
    clear;


% --------------------------------------------------------------------
function varargout = done_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    % Retrieve old results data structure
    global TimeLag Tol_val timevary_PSTH smooth_rt DS Std_val
    if isempty(DS)
        errordlg('You need choose input data first.',...
            'Incorrect Selection', 'modal')
        return
    elseif length(DS) == 1
        warndlg('Since only one data set is chosen, Std_val has to be zero.');
        Std_val = [0];
    end
    if ~isempty(TimeLag) & ~isempty(Tol_val) & ~isempty(timevary_PSTH)& ~isempty(smooth_rt)
        delete(handles.figure1);
        clear;
    else % Set up the results data structure
        errordlg('You need set all the parameters.',...
            'Incorrect Selection', 'modal')
        return
    end

    % 11/11/2005
    % Remove all the existing results
    global outputPath
    outdatafiledir = outputPath;
    if ~isempty(outputPath)
         % Junli: 11/11/2005
          % existing directory
            tt = dir(outdatafiledir);
            if length(tt) >2    % more file there
                anw = warndlg(['You modified calc parameters. You will override the calculation results.'],...
                    'Warning prediction', 'modal');
                %uiwait(anw);
                current_dir = pwd;
                cd(outdatafiledir)
                delete Global_*.mat;
                delete stim_avg.mat;
                delete display_*.mat;
                delete info_r*.mat;
                delete pred*.mat;
                delete strfResult*.mat;
                delete strfH*.mat;
                cd(current_dir);
            end
    end
                       
               

% --------------------------------------------------------------------
function varargout = nband_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
if 0
    newStrValue = get(h, 'String');
    NewVal = str2double(newStrValue);

    % Check that the entered value falls within the allowable range
    if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
        errordlg({'Please enter valid spatio_size value.(positive ',...
                 'number only)'}, 'Input Error', 'modal')
        return;
    end

    handles.setnband = 1;
    guidata(h, handles);

    % Assign the global variable NBAND
    global NBAND 
    NBAND = NewVal;
end
    

% --------------------------------------------------------------------
function varargout = twindowhelp_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    % TWINOW HELP
    helpdlg({'Time-lag: ',...
      ' ',...
      '-- is the time lag used for calculating the spike triggered',...
      ' average and the stimulus auto-correlation.  Time-lag should',...
      ' be long enough to capture the memory of the neuron and the', ...
      ' correlation time of the stimulus. The unit of Time-lag is ms.',...
      ' '},...
      'Time-lag Help');
     
% --------------------------------------------------------------------
function varargout = tolvalhelp_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    % Tol Value Help
    helpdlg({'Tol. Values: ',...
            ' ',...
            ' -- a regularization hyperparameter used to estimate the STRF.',...
            ' Larger tol values yield smooth STRFs. Lower tol values remove',...
            ' more of the stimulus correlation but can amplify noise.', ...
            ' Use a range of them to find the best STRF for your specific',...
            ' stimulus and neuron. You can separate them by ",", "space"',...
            ' or "tab". For example: Tol_Val = [0.1 0.05 0.001 0.0005].',...
            ' '}, 'Tol Value Help');
    




% --------------------------------------------------------------------
function varargout = yesRadio_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.yesRadio.

    global timevary_PSTH
    newVal = get(h, 'Value');
    if newVal == 1
        timevary_PSTH = 1;
        set(handles.noRadio, 'Value', 0);
    end


% --------------------------------------------------------------------
function varargout = noRadio_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.noRadio.

    global timevary_PSTH
    newVal = get(h, 'Value');
    if newVal == 1
        timevary_PSTH = 0;
        set(handles.yesRadio, 'Value', 0);
    end


% --- Executes on button press in stdvals.
function stdvals_Callback(hObject, eventdata, handles)
% hObject    handle to stdvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
     helpdlg({'Std Values: ',...
             ' ',...
            ' -- a second hyperparameter in the STRF estimation that enforces',...
            ' sparseness. High Std values limits the number of significant',...
            ' pixels in the STRF. Use a range of them to find the best STRF',...
            ' for your specific stimulus and neuron. You can separate them',...
            ' by ",", "space" or "tab". For example: ',...
            ' Std_Val = [0 0.1 2 4].',...
            ' '}, 'Tol Value Help');

% --- Executes during object creation, after setting all properties.
function stdval_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stdval_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

    set(hObject,'BackgroundColor','white');


function stdval_list_Callback(hObject, eventdata, handles)
% hObject    handle to stdval_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stdval_list as text
%        str2double(get(hObject,'String')) returns contents of stdval_list as a double

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

    if length(DS) == 1
        Std_val = [0];
        set(handles.stdval_list, 'String', num2str(Std_val, 3));
    end
    guidata(hObject, handles);


% --- Executes on button press in modeldes.
function modeldes_Callback(hObject, eventdata, handles)
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


% --- Executes during object creation, after setting all properties.
function smoothmeanrate_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor','white');



function smoothmeanrate_Callback(hObject, eventdata, handles)
    %  Get smooth_rt from fig
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
% =====================================================================
%  END of calc_parameterGUI.m
% =====================================================================


% --- Executes on button press in calalg.
function calalg_Callback(hObject, eventdata, handles)
% hObject    handle to calalg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

helpdlg({'Calculation Algorithm: ',...
    ' ',...
    '    STRFPAK provides two options for estimating the STRF: ',...
    '    space-time separable algorithm (fast) and ',...
    '    space-time nonseparable algorithm (more general).     ',...
    '                                                          ',...
    '                                                          '},...
    'Calculation Algorithm Help');

% --- Executes on button press in separable.
function separable_Callback(hObject, eventdata, handles)
% hObject    handle to separable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of separable

 global setSep
    newVal = get(hObject, 'Value');
    if newVal == 1
        setSep = 1;
        set(handles.nonseparable, 'Value', 0);
    end

% --- Executes on button press in nonseparable.
function nonseparable_Callback(hObject, eventdata, handles)
% hObject    handle to nonseparable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nonseparable

  global setSep
    newVal = get(hObject, 'Value');
    if newVal == 1
        setSep = 0;
        set(handles.separable, 'Value', 0);
    end

