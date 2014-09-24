function varargout = get_samp_rates(varargin)
% GET_SAMP_RATES Application M-file for songwave_specgram.fig
%    FIG = GET_SAMP_RATES launch songwave_specgram GUI.
%    GET_SAMP_RATES('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 30-Mar-2006 14:54:40
if nargin == 0  % LAUNCH GUI

    fig = openfig(mfilename,'reuse');

    % for resize property
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
        hUIControls],'units','normalized','fontunits','normalized');

    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    handles.index = 1;
    handles.outSpectrum = {};
    guidata(fig, handles);

    % Check whether the input is right for preprocessing
    global rawData rawDS sDim
    if isempty(rawDS)
        errordlg('You need select data files first.',...
            'Selection Error', 'modal')
        return;
    end

    global ampsamprate respsamprate
    if length(ampsamprate) == 0
        ampsamprate = 1000;
    end
    if length(respsamprate) == 0
        respsamprate = 1000;
    end

    set(handles.ampsamprate, 'String', ampsamprate);
    set(handles.respsamprate, 'String', respsamprate);
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
function varargout = ampsamprate_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global ampsamprate 
newStr = get(h, 'String');
NewVal = str2double(newStr);

% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid frequency bandwith (positive number only).',...
        'Variable Error', 'modal')
    return;
end

ampsamprate = NewVal;

% --------------------------------------------------------------------
function varargout = respsamprate_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
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
function varargout = close_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global ampsamprate respsamprate outputPath
save(fullfile(outputPath,'STRFPAK_script_parameters.mat'), 'ampsamprate' ,  'respsamprate' ,'-APPEND');

delete(handles.figure1);



% --- Executes during object creation, after setting all properties.
function upfreq_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');


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


