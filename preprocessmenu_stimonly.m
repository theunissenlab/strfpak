function varargout = preprocessmenu_stimonly(varargin)
% PREPROCESSMENU_STIMONLY Application M-file for preprocessmenu_stimonly.fig
%    FIG = PREPROCESSMENU_STIMONLY launch preprocessmenu_stimonly GUI.
%    PREPROCESSMENU_STIMONLY('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 03-Mar-2006 10:57:27

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
	%set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end
    
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
           hUIControls],'fontname', 'Times New Roman','units','normalized','fontunits','normalized');

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end


% --------------------------------------------------------------------
function varargout = wavespectrum_Callback(h, eventdata, handles, varargin)
   % call songwave_spectrum to transform real data
global preprocessOption

tt=songwave_spectrum_stimonly;
uiwait(tt);
preprocessOption = 'SongWave->Spectrogram(STFFT)';
delete(handles.figure1);
strfFirstGUI;

% --------------------------------------------------------------------
function varargout = wavelet1d_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global preprocessOption

tt=wavelet1d_scalogram_stimonly;
uiwait(tt);
preprocessOption = 'SongWave->Scalogram(wavelet)';
delete(handles.figure1);
strfFirstGUI;
% --------------------------------------------------------------------
function varargout = pfft2d_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
   
global preprocessOption

tt=movie_pfft_stimonly;
uiwait(tt);
preprocessOption = 'Fourier Power transform';
delete(handles.figure1);
strfFirstGUI;

% --------------------------------------------------------------------
function varargout = spacetime_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global preprocessOption

tt=movie_spacetime_stimonly;
uiwait(tt);
preprocessOption = 'Space-Time transform';
delete(handles.figure1);
strfFirstGUI;

% --------------------------------------------------------------------
function varargout = lyon_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
   msgbox('Coming soon...', 'modal');

% --------------------------------------------------------------------
function varargout = psfft2D_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
   msgbox('Coming soon...', 'modal');

% --------------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    % Close button function
    delete(handles.figure1)
    
% --- Executes on button press in pp_help.
function pp_help_Callback(hObject, eventdata, handles)
ttlStr = get(handles.figure1, 'Name');
hlpStr = [...
    '     STRFPAK: Input Preprocessing Menu (Help window)              '
    '                                                                  '
    '   STRFPAK provides different options for preprocessing           '
    ' time-series data. 1-D option is for audition data; 2-D option is '
    ' for vision data.                                                 '
    ' For 1-D option:                                                  '
    '   1. Short-time Fourier Transform:                               '
    '      First create a bank of band-passed filters whose impulse    '
    '      response consists of a carrier frequency (the center        '
    '      frequency of the filter) modulated by a short-time window.  '
    '      The Fourier transform of this time window is the frequency  '
    '      gain curve of the band-passed filter.                       '
    '   2. Wavelet Transform:                                          '
    '      This transform uses the Morelet wavelet filter banks to get '
    '      approximately spaced logarithmically in Frequency domain.   '   
    '                                                                  '
    '   3. Lyon Model:    Lyon model includes not only the approximate '
    '      logarithmic spacing of filter center frequencies (log at    '
    '      high frequencies and more linear at low frequencies), but   '
    '      also the rectification performed by inner hair cells and,   '
    '      optionally, adaptive gain control.                          '
    '                                                                  '
    ' Here is the reference for the above methods:                     '
    '                                                                  '
    '      Patrick Gill, Junli Zhang, Sarah M. N. Woolley,             '
    '      Thane Fremouw, Frederic E. Theunissen (2005). Analysis of   '
    '      Signal Preprocessing Methods in Spectro-temporal Receptive  '
    '      Field Estimation, Journal of Comp.                          '
    '                                                                  '
    ' For 2-D options, STRFPAK has the following choices:              '
    '                                                                  '
    '    1. Linear Model:                                              '
    '    2. Fourier Power Transform:                                   '
    '    3. Phase-separated Fourier Power Transform:                   '
    '                                                                  '
    ' For the reference:                                               '
    '  1. David SV, Vinje WE and Gallant JL (2004). Natural sitmulus   '
    '     statistics alter the receptive field  structure of V1 neurons'
    '     Journal of Neuroscience, 24, 6991-7006.                      '
    '  2. David, SV, Gallant JL (2005). Predicting neural responses    '
    '     during natural vision, Network: Comp. in Neural Systems,     '
    '     16, 239-260.                                                 '
    '                                                                  '];
    myFig = handles.figure1;
    helpwin(hlpStr, ttlStr);
    


% --- Executes on button press in centerpp.
function centerpp_Callback(hObject, eventdata, handles)
% hObject    handle to centerpp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Coming soon...', 'modal');

% --- Executes on button press in spherizepp.
function spherizepp_Callback(hObject, eventdata, handles)
% hObject    handle to spherizepp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Coming soon...', 'modal');

% --- Executes on button press in downsamplepp.
function downsamplepp_Callback(hObject, eventdata, handles)
% hObject    handle to downsamplepp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Coming soon...', 'modal');

% --- Executes on button press in croppp.
function croppp_Callback(hObject, eventdata, handles)
% hObject    handle to croppp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Coming soon...', 'modal');

% --- Executes on button press in fft2pp.
function fft2pp_Callback(hObject, eventdata, handles)
% hObject    handle to fft2pp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Coming soon...', 'modal');

% --- Executes on button press in pcapp.
function pcapp_Callback(hObject, eventdata, handles)
% hObject    handle to pcapp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Coming soon...', 'modal');

% --- Executes on button press in ownpp.
function ownpp_Callback(hObject, eventdata, handles)
% hObject    handle to ownpp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Coming soon...', 'modal');

% ====================================================================
%   END of preprocessmenu_stimonly.m
% ====================================================================
 
