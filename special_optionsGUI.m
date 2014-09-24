function varargout = special_optionsGUI(varargin)
% SPECIAL_OPTIONSGUI M-file for special_optionsGUI.fig
%      SPECIAL_OPTIONSGUI, by itself, creates a new SPECIAL_OPTIONSGUI or raises the existing
%      singleton*.
%
%      H = SPECIAL_OPTIONSGUI returns the handle to a new SPECIAL_OPTIONSGUI or the handle to
%      the existing singleton*.
%
%      SPECIAL_OPTIONSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECIAL_OPTIONSGUI.M with the given input arguments.
%
%      SPECIAL_OPTIONSGUI('Property','Value',...) creates a new SPECIAL_OPTIONSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before special_optionsGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to special_optionsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help special_optionsGUI

% Last Modified by GUIDE v2.5 18-Jan-2007 12:02:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @special_optionsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @special_optionsGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

if nargin == 0
    fig = openfig(mfilename,'reuse');
    handles = guihandles(fig);
    empty_file_string = 'Click here to choose subspace file.';
    global alien_space_file allow_negative_rates use_alien_space

    if isempty(alien_space_file)
        set(handles.choose_subspace_file,'String',empty_file_string);
    else
        set(handles.choose_subspace_file,'String',alien_space_file);
    end

    if allow_negative_rates == 1
        set(handles.allow_negative_rates,'Value',1);
    end

    if use_alien_space == 1
        set(handles.alien_space,'Value',1);
    end

end
% --- Executes just before special_optionsGUI is made visible.
function special_optionsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to special_optionsGUI (see VARARGIN)

% Choose default command line output for special_optionsGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes special_optionsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = special_optionsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in allow_negative_rates.
function allow_negative_rates_Callback(hObject, eventdata, handles)
% hObject    handle to allow_negative_rates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global allow_negative_rates outputPath
allow_negative_rates = get(hObject,'Value');
save(fullfile(outputPath,'STRFPAK_script_parameters.mat'),'allow_negative_rates','-APPEND');
% Hint: get(hObject,'Value') returns toggle state of allow_negative_rates


% --- Executes on button press in alien_space.
function alien_space_Callback(hObject, eventdata, handles)
% hObject    handle to alien_space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global use_alien_space alien_space_file outputPath
use_alien_space = get(hObject,'Value');
save(fullfile(outputPath,'STRFPAK_script_parameters.mat'),'use_alien_space','-APPEND');
if isempty(alien_space_file) & use_alien_space
    choose_subspace_file_Callback;
end
% Hint: get(hObject,'Value') returns toggle state of alien_space


% --- Executes on button press in choose_subspace_file.
function choose_subspace_file_Callback(hObject, eventdata, handles)
% hObject    handle to choose_subspace_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%disp('I''m here!!!');
global alien_space_file use_alien_space
if isempty(alien_space_file)
    [fname,fpath] = uigetfile('*.mat','Select subspace file (might be called "subspace.mat")');
else
    if exist(alien_space_file,'file')
        [fname,fpath] = uigetfile('*.mat','Select subspace file (might be called "subspace.mat")',alien_space_file);
    else
        [fname,fpath] = uigetfile('*.mat','Select subspace file (might be called "subspace.mat")');
    end
end
if ~fname
    msgbox(['Trouble getting the subspace file.  If you''ve changed your mind, OK.  ' char(10) ...
        'If you''re having trouble with the "uigetfile" dialog on Ubuntu, you''re not alone;' char(10) ...
        'try typing in the full path of the subspace file you want to use.']);

    return
end

trial_alien_space_file = fullfile(fpath,fname);


cached_dir = dir_of_caches;
loaded_alien_pointer = load(trial_alien_space_file);
if ~isfield(loaded_alien_pointer,'original_subspace_checksum') | ~isfield(loaded_alien_pointer,'original_subspace_tol_vals')
    msg = ['Error: the file ' fname ' is not a valid subspace pointer file.'];
    tt = errordlg(msg);
    uiwait(tt);
    return;
else
    if ~exist(fullfile(cached_dir,[loaded_alien_pointer.original_subspace_checksum,'.mat']))
        msg = ['Error: the file ' fname ' is well-formed, but the cached file it points to is missing.']
        errordlg(msg);
        return;
        if isfield(loaded_alien_pointer,'original_time_subspace_checksum') % Used only for separable STRFs
            if ~exist(fullfile(cached_dir,[loaded_alien_pointer.original_time_subspace_checksum,'.mat']))
                msg = ['Error: the file ' fname ' is well-formed, but the cached file it points to is missing.']
                errordlg(msg);
                return;
            end
        end
    else
        use_alien_space = 1;
        alien_space_file = trial_alien_space_file;
    end
end

global Tol_val
Tol_val = loaded_alien_pointer.original_subspace_tol_vals;
special_optionsGUI;

% --- Executes on button press in finished_special_options.
function finished_special_options_Callback(hObject, eventdata, handles)
% hObject    handle to finished_special_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = openfig(mfilename,'reuse');
handles = guihandles(fig);
delete(handles.figure1);
return
% --- Executes on button press in Show_neg_rates_help.
function Show_neg_rates_help_Callback(hObject, eventdata, handles)
% hObject    handle to Show_neg_rates_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ttlStr = get(handles.figure1, 'Name');
fid = fopen('negative_firing_rate_help.txt','r');
tdump = char(fread(fid,'char')');
tdump = tdump(4:end);  %Removes strange charachters at the beginning.
fclose(fid);
helpwin(tdump, ttlStr);


% --- Executes on button press in help_on_subspacecs.
function help_on_subspacecs_Callback(hObject, eventdata, handles)
% hObject    handle to help_on_subspacecs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ttlStr = get(handles.figure1, 'Name');
fid = fopen('alien_subspace_help.txt','r');
tdump = char(fread(fid,'char')');
tdump = tdump(4:end);  %Removes strange charachters at the beginning.
fclose(fid);

helpwin(tdump, ttlStr);

