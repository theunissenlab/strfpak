function varargout = get_filesGUI(varargin)
% GET_FILESGUI Application M-file for get_filesGUI.fig
%    FIG = GET_FILESGUI launch get_filesGUI GUI.
%    GET_FILESGUI('callback_name', ...) invoke the named callback.
%
%             STRFPAK: STRF Estimation Software
% Copyright ?2003. The Regents of the University of California (Regents).
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

% Created by JXZ, 2002.
% Modified by JXZ, 8/23/2005
%  -- remove psth_smoothconst since so far psth_smooth is only for
%  displaying. We only need psth_smoothconst for validation.
% Modified by JXZ, 2/1/2006
%
if nargin == 0  % LAUNCH GUI
    initialize_outputPath;
    if exist(outputPath,'dir')
        initialize = 'done';
        save(fullfile(outputPath,'STRFPAK_script_dataset.mat'),'initialize');
        save(fullfile(outputPath,'STRFPAK_script_parameters.mat'),'initialize');
        global rawDS
        if isempty(rawDS)
            rawDS = guess_rawDS;
        end


        fig = openfig(mfilename,'reuse');
        %set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

        % For resize property
        hAxes = findall(fig,'type','axes');
        hText  = findall(hAxes,'type','text');
        hUIControls = findall(fig,'type','uicontrol');
        set([hAxes; hText;...
            hUIControls],'fontname', 'Times New Roman','units','normalized','fontunits','normalized');

        % Generate a structure of handles to pass to callbacks, and store it.
        handles = guihandles(fig);
        handles.orgpath = pwd;
        % ==============================================
        %  get initial path
        % ==============================================
        initial_dir = pwd;
        handles.orgpath = pwd;
        if exist(fullfile(initial_dir,'DemoData'), 'file')
            initial_dir = fullfile(initial_dir, 'DemoData');
        end

        handles.stimpath = initial_dir;
        handles.resppath = initial_dir;
        handles.stimfiles_selected = {};
        handles.respfiles_selected = {};
        guidata(fig, handles);

        % ==============================================
        %  get initial list of files at the current dir
        % ==============================================
        load_initlistbox(initial_dir, handles);
        cd (handles.orgpath)

        global sDim
        if ~isempty(sDim)
            if strcmp(sDim, '0-D')
                set(handles.spatialdomain, 'value', 2);
            elseif strcmp(sDim, '1-D')
                set(handles.spatialdomain, 'value', 3);
            elseif strcmp(sDim, '2-D')
                set(handles.spatialdomain, 'value', 4);
            end
        end
        global rawData
        if ~isempty(rawData)
            if rawData ==1
                set(handles.rawdata, 'Value', rawData);
                set(handles.preprocessedData, 'Value', 0);
            else
                set(handles.rawdata, 'Value', rawData);
                set(handles.preprocessedData, 'Value', 1);
            end
        end
        if ~isempty(rawDS)
            ResultsStr = {};
            for ii = 1:length(rawDS)
                %ResultsStr = {};
                [p, n, e, r] = fileparts(rawDS{ii}.stimfiles);
                sfilename = [n,e];

                [p, n, e, r] = fileparts(rawDS{ii}.respfiles);
                rfilename = [n, e];
                ResultsStr = [ResultsStr; {[sfilename, '   ',...
                    rfilename]}];

                % reset handles.ResultsData for restarting GET_FILES window
                % modified by Jxz, 05/02/2003
                handles.ResultsData(ii).Stimfile = rawDS{ii}.stimfiles;
                handles.ResultsData(ii).Respfile = rawDS{ii}.respfiles;
            end
            set(handles.list_datapairs, 'String', ResultsStr);
            load_initlistbox(p, handles);
            cd (handles.orgpath)

        end

        if nargout > 0
            varargout{1} = fig;
        end
    end
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    try
        [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    catch
        disp(lasterr);
    end

end


% --------------------------------------------------------------------
function varargout = list_datapairs_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% display list_datapairs

% --------------------------------------------------------------------
function varargout = list_stimfiles_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
currentPath = handles.stimpath;
list_entries = get(handles.list_stimfiles,'String');
index_selected = get(handles.list_stimfiles,'Value');

handles.stimfiles_selected = {};
for ii = 1:length(index_selected)
    filename = list_entries{index_selected(ii)};
    %Check whether the selected file is directory
    newDir = fullfile(currentPath, filename);
    if isdir(newDir)
        %newDir = fullfile(currentPath, filename);
        handles.stimpath = newDir;
        load_stimlistbox(newDir, handles);
    else
        handles.stimfiles_selected{ii} = filename;
    end
end
cd (handles.orgpath)
guidata(h,handles)


% --------------------------------------------------------------------
function varargout = list_respfiles_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
currentPath = handles.resppath;
list_entries = get(handles.list_respfiles,'String');
index_selected = get(handles.list_respfiles,'Value');

handles.respfiles_selected = {};
for ii = 1:length(index_selected)
    filename = list_entries{index_selected(ii)};

    % Check whether the selected file is directory
    newDir = fullfile(currentPath, filename);
    if isdir(newDir)
        %newDir = fullfile(currentPath, filename);
        handles.resppath = newDir;
        load_resplistbox(newDir, handles);
    else
        handles.respfiles_selected{ii} = filename;
    end
end
cd (handles.orgpath)
guidata(h,handles)


% --------------------------------------------------------------------
function varargout = cancel_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
cd (handles.orgpath)
delete(handles.figure1);
clear;


% --------------------------------------------------------------------
function varargout = done_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------

% Check whether data pairs are selected and data type is specified
% before clicking this button.
global rawData rawDS sDim

if isempty(rawDS)
    errordlg('Please select data before exist.', 'Data Missing', 'modal');
    return;
end
if isempty(sDim)
    errordlg('Please specifiy the dimension of stim.','Dimension Missing','modal');
    return;
end
if isempty(rawData)
    errordlg('Please specify the data type (raw or preprocessed).',...
        'Data type missing', 'modal')
    return;
end

temprawData = rawData;
temprawDS = rawDS;
tempsDim = sDim;
global rawData rawDS sDim
rawData = temprawData;
rawDS = temprawDS;
sDim = tempsDim;
if rawData == 0   % User loads already preprocessed data

    %  The following marker %^3 $^3 delimits where STRFPAK should start
    %  cutting and pasting code for the STRFPAK_script.m file.
    %%%$$$ begin preprocess

    global DS preprocessOption
    DS = rawDS;
    preprocessOption = 'Preprocessed data from user';
    % assign global variable for already preprocessed data
    for ii = 1:length(DS)

        % Assign the stimuli data file
        sfile = DS{ii}.stimfiles;

        if 1%ii ==  1  % get value for NumBand
            % Check whether data file type is OK
            [path,name,ext,ver] = fileparts(sfile);
            switch ext
                case {'.dat', '.txt', ''}
                    t = load(sfile);
                    global NBAND
                    NBAND = min(size(t));
                    DS{ii}.nlen = size(t, 2);
                    clear t
                case {'.mat'}
                    smat = load(sfile);
                    flds = fieldnames(smat);
                    stim = getfield(smat, flds{1});
                    global NBAND
                    NBAND = min(size(stim));
                    DS{ii}.nlen = size(stim, 2);
                    clear smat

                otherwise
                    errordlg(['Preprocessed data requirement:Only ASCII and MAT binary fileformats with no extension ',...
                        'or .dat, .txt, and .mat. work now.'], 'File Type Error', 'modal')
                    return;
            end

        end


        % Now Check if RESPONSE FILE is valid
        rfile = DS{ii}.respfiles;
        %real_rfile = [get(handles.resp_filepath, 'String') '/' rfile];

        [path,name,ext,ver] = fileparts(rfile);
        switch ext
            case {'.dat', '.txt', ''}
                t = load(rfile);

                % if stimuli and reponse has different length, we choose min value
                %DS{ii}.nlen = min(size(t, 2), tmpLength);
                DS{ii}.nlen = max(size(t));
                DS{ii}.ntrials = min( size(t));
                clear t;

            case {'.mat'}
                rmat = load(rfile);
                flds = fieldnames(rmat);
                resp = getfield(rmat, flds{1});
                if iscell(resp)  % Check if resp is in spike arrival time format
                    DS{ii}.ntrials = length(resp);
                else
                    DS{ii}.nlen = max(size(resp));
                    DS{ii}.ntrials = min(size(resp));
                end
                clear resp;

            otherwise
                errordlg(['Preprocessed data requirment:Only ASCII and MAT binary fileformats with no extension ',...
                    'or .dat, .txt, and .mat. work now.'], 'File Type Error', 'modal')
                return;
        end

    end
    global originalDS
    originalDS = DS;
    %%%$$$ end preprocess
    add_to_STRFPAK_script('get_filesGUI.m','preprocess');
    global outputPath
    initialize_outputPath;
    save(fullfile(outputPath,'STRFPAK_script_dataset.mat'),'rawDS','outputPath','-APPEND');

else  %  If the stim needs preprocessing, save rawDS instead.
    global outputPath
    initialize_outputPath;
    save(fullfile(outputPath,'STRFPAK_script_dataset.mat'),'rawDS','outputPath','-APPEND');
end
global ampsamprate initialFreq NBAND endFreq respsamprate outputPath
save(fullfile(outputPath,'STRFPAK_script_parameters.mat'),'rawData','sDim','ampsamprate', 'initialFreq', 'NBAND', 'endFreq', '-APPEND');

cd (handles.orgpath)

if rawData == 0
    tt = get_samp_rates;
    uiwait(tt);
end
delete(handles.figure1)
delete(fullfile(outputPath,'*_mean_removed.mat'));

strfFirstGUI;

% --------------------------------------------------------------------
function varargout = stim_samprate_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
%  Callback for stim_samprate_Callback
newStr_stimsamprate = get(h, 'String');
NewVal = str2double(newStr_stimsamprate);

% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for stim sample rate
    errordlg('Please enter valid sampling rate (positive number only).',...
        'Variable Error', 'modal')
    return;
end

% Assign global variable to ampsamprate
global ampsamprate;
ampsamprate = NewVal;


% --------------------------------------------------------------------
function varargout = initialfreq_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
%  Callback for resp_samprate_Callback

newStr_respsamprate = get(h, 'String');
NewVal = str2double(newStr_respsamprate);

% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid sampling rate (positive number only).',...
        'Variable Error', 'modal')
    return;
end


% Assign global variable to ampsamprate
global initialFreq;
initialFreq = NewVal;

% --------------------------------------------------------------------
function varargout = select_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------

% Check whether choose the same number of stimfile and respfile
numSelected = length(handles.stimfiles_selected);
if ( numSelected ~= length(handles.respfiles_selected))
    errordlg(['You must select the same numbers of',...
        ' stimuli and response data.'],...
        'Incorrect Selection', 'modal')
    return;
end

% Retrieve old results data structure
if isfield(handles,'ResultsData') & ~isempty(handles.ResultsData)
    ResultsData = handles.ResultsData;
    hadNum = length(ResultsData);
else % Set up the results data structure
    ResultsData = struct('Stimfile',[],'Respfile',[]);
    hadNum = 0;
end

% Build the new results list string for the list box
ResultsStr = get(handles.list_datapairs, 'String');
if (hadNum == 0)
    ResultsStr = {};
end

% Assign global Variable rawDS
global rawDS
for ii = hadNum+1:hadNum+numSelected

    % Assign the stimuli data file
    sfile = handles.stimfiles_selected{ii-hadNum};
    %real_sfile = [get(handles.stim_filepath, 'String') '/' sfile];
    real_sfile = fullfile(get(handles.stim_filepath, 'String'), sfile);
    ResultsData(ii).Stimfile = sfile;

    rfile = handles.respfiles_selected{ii-hadNum};
    %real_rfile = [get(handles.resp_filepath, 'String') '/' rfile];
    real_rfile = fullfile(get(handles.resp_filepath, 'String'), rfile);

    ResultsData(ii).Respfile = rfile;

    % Assign the global Variable rawDS
    rawDS{ii}.stimfiles = real_sfile;
    rawDS{ii}.respfiles = real_rfile;

    ResultsStr = [ResultsStr; {[ResultsData(ii).Stimfile, '   ',...
        ResultsData(ii).Respfile]}];
end

set(handles.list_datapairs, 'String', ResultsStr);
handles.ResultsData = ResultsData;
guidata(h, handles)


% --------------------------------------------------------------------
function varargout = remove_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
currentVal = get(handles.list_datapairs,'Value');
resultsStr = get(handles.list_datapairs,'String');
numResults = size(resultsStr,1);

% Remove the data and list entry for the selected value
resultsStr(currentVal) =[];
handles.ResultsData(currentVal)=[];

% Update global variable rawDS and NBAND
global rawDS NBAND;
rawDS(currentVal) = [];
%     if ~isempty(rawDS)
%         stim_env = Check_And_Load(rawDS{1}.stimfiles);
%         NBAND = size(stim_env,1);
%     end

% If there are no other entries, disable the Remove and Plot button
% and change the list sting to <empty>
if isequal(numResults,length(currentVal)),
    resultsStr = {'<empty>'};
    currentVal = 1;
end

% Ensure that list box Value is valid, then reset Value and String
currentVal = min(currentVal,size(resultsStr,1));
set(handles.list_datapairs,'Value',currentVal,'String',resultsStr)

% Store the new ResultsData
guidata(h,handles)


% --------------------------------------------------------------------
function varargout = update_stimpath_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------

outputdir = uigetdir(pwd, 'Pick the stimulus directory');

datafiledir = outputdir;
if not(exist(datafiledir,'dir'))
    errordlg('Directory not found, Please try again.','Input Error', 'modal')
    return
end
handles.stimpath = datafiledir;
load_stimlistbox(datafiledir, handles);
cd (handles.orgpath)

% --------------------------------------------------------------------
function varargout = update_resppath_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------

outputdir = uigetdir(pwd, 'Pick the response data directory');

datafiledir = outputdir;
if not(exist(datafiledir,'dir'))
    errordlg('Directory not found. Please try again.','Input Error', 'modal')
    return
end
handles.resppath = datafiledir;
load_resplistbox(datafiledir, handles);
cd (handles.orgpath)

% --------------------------------------------------------------------
function load_initlistbox(dir_path, handles)
% --------------------------------------------------------------------
%  Read the current directory and sort the names
%currentPath = pwd;
cd (dir_path)
dir_struct = dir(dir_path);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.stimfile_names = sorted_names;
handles.stimis_dir = [dir_struct.isdir];
handles.stimsorted_index = [sorted_index];
handles.respfile_names = sorted_names;
handles.respis_dir = [dir_struct.isdir];
handles.respsorted_index = [sorted_index];

guidata(handles.figure1, handles)
set(handles.list_stimfiles,'String',handles.stimfile_names,...
    'Value', 1)
set(handles.stim_filepath, 'String',pwd)
set(handles.list_respfiles,'String',handles.respfile_names,...
    'Value', 1)
set(handles.resp_filepath, 'String',dir_path)
%cd (currentPath)

% --------------------------------------------------------------------
function load_stimlistbox(dir_path, handles)
% --------------------------------------------------------------------
%  Read the current directory and sort the names
%currentPath = pwd;
cd (dir_path)
dir_struct = dir(dir_path);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.stimfile_names = sorted_names;
handles.stimis_dir = [dir_struct.isdir];
handles.stimsorted_index = [sorted_index];

guidata(handles.figure1, handles)
set(handles.list_stimfiles,'String',handles.stimfile_names,...
    'Value', 1)
set(handles.stim_filepath, 'String',pwd)
%cd (currentPath)

% --------------------------------------------------------------------
function load_resplistbox(dir_path, handles)
% --------------------------------------------------------------------
%  Read the current directory and sort the names
%currentPath = pwd;
cd (dir_path)
dir_struct = dir(dir_path);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.respfile_names = sorted_names;
handles.respis_dir = [dir_struct.isdir];
handles.respsorted_index = [sorted_index];

guidata(handles.figure1,handles)

set(handles.list_respfiles,'String',handles.respfile_names,...
    'Value', 1)
set(handles.resp_filepath, 'String',pwd)
%cd (currentPath);


% --------------------------------------------------------------------
function varargout = endfreq_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
newStr_respsamprate = get(h, 'String');
NewVal = str2double(newStr_respsamprate);

% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid sampling rate (positive number only).',...
        'Variable Error', 'modal')
    return;
end


% Assign global variable to ampsamprate
global endFreq;
endFreq = NewVal;


% --------------------------------------------------------------------
function varargout = reset_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
%clear rawDS initialFreq endFreq ampsamprate sDim
clear global rawDS initialFreq endFreq ampsamprate sDim NBAND
handles.ResultsData = [];
guidata(h, handles);
set(handles.initialfreq, 'String', ' ');
set(handles.endfreq, 'String', ' ');
set(handles.stim_samprate, 'String', ' ');
%v = get(handles.dimension, 'String');
set(handles.dimension, 'value',1);
set(handles.list_datapairs, 'String', '<empty>');

% --------------------------------------------------------------------
function varargout = help_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
ttlStr = get(handles.figure1, 'Name');
hlpStr = [...
    '                                                                            '
    '  STRFPAK: Load Input Window                                                '
    'This window lets users to load input files and specify the type of input    '
    'files.                                                                      '
    '                                                                            '
    '   Choose Stimulus:                                                         '
    '        Top panel shows the directory of stimulus or the user can type the  '
    '                   new directory.                                           '
    '        Middle panel lists all the files of the current directory.          '
    '        Bottom panel is the button. When it is clicked, the user can browser'
    '                   the directory.                                           '
    '                                                                            '
    '                                                                            '
    '   Choose Response:                                                         '
    '        Top panel shows the directory of response or the user can type the  '
    '                   new directory.                                           '
    '        Middle panel lists all the files of the current directory.          '
    '        Bottom panel is the button. When it is clicked, the user can browser'
    '                   the directory.                                           '
    '                                                                            '
    '   Select/Remove:                                                           '
    '        To select the data set, the user need click the file names from     '
    '        the middle panels of CHOOSE STIMULUS and CHOOSE REPONSE and then    '
    '        click SELECT button. Similiarly for removing.  For multiple         '
    '        selection, please click Ctrl key and the file names at the same     '
    '        time. Again it is similiar for removing.                            '
    '                                                                            '
    '   Show Selected Data Sets:                                                 '
    '        The panel shows the stimulus name and its corresponding response.   '
    '                                                                            '
    '   Data Type Options:                                                       '
    '        Raw Data: It means the stimulus and the response are in raw format. '
    '        For example: raw auditory stimulus is song wave file with the file  '
    '        ext as .wav and the                                                 '
    '        response data is in spike arrival time or in text format with the   '
    '        number of spikes at the particular time location. For raw vision    '
    '        stimulus, it could be in movie format, X x Y x T. Here X and Y are  '
    '        pixel space position and T is for time.                             '
    '                                                                            '
    '        Preprocessed Data: It means the stimulus is in X x T and the        '
    '        response is in Trial x T.                                           '
    '                                                                            '
    '        HELP - click help button to get this help window.                   '
    '        Cancel - clear all variables and close the window.                  '
    '        DONE - assign the global variables with selected data files and     '
    '         the type of the data sets and close the window.                    '
    '                                                                            '
    ' Updated by JXZ, May 2006 .                                                 '];
myFig = handles.figure1;
helpwin(hlpStr, ttlStr);



% --- Executes on button press in rawdata.
function rawdata_Callback(hObject, eventdata, handles)
global rawData
newVal = get(hObject, 'Value');
if newVal == 1
    rawData = 1;
    set(handles.preprocessedData, 'Value', 0);
end


% --- Executes on button press in preprocessedData.
function preprocessedData_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessedData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of preprocessedData
global rawData
newVal = get(hObject, 'Value');
if newVal == 1
    rawData = 0;
    set(handles.rawdata, 'Value', 0);
end



% --- Executes on button press in rawdatahelp.
function rawdatahelp_Callback(hObject, eventdata, handles)
helpdlg({'Raw data help: ',...
    ' ',...
    '-- Raw data means raw stimulus and raw response data.',...
    '   For example, raw auditory stimulus is .wav file and ',...
    '   raw vision stimulus is movie format. The raw response ',...
    '   data could be raw spike arrival times or numbers of spikes.',...
    ' '},...
    'Raw Data Help');

% --- Executes on button press in ppdatahelp.
function ppdatahelp_Callback(hObject, eventdata, handles)
helpdlg({'Preprocessed data help: ',...
    ' ',...
    '-- Preprocessed data means the stimulus is in the representation appropriate for STRF calculation ',...
    '   and the response is a raster of 0s and 1s.',...
    ' '},...
    'Preprocessed Data Help');

% ============================================================
%  END of get_filesGUI.m
% ============================================================


% --- Executes on selection change in spatialDomain.
function spatialDomain_Callback(hObject, eventdata, handles)
% hObject    handle to spatialDomain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns spatialDomain contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spatialDomain


% --- Executes during object creation, after setting all properties.
function spatialDomain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spatialDomain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in spatialdomain.
function spatialdomain_Callback(hObject, eventdata, handles)
v = get(handles.spatialdomain, 'value');
global sDim;
if v == 1
    errordlg('Please choose stim dimension size!', 'Dim size missing', 'modal')
    clear global sDim
    return;
elseif v == 2
    sDim = '0-D';
    global rawData
    if rawData == 1
        warndlg('For current version of STRFpak, 0-D is no need to preprocess.',...
            '0-D limitation', 'modal');
        rawData = 0;
        set(handles.rawdata, 'Value', rawData);
        set(handles.preprocessedData, 'Value', 1);

    end
elseif v == 3
    sDim = '1-D';
elseif v == 4
    sDim = '2-D';
end


% --- Executes during object creation, after setting all properties.
function spatialdomain_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in spatialdomainhelp.
function spatialdomainhelp_Callback(hObject, eventdata, handles)
helpdlg({'Stim Dimension help: ',...
    ' ',...
    '-- Stim dimension means the number of spatial (not temporal) dimensions the stimulus has.',...
    '   For example, spectrograms are 1-D since ',...
    '   they have one spatial dimension: frequency. The movie file is 2-D ',...
    '   feature dimension because scenes have two spatial dimensions. ',...
    '   If the data is only time series data with no spatial features, ',...
    '   it can be specified as 0-D data. ',...
    ' '},...
    'Stim Dimension Help');


% --- Executes on button press in Auto_select.
function Auto_select_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rawDS
rawDS = guess_rawDS(handles.stimpath,handles.resppath);
if ~isempty(rawDS)
    ResultsStr = {};
    for ii = 1:length(rawDS)
        %ResultsStr = {};
        [p, n, e, r] = fileparts(rawDS{ii}.stimfiles);
        sfilename = [n,e];

        [p, n, e, r] = fileparts(rawDS{ii}.respfiles);
        rfilename = [n, e];
        ResultsStr = [ResultsStr; {[sfilename, '   ',...
            rfilename]}];

        % reset handles.ResultsData for restarting GET_FILES window
        % modified by Jxz, 05/02/2003
        handles.ResultsData(ii).Stimfile = rawDS{ii}.stimfiles;
        handles.ResultsData(ii).Respfile = rawDS{ii}.respfiles;
    end
    set(handles.list_datapairs, 'String', ResultsStr);
    load_initlistbox(p, handles);
    cd (handles.orgpath)

end

%select_button_Callback;
