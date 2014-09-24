function varargout = get_predfilesGUI(varargin)
% GET_PREDFILESGUI Application M-file for get_predfilesGUI.fig
%    FIG = GET_PREDFILESGUI launch get_predfilesGUI GUI.
%    GET_PREDFILESGUI('callback_name', ...) invoke the named callback.
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

    fig = openfig(mfilename);
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
    
    global rawpredData
    if ~isempty(rawpredData)
        if rawpredData ==1
            set(handles.rawdata, 'Value', rawpredData);
            set(handles.preprocessedData, 'Value', 0);
        else
            set(handles.rawdata, 'Value', rawpredData);
            set(handles.preprocessedData, 'Value', 1);
        end
    end
    global rawpredDS
    if ~isempty(rawpredDS)
        ResultsStr = {};
        for ii = 1:length(rawpredDS)
            %ResultsStr = {};
            [p, n, e, r] = fileparts(rawpredDS{ii}.stimfiles);
            sfilename = [n,e];
            
            [p, n, e, r] = fileparts(rawpredDS{ii}.respfiles);
            rfilename = [n, e];
            ResultsStr = [ResultsStr; {[sfilename, '   ',...
                    rfilename]}];
    
            % reset handles.ResultsData for restarting GET_FILES window
            % modified by Jxz, 05/02/2003
            handles.ResultsData(ii).Stimfile = rawpredDS{ii}.stimfiles;
            handles.ResultsData(ii).Respfile = rawpredDS{ii}.respfiles;
        end
        set(handles.list_datapairs, 'String', ResultsStr);
        load_initlistbox(p, handles);
        cd (handles.orgpath)
        
    end
    
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
global rawpredData rawpredDS
if isempty(rawpredData) || isempty(rawpredDS)
    errordlg(['Please select the data sets AND assign ',...
        'data type ( raw or preprocessed) '],...
        'Incorrect Selection', 'modal')
    return;
end
temprawpredData = rawpredData;
temprawpredDS = rawpredDS;
clear global *pred*
global rawpredData rawpredDS
rawpredData = temprawpredData;
rawpredDS = temprawpredDS;
if rawpredData == 0   % User loads already preprocessed data
    global predDS preprocessOption_pred
    predDS = rawpredDS;
    preprocessOption_pred = 'Preprocessed data from user';
    % assign global variable for already preprocessed data
    for ii = 1:length(predDS)

        % Assign the stimuli data file
        sfile = predDS{ii}.stimfiles;

        if ii ==  1  % get value for NumBand
            % Check whether data file type is OK
            [path,name,ext,ver] = fileparts(sfile);
            switch ext
                case {'.dat', '.txt', ''}
                    t = load(sfile);
                    global NBAND
                    NBAND = min(size(t));
                    %tmpLength = size(t, 2);
                    clear t
                case {'.mat'}
                    smat = load(sfile);
                    flds = fieldnames(smat);
                    stim = getfield(smat, flds{1});
                    global NBAND
                    NBAND = min(size(stim));
                    %tmpLength = size(stim, 2);
                    clear smat stim

                otherwise
                    errordlg(['Preprocessed data requirement:Only ASCII and MAT binary fileformats with no extension ',...
                        'or .dat, .txt, and .mat. work now.'], 'File Type Error', 'modal')
                    return;
            end

        end

        % Now Check if RESPONSE FILE is valid
        rfile = predDS{ii}.respfiles;
        %real_rfile = [get(handles.resp_filepath, 'String') '/' rfile];
        
        [path,name,ext,ver] = fileparts(rfile);
        switch ext
            case {'.dat', '.txt', ''}
                t = load(rfile);

                % if stimuli and reponse has different length, we choose min value
                %predDS{ii}.nlen = min(size(t, 2), tmpLength);
                predDS{ii}.nlen = max(size(t));
                predDS{ii}.ntrials = min( size(t));
                clear t;

            case {'.mat'}
                rmat = load(rfile);
                flds = fieldnames(rmat);
                resp = getfield(rmat, flds{1});
                predDS{ii}.nlen = max(size(resp));
                predDS{ii}.ntrials = min(size(resp));
                clear resp;

            otherwise
                errordlg(['Preprocessed data requirment:Only ASCII and MAT binary fileformats with no extension ',...
                    'or .dat, .txt, and .mat. work now.'], 'File Type Error', 'modal')
                return;
        end

    end
end

cd (handles.orgpath)
delete(handles.figure1)
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
function varargout = dimension_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    v = get(handles.dimension, 'value');
    option = get(handles.dimension, 'String');
    option = deblank(option(v,:));
    
    global sDim;
    if v == 1
        errordlg('Invalid spatio-dimension', 'Incorrect Selection', 'modal')
        return;
    elseif v == 2
        sDim = '1-D';
    elseif v == 3
        sDim = '2-D';
    else
        errordlg('This option is still under construction.');
        set(handles.dimension, 'value', 1);
        return;
    end
    

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
             
    % Assign global Variable rawpredDS
    global rawpredDS
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

        % Assign the global Variable rawpredDS
        rawpredDS{ii}.stimfiles = real_sfile;
        rawpredDS{ii}.respfiles = real_rfile;

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
    
    % Update global variable rawpredDS and NBAND
    global rawpredDS NBAND;
    rawpredDS(currentVal) = [];
%     if ~isempty(rawpredDS)
%         stim_env = Check_And_Load(rawpredDS{1}.stimfiles);
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
    %clear rawpredDS initialFreq endFreq ampsamprate sDim
    clear global rawpredDS initialFreq endFreq ampsamprate sDim NBAND outputPath
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
            '  STRFPAK: Load Validated Input Window                                      '
            'This window lets users to input data files and any related information.     '
            'If the preprocessing procedure is used through STRFPAK, this window is not  '
            'necessary to fill out since it will be automatically filled. However, if    '
            'the data has already been preprocessed somewhere else or sometime earlier,  '
            'the preprocssed data can be input from here.                                '
            'All related information includes:                                           '
            '                                                                            '
            '        STIM_FILES       - one or multiple stimulus files                   '
            '        RESP_FILES       - one or multiple response files                   '
            '        DataSAMPRATE      - preprocessed stimulus/response sampling rate    '
            '        SpatioDimension  - dimension size of data spatio domain.            '
            '                                                                            '
            'Where the requirement for selected data files are:                          '
            '        STIM_FILES must be two-dimensional: Spatio X Time. If the spatio    '
            '        domain has 2/more dimensions, please wrapper them into big vectors. '
            '        RESP_FILES must be two-dimensional: Trials X Time (if PSTH file is  '
            '            in 1 X TIME, it is also OK.)  The data fileformat must be ascii '
            '            or MAT format. They have ext like: .txt, .dat, .mat.            '
            '        DataSamping_Rate is in Hz.                                          '
            '        SpationDimension is the number of dimension in spatio domain.       '
            '               It the data has 2-D spatio domain, displaying will reshape   '
            '               the data into X x Y x Time.                                  '
            '                                                                            '
            'The window also provide the following actions:                              '
            '                                                                            '
            '        Browse_STIM_PATH - click the button to update stim data path        '
            '        Browse_RESP_PATH - click the button to update resp data path        '
            '        SELECT - select the data pair of selected stim file and resp files  '
            '        To select multiple files, you hold Ctr and press left mouse button. '
            '        The selcted data pairs appear in the Stim Response list box.        '
            '        You can remove it if you select the pairs from SRlistbox.           '
            '        HELP - click help button to get this help window.                   '
            '        Cancel - clear all variables and close the window.                  '
            '        DONE - assign the global variables with selected data files and     '
            '         sampling rates and dimension selction and close the window.        '
            '                                                                            '
            ' Updated by JXZ, Sept. 2005.                                                '];
    myFig = handles.figure1;
    helpwin(hlpStr, ttlStr);



% --- Executes on button press in rawdata.
function rawdata_Callback(hObject, eventdata, handles)
% hObject    handle to rawdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rawdata
 global rawpredData
    newVal = get(hObject, 'Value');
    if newVal == 1
        rawpredData = 1;
        set(handles.preprocessedData, 'Value', 0);
    end


% --- Executes on button press in preprocessedData.
function preprocessedData_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessedData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of preprocessedData
global rawpredData
    newVal = get(hObject, 'Value');
    if newVal == 1
        rawpredData = 0;
        set(handles.rawdata, 'Value', 0);
    end

% ============================================================
%  END of get_predfilesGUI.m
% ============================================================

