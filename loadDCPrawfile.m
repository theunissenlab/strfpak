function varargout = loaddcprawfile(varargin)
% LOADDCPRAWFILE Application M-file for loaddcprawfile.fig
%    FIG = LOADDCPRAWFILE launch loaddcprawfile GUI.
%    LOADDCPRAWFILE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 20-Aug-2003 14:17:04
%             STRFPAK: STRF Estimation Software
% Copyright ©2003. The Regents of the University of California (Regents).
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
%


if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	%set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
           hUIControls],'units','normalized','fontunits','normalized');


	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    % ==============================================
    %  get initial path
    % ==============================================
    initial_dir = pwd;
    handles.orgpath = pwd;
    if exist(fullfile(initial_dir,'DemoData'), 'file')
        initial_dir = fullfile(initial_dir, 'DemoData');
    end
    handles.filepath = initial_dir;
    handles.psthfilepath = initial_dir;
    guidata(fig, handles);

    % ==============================================
    %  get initial list of files at the current dir
    % ==============================================
    load_initlistbox(initial_dir, handles);
    cd (handles.orgpath)
    
    global rawDS
    if ~isempty(rawDS)
        ResultsStr = {};
        for ii = 1:length(rawDS)
            %ResultsStr = {};
            [p, n, e, r] = fileparts(rawDS{ii}.stimfiles);
            sfilename = [n,e];
            
            ResultsStr = [ResultsStr; {sfilename}];
    
            % reset handles.ResultsData for restarting GET_FILES window
            % modified by Jxz, 05/02/2003
            handles.ResultsData(ii).Stimfile = rawDS{ii}.stimfiles;
            
        end
     
        set(handles.selectstim, 'String', ResultsStr);
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


% ====================================================================
%   Callback functions
% ====================================================================
% --------------------------------------------------------------------
function varargout = stimlist_Callback(h, eventdata, handles, varargin)
    % list file under current file directory
    currentPath = handles.filepath;
    list_entries = get(handles.stimlist,'String');
    index_selected = get(handles.stimlist,'Value');

    handles.stimfiles_selected = {};
    for ii = 1:length(index_selected)
        filename = list_entries{index_selected(ii)};
        %Check whether the selected file is directory
        newDir = fullfile(currentPath, filename);
        if isdir(newDir)
            %newDir = fullfile(currentPath, filename);
            handles.filepath = newDir;
            load_stimlistbox(newDir, handles);
        end
        handles.stimfiles_selected{ii} = filename;
    end
    cd (handles.orgpath)
    guidata(h,handles)
    
% --------------------------------------------------------------------
function varargout = resplist_Callback(h, eventdata, handles, varargin)
   % list file under current file directory
    currentPath = handles.psthfilepath;
    list_entries = get(handles.resplist,'String');
    index_selected = get(handles.resplist,'Value');

    handles.respfiles_selected = {};
    for ii = 1:length(index_selected)
        filename = list_entries{index_selected(ii)};
        %Check whether the selected file is directory
        newDir = fullfile(currentPath, filename);
        if isdir(newDir)
            %newDir = fullfile(currentPath, filename);
            handles.psthfilepath = newDir;
            load_resplistbox(newDir, handles);
        end
        handles.respfiles_selected{ii} = filename;
    end
    cd (handles.orgpath)
    guidata(h,handles)

% --------------------------------------------------------------------
function varargout = selectstim_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = select_Callback(h, eventdata, handles, varargin)
    % Check whether choose the same number of stimfile and respfile
    numSelected = length(handles.stimfiles_selected);
    
    % Retrieve old results data structure
    if isfield(handles,'ResultsData') & ~isempty(handles.ResultsData) 
        ResultsData = handles.ResultsData;
        hadNum = length(ResultsData);
    else % Set up the results data structure
        ResultsData = struct('Stimfile',[],'Respfile',[]);
        hadNum = 0;
    end
    
    % Build the new results list string for the list box
    ResultsStr = get(handles.selectstim, 'String');
    if (hadNum == 0)
        ResultsStr = {};
    end
             
    % Assign global Variable DS
    global rawDS;
    for ii = hadNum+1:hadNum+numSelected
       
        % Assign the stimuli data file 
        sfile = handles.stimfiles_selected{ii-hadNum};
        %real_sfile = [get(handles.stim_filepath, 'String') '/' sfile];
        real_sfile = fullfile(get(handles.stimpath, 'String'), sfile);
         
        ResultsData(ii).Stimfile = sfile;
        
        % Assign the response data file 
        rfile = handles.respfiles_selected{ii-hadNum};
       
        real_rfile = fullfile(get(handles.resppath, 'String'), rfile);
         
        ResultsData(ii).Stimfile = sfile;
        ResultsData(ii).Respfile = rfile;
        
        % Assign the global Variable DS
        rawDS{ii}.stimfiles = real_sfile;
        rawDS{ii}.respfiles = real_rfile;
        
        ResultsStr = [ResultsStr; {[ResultsData(ii).Stimfile, '       ',...
                    ResultsData(ii).Respfile]}];
    end
    
    set(handles.selectstim, 'String', ResultsStr);
    handles.ResultsData = ResultsData;
    guidata(h, handles)
    


% --------------------------------------------------------------------
function varargout = remove_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    currentVal = get(handles.selectstim,'Value');
    resultsStr = get(handles.selectstim,'String');
    numResults = size(resultsStr,1);

    % Remove the data and list entry for the selected value
    resultsStr(currentVal) =[];
    handles.ResultsData(currentVal)=[];
    
    % Update global variable rawDS and NBAND
    global rawDS
    rawDS(currentVal) = [];
  
    % If there are no other entries, disable the Remove and Plot button
    % and change the list sting to <empty>
    if isequal(numResults,length(currentVal)),
	    resultsStr = {'<empty>'};
	    currentVal = 1;	
    end

    % Ensure that list box Value is valid, then reset Value and String
    currentVal = min(currentVal,size(resultsStr,1));
    set(handles.selectstim,'Value',currentVal,'String',resultsStr)

    % Store the new ResultsData
    guidata(h,handles)    


% --------------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
   delete(handles.figure1);

% --------------------------------------------------------------------
function varargout = updatepath_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Stub for Callback of the uicontrol handles.update_stimpath.
    prompt={['Enter the path of input files:']};
    currentpwd = pwd;
    % def={[currentpwd '/Data']};
    def = {currentpwd};
    dlgTitle='Please input stimfile path ';
    lineNo=1;
    datadir =inputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(datadir)
        errordlg('You dont enter the path. Bye.')
        return
    end
    datafiledir = datadir{1};
    if not(exist(datafiledir,'dir'))
         errordlg('Directory not found, Please try again.','Input Error', 'modal')
         return
    end
    handles.filepath = datafiledir;
    load_stimlistbox(datafiledir, handles);
    cd (handles.orgpath)
% --------------------------------------------------------------------
function varargout = updateresppath_Callback(h, eventdata, handles, varargin)
    prompt={['Enter the path of input files:']};
    currentpwd = pwd;
    % def={[currentpwd '/Data']};
    def = {currentpwd};
    dlgTitle='Please input stimfile path ';
    lineNo=1;
    datadir =inputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(datadir)
        errordlg('You dont enter the path. Bye.')
        return
    end

    datafiledir = datadir{1};
    if not(exist(datafiledir,'dir'))
         errordlg('Directory not found, Please try again.','Input Error', 'modal')
         return
    end
    handles.psthfilepath = datafiledir;
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
    
    set(handles.stimlist,'String',handles.stimfile_names,...
            'Value', 1)
    set(handles.resplist,'String',handles.respfile_names,...
            'Value', 1) 
    set(handles.stimpath, 'String',dir_path)
    set(handles.resppath, 'String',dir_path)

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
    set(handles.stimlist,'String',handles.stimfile_names,...
            'Value', 1)
    set(handles.stimpath, 'String',dir_path)
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
    
    guidata(handles.figure1, handles)
    set(handles.resplist,'String',handles.respfile_names,...
            'Value', 1)
    set(handles.resppath, 'String',dir_path)
% ====================================================================
%   END of loadrawfile.m
% ====================================================================









