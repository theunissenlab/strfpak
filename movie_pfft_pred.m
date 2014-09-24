function varargout = movie_pfft_pred(varargin)
% MOVIE_PFFT_PRED Application M-file for songwave_specgram.fig
%    FIG = MOVIE_PFFT_PRED launch songwave_specgram GUI.
%    MOVIE_PFFT_PRED('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 03-Apr-2006 14:01:47

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
	guidata(fig, handles);
  
    % Check whether the input is right for preprocessing
    global rawpredData rawpredDS sDim
    if rawpredData == 0
        tt=errordlg('This option is only for raw data.',...
            'Data Selection Error','modal');
        uiwait(tt);
        delete(handles.figure1);
        return;
    else
        if isempty(rawpredDS)
            errordlg('You need select data files first.',...
                'Selection Error', 'modal')
            return;
        end
    end
    
    if ~isempty(sDim)
        if sDim == '1-D'
            tt=errordlg('This option is only for 2-D (movie) files.',...
                'Data Selection Error','modal');
            uiwait(tt);
            delete(handles.figure1);
            return;
        end
    end
 
    if ~isempty(rawpredDS)
        [path, name, ext] = fileparts(rawpredDS{handles.index}.stimfiles); 
        set(handles.filename, 'String',[name ext]);
        [path, name, ext] = fileparts(rawpredDS{handles.index}.respfiles); 
        set(handles.rfilename, 'String',[name ext]);
        set(handles.inputdirshow, 'String', path);
    end
    
    global startidx_v stopidx_v bwindow_v btakepower_v bzdc_v btempnl_v
    if isempty(startidx_v) | isempty(stopidx_v) | isempty(bwindow_v)...
            | isempty(btakepower_v) | isempty(bzdc_v) | isempty(btempnl_v)

        startidx_v  = 0;
        stopidx_v = 0;
        bwindow_v = 0;
        btakepower_v = 1;
        bzdc_v = 0;
        btempnl_v = 0;
    end
    set(handles.startidx, 'String', startidx_v);
    set(handles.stopidx, 'String', stopidx_v);
    set(handles.bwindow, 'String', bwindow_v);
    set(handles.btakepower, 'String', btakepower_v);
    set(handles.btempnl, 'String', btempnl_v);
    set(handles.bzdc, 'String', bzdc_v);
    
    % Set values for dataset
    if length(rawpredDS) ==1
        set(handles.datasetSlider, 'min', 0, 'max', length(rawpredDS), 'value',0,'sliderstep', [1 1]);
    elseif length(rawpredDS) > 1
        set(handles.datasetSlider, 'Min', 1, 'Max', length(rawpredDS),'value',1,...
            'sliderstep', [1/(length(rawpredDS)-1) 1/(length(rawpredDS)-1)]);
    end
    set(handles.DataSet, 'string', num2str(handles.index));
    
    global outputPath
    if ~isempty(outputPath)
        set(handles.outdirshow, 'visible', 'on', 'String', outputPath);
        set(handles.outdirBrowser, 'visible', 'on');
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
    global rawpredDS
    if ~isempty(rawpredDS)
        [path, name, ext] = fileparts(rawpredDS{handles.index}.stimfiles); 
        set(handles.filename, 'String',[name ext]);
        [path, name, ext] = fileparts(rawpredDS{handles.index}.respfiles); 
        set(handles.rfilename, 'String',[name ext]);
    end
    guidata(h, handles);
% --------------------------------------------------------------------
function varargout = loadresp_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
   %msgbox('Coming soon...', 'modal');
   loadrawfile;

% --------------------------------------------------------------------
function varargout = compsave_Callback(h, eventdata, handles, varargin)
set(handles.stim, 'visible', 'off')
child = get(handles.stim, 'Children');
set(child, 'Visible', 'off')
set(handles.psth, 'visible', 'off');
child = get(handles.psth, 'Children');
set(child, 'Visible', 'off')

set(handles.smoothpsth,'visible','off');
%set(handles.mstag, 'visible', 'off');
set(handles.smoothconst,'visible','off');

global outputPath
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
            errordlg('Even upper directory not found. existing...');
            return
        end
        cd (p)
        mkdir(n)
        % Junli: 11/11/2005
    else   % existing directory
        tt = dir(outdatafiledir);
        if length(tt) >2    % more file there

            anw = questdlg({'The directory is not empty. What do you want to do?'},...
                'Warning Message', 'Empty the directory', 'Create new directory',...
                'Create new directory');
            switch anw
                case 'Empty the directory'
                    current_dir = pwd;
                    cd(outdatafiledir)
                    delete *_*.mat;
                    cd(current_dir);

                case 'Create new directory'

                    newDir = 0;
                    while newDir == 0
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
                                errordlg('Even upper directory not found. existing...');
                                return
                            end
                            cd (p)
                            mkdir(n)
                            newDir == 1;
                        else
                            tt =warndlg('Directory is exsiting. Please create new one.', 'Warning', 'modal');
                            uiwait(tt);
                        end
                    end
            end
        end
    end

    outputPath = outdatafiledir;
else   % not empty outputPath
    outdatafiledir = outputPath;
    tt = dir(outdatafiledir);
    if length(tt) >2    % more file there

        anw = questdlg({'The directory is not empty. What do you want to do?'},...
            'Warning Message', 'Empty the directory', 'Create new directory',...
            'Create new directory');
        switch anw
            case 'Empty the directory'
                current_dir = pwd;
                cd(outdatafiledir)
                delete *_*.mat;
                cd(current_dir);

            case 'Create new directory'
                newDir = 0;
                while newDir == 0
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
                            errordlg('Even upper directory not found. existing...');
                            return
                        end
                        cd (p)
                        mkdir(n)
                        newDir = 1;
                    else
                        tt =warndlg('Directory is exsiting. Please create new one.', 'Warning', 'modal');
                        uiwait(tt);
                    end
                end
        end
    end
    %uiwait(anw);
    outputPath = outdatafiledir;
end
set(handles.outdirshow, 'visible', 'on');
set(handles.outdirshow, 'String', outputPath);
set(handles.outdirBrowser, 'visible', 'on');

% Now go through all the data
global rawpredDS predDS
numfiles = length(rawpredDS);
global NBAND num_trials
% axes(handles.waitaxes)
% ttH=plot(ones(1, numfiles), 'k');
% set(get(ttH,'Parent'),'YTickLabel',[]);
% set(get(ttH,'Parent'),'YTick',[]);
% set(get(ttH,'Parent'),'XTickLabel',[]);
% set(get(ttH,'Parent'),'XTick',[]);
% set(handles.inprogress, 'visible', 'on');
tempWait = waitbar(0,...
     sprintf('Preprocessing 2-D movie data using Fourier Power Model'));
 for ii=1:numfiles
     %     patch([0 ii ii 0], [0 0 1 1], 'r');
     waitbar(ii/numfiles, tempWait);
     [path,name,ext,ver] = fileparts(rawpredDS{ii}.stimfiles);
     switch ext
         case {'.wav'}
             errordlg('This option is only for movie. ','Data Type Error','modal');
             close(tempWait);
             return;
         case {'.mat'}  % 2-D raw data (movie)
             stimstim = load(rawpredDS{ii}.stimfiles);
             flds = fieldnames(stimstim);
             if (length(flds) == 1)
                 stimstim = getfield(stimstim, char(flds{1}));
             end

         case {'.dat', '.txt'}
             stimstim = dlmread(rawpredDS{ii}.stimfiles);

     end
     if length(size(stimstim))>= 3
         sDim = '2-D';
     else
         sDim = '1-D';
         errordlg('This option is only for movie. ','Data Type Error','modal');
         close(tempWait);
         return;
     end

     % Reshaping the 2-D movie into one big vector
     [xsize, ysize, framesize] = size(stimstim);
     global NBAND
     global startidx_v stopidx_v bwindow_v btakepower_v bzdc_v btempnl_v

     % Take fourier power of movie
     stim=movpower(stimstim, startidx_v, stopidx_v, bwindow_v,...
         btakepower_v, bzdc_v, btempnl_v);
     [NBAND, predDS{ii}.nlen] = size(stim);

     save(fullfile(outputPath,[name,'_Stim_',num2str(ii),'.mat']), 'stim');
     % Assign values to global variable predDS
     predDS{ii}.stimfiles = fullfile(outputPath,[name,'_Stim_',num2str(ii),'.mat']);
     %predDS{ii}.nlen = framesize;

     % 3.2. Take care of Response file
     %    If you have multiple trial data, calculate psth first
     %    Then resample it using new amp_samp_rate

     %rawResp = load(rawpredDS{ii}.respfiles);
     % Modified by Junli, 2003 to read new spike arrivial time file
     %
     [rpath,rname,rext,rver] = fileparts(rawpredDS{ii}.respfiles);
     switch rext
         case {'.dat', '.txt', ''}
             newpsth = load(rawpredDS{ii}.respfiles);

         case {'.mat'}
             respMat = load(rawpredDS{ii}.respfiles);

             % Validate the MAT-file
             flds = fieldnames(respMat);

             % Check if response data is in spike arrival time
             % or already preprocessed.
             if (length(flds) == 1)
                 rawResp = getfield(respMat, char(flds{1}));
                 if iscell(rawResp)
                     spiketrain = zeros(rawpredDS{ii}.ntrials,rawpredDS{ii}.nlen);
                     for trial_ind =1:rawpredDS{ii}.ntrials

                         spiketrain(trial_ind, rawResp{trial_ind}) = ones(1,length(rawResp{trial_ind}));
                     end
                     newpsth = resample(spiketrain', predampsamprate, respsamprate);

                     newpsth = newpsth'; % make sure new response data is trials x T.
                     newpsth(find(newpsth < 0)) = 0;
                 else
                     newpsth = rawResp;
                 end
                 [xs, ys] = size(newpsth);
                 if xs > ys   % make sure response data is trials X TimeFrame
                     newpsth = newpsth';
                 end
                 rgoodidx = find(~isnan(newpsth));
                 newpsth = newpsth(rgoodidx);

             end
         otherwise
             errordlg(['Only ASCII and MAT binary fileformats work',...
                 ' in this version.'], 'File Type Error', 'modal')
             return;
     end

     % save to the file for each data pair
     save(fullfile(outputPath,[name,'_Spike_time_',num2str(ii),'.mat']), 'newpsth');

     % Assign values to global variable predDS
     predDS{ii}.respfiles = fullfile(outputPath,[name,'_Spike_time_',num2str(ii),'.mat']);
     predDS{ii}.ntrials = size(newpsth, 1);

 end
 handles.compsaveflag = 1;
guidata(handles.figure1, handles);
close(tempWait);

% --------------------------------------------------------------------
function varargout = display_Callback(h, eventdata, handles, varargin)
    % 1. Load display result
    %
    
    handles.index = 1;
    guidata(h, handles);
    displayspecgram(handles);

% --------------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
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

    global rawpredDS;
    index = handles.index;
    % update index
	i = index - 1;
	% If the index is less then one then set it equal to the index of the
	% last element in the Addresses array
	if i < 1
	    i = length(rawpredDS);
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
	if i > length(rawpredDS)
	    i = 1;
	end

	handles.index = i;
    guidata(h, handles);
    
    displayspecgram(handles)

% --------------------------------------------------------------------
function varargout = reset_Callback(h, eventdata, handles, varargin)

    % 5/24/2005: Also delete the saved files
    current_dir = pwd;
    global outputPath
    cd(outputPath)
    delete *_Spike_time*;
    delete *_Stim_*;
    cd(current_dir);

    %clear predDS predinitialFreq predendFreq StimSampRate sDim
    clear global startidx_v stopidx_v bwindow_v btakepower_v bzdc_v btempnl_v
    clear global predDS NBAND num_trials  sDim
   
    set(handles.startidx, 'String', ' ');
    set(handles.stopidx, 'String', ' ');
    set(handles.bwindow, 'String', ' ');
    set(handles.btakepower, 'String', ' ');
    set(handles.smoothconst, 'String', ' ');
    set(handles.bzdc, 'String', ' ');
    set(handles.btempnl, 'String',' ');
    
    %v = get(handles.dimension, 'String');
   
    set(handles.stim, 'visible', 'off');
    child = get(handles.stim, 'Children');
    set(child, 'Visible', 'off')
    set(handles.psth, 'visible', 'off');
    child = get(handles.psth, 'Children');
    set(child, 'Visible', 'off')
    
    set(handles.smoothpsth,'visible','off');
    %set(handles.mstag, 'visible', 'off');
    set(handles.smoothconst,'visible','off');
    set(handles.outdirshow, 'visible', 'off', 'string', '');
    set(handles.outdirBrowser, 'visible', 'off')
    

    % Save all handles to figures
    guidata(h,handles);
% --------------------------------------------------------------------
function varargout = samprate_Callback(h, eventdata, handles, varargin)
    
    
% --------------------------------------------------------------------  
function displayspecgram(handles)
if ~isfield(handles, 'compsaveflag')
    errordlg('Please do compution first', 'Order error', 'modal');
    return;
end
    
    %  Prepare for axes
    %
    global predDS psth_smooth
    global outputPath
    
    if isempty(psth_smooth)
        psth_smooth = 21;
    end
    
    forward = Check_And_Load(predDS{handles.index}.stimfiles);
    newpsth = Check_And_Load(predDS{handles.index}.respfiles);
    newpsth(find(newpsth < 0)) = 0;
    
    [xm nt] = size(forward);
    
    maxforward = max(max(forward));
    minforward = min(min(forward));
    absforward = max(abs(minforward),abs(maxforward));
    
    % Now display input stimulus file
    
    axes(handles.stim);
    %imagesc(xlabelrange, faxis/1000, 3*forward(:,1:nt), [-absforward absforward]);
    imagesc(forward,[-absforward absforward] );
    axis xy;
    title('Stimulus')
    
    set(handles.smoothconst, 'String', psth_smooth);
    set(handles.stim, 'visible', 'on');
    set(handles.smoothpsth, 'visible','on');
    set(handles.smoothconst, 'visible','on');
    %set(handles.mstag, 'visible', 'on');
   

    % Display raw psth
    if size(newpsth) > 1
        spikeavg = mean(newpsth);
        %plot(xlabelrange, mean(newpsth)*bwindow);
        %plot(xlabelrange, mean(newpsth));
    else
        spikeavg = newpsth;
    end
    
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
    halfwinsize = floor(psth_smooth/2);
    wind1 = hanning(psth_smooth)/sum(hanning(psth_smooth));
    svagsm=conv(spikeavg(1:nt),wind1);
    axes(handles.psth);
    plot(svagsm(halfwinsize+cutsize:length(svagsm)-halfwinsize));
    axis xy;
    ylabel('PSTH (Spikes)')
    %xlabel(sprintf('Time ( in %2.2f ms)', 1000/bwindow))
    xlabel('Frames ')
    p_axis = axis;
    
    
   
    global rawpredDS
    if ~isempty(rawpredDS)
        [path, name, ext] = fileparts(rawpredDS{handles.index}.stimfiles); 
        set(handles.filename, 'String',[name ext]);
        [path, name, ext] = fileparts(rawpredDS{handles.index}.respfiles); 
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
        '   Fourier Power Model (Help window)                                        '
        '                                                                            '
        '  Load files: To load data sets for short-time fourier transform. Here one  '
        '          data set need include stimulus and its corresponding response     '
        '          data, although the response data may not be changed by STFT.      '
        '  stimulus file:  shows the filename of one stimulus selected.              '
        '  response file:  shows the filename of the response selected.              '
        '                                                                            '
        '  Parameters:                                                               '
        '      filter_width (Hz): the width of the filter in Hz. It defines window   ' 
        '              length of filter (here we use Gaussian filter).               '
        '      filter_width (ms): the width of the filter in ms                      '
        '      amp_samprate (Hz): the sampling rate that we want for the amplitude   '
        '                  envelope. By default, it is 10 times of filter_width (Hz).'
        '      resp_samprate (Hz): the sampling rate of spike data.                  '
        '      upper frequency (Hz): the upper frequency you want to study.          '
        '                           The max frequency is limited to stim_samprate/2. '
        '      lower frequency (Hz): the lower frequency you want to study (>= 0Hz). '
        '      stim_samprate (Hz): the sampling rate of stimulus in Hz.              '
        '      nbands: the numbers of frequency bands covered for the calculation.   '
        '      Scale option: the choice of linear scale or logarithmic scale         '
        '                                                                            '
        '  Compute and Save: Compute the spectrogram of the signal and save          '
        '       the results into the directory where you will be asked to input.     '
        '       The computing status bar also shows up so that you can know progress.'
        '  Display: Graphically display the spectrogram of the stimulus and          '
        '      the smoothed psth with  Smooth_PSTH window size. The smoothing width  '
        '      for psth is shown here. You can modify it by typing different number. '
        '       If more than one data sets are chosen, \fbox{Next} and               '
        '       Prev buttons show up so that you can click to see next data sets.    '
        '  Reset: Reset all the parameters and the data sets chosen.                 '
        '  Close: Close this window and save all the parameters and all the result   '
        '                                                                            '
        ' FOR DETAILED INFORMATION, PLEASE REFER TO STRFPAK DOCUMENTATION.           '
        '                                                                            '
        ' Updated by Junli, Sept. 2005.                                              '
    ];
    myFig = handles.figure1;
    helpwin(hlpStr, ttlStr);


% --- Executes during object creation, after setting all properties.
function bzdc_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor','white');



function bzdc_Callback(hObject, eventdata, handles)
newStr_respsamprate = get(hObject, 'String');
NewVal = str2double(newStr_respsamprate);

% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid high freq (positive number only).',...
        'Variable Error', 'modal')
    return; 
end

% Assign global variable to bwindow
global predendFreq;
if isempty(predendFreq)
    predendFreq = NewVal;
else
    predendFreq = NewVal;
    compsave_Callback(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function btempnl_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');


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
global predinitialFreq;
if isempty(predinitialFreq)
    predinitialFreq = NewVal;
else
    predinitialFreq = NewVal;
    compsave_Callback(hObject, eventdata, handles);
end

% --- Executes on slider movement.
function datasetSlider_Callback(hObject, eventdata, handles)
newSlider = round(get(hObject, 'Value'));
if newSlider <= 0
    set(handles.DataSet, 'String', handles.index);
else
    set(handles.DataSet, 'String', newSlider);
    handles.index = newSlider;
end
guidata(handles.figure1, handles);
displayspecgram(handles)

% --- Executes during object creation, after setting all properties.
function datasetSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function DataSet_Callback(hObject, eventdata, handles)
global rawpredDS
maxlength = length(rawpredDS);
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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in outdirBrowser.
function outdirBrowser_Callback(hObject, eventdata, handles)
global outputPath
outputdir = uigetdir(pwd, 'Pick the response data directory');
if not(exist(outputdir,'dir'))
     errordlg('There is no selected directory', 'Output Dir Error', 'modal');
     return;
end
set(handles.outdirshow, 'String', outputdir);
outputPath = outputdir;
compsave_Callback(hObject, eventdata, handles);;


function outdirshow_Callback(hObject, eventdata, handles)
 global outputPath
newVal = get(hObject, 'String');
if not(exist(newVal,'dir'))
    tt=warndlg('Directory not found. Creating new directory.','dir error', 'modal');
    uiwait(tt);
    [p, n, e] = fileparts(newVal);
    if not(exist(p, 'dir'))
        errordlg('Even upper directory not found. existing...','Dir error','modal');
        return
    end
    cd (p)
    mkdir(n)
    set(handles.outdirshow, 'String', newVal);
    outputPath = newVal;
    compsave_Callback(hObject, eventdata, handles);
else
    outputPath = newVal;
    compsave_Callback(hObject, eventdata, handles);
end


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


