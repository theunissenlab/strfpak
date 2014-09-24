function varargout = songwave_specgram(varargin)
% SONGWAVE_SPECTRUM_PRED Application M-file for songwave_specgram.fig
%    FIG = SONGWAVE_SPECTRUM_PRED launch songwave_specgram GUI.
%    SONGWAVE_SPECTRUM_PRED('callback_name', ...) invoke the named callback.

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
    global rawData rawpredDS sDim
    if rawData == 0
        tt=errordlg('This option is only for raw data.',...
            'Data Selection Error','modal');
        uiwait(tt);
        delete(handles.figure1);
        return;
    else
        if isempty(rawpredDS)
            errordlg('You need select raw validation data files first.',...
                'Selection Error', 'modal')
            return;
        end
    end
    
    if ~isempty(sDim)
        if sDim == '2-D'
            tt=errordlg('This option is only for sound files.',...
                'Data Selection Error','modal');
            uiwait(tt);
            delete(handles.figure1);
            return;
        end
    end
    global fwidthHz stimsamprate fstep filteroption
    global predampsamprate respsamprate psth_smooth NBAND
    if ~isempty(fwidthHz)
        set(handles.fwidthHz, 'String', fwidthHz);
        set(handles.fwidthms, 'String', 1000/(2*pi*fwidthHz));
        set(handles.predampsamprate, 'String', predampsamprate);
        set(handles.respsamprate, 'String', respsamprate);
        set(handles.samprate, 'String', stimsamprate);
        set(handles.smoothconst, 'String', psth_smooth);
        set(handles.fstep, 'String', NBAND);
    else
        fwidthHz = 125;
        ampsamprate = 1000;
        respsamprate = 1000;
        
        set(handles.fwidthHz, 'String', fwidthHz);
        set(handles.fwidthms, 'String', 1000/(2*pi*fwidthHz));
        set(handles.ampsamprate, 'String', ampsamprate);
        set(handles.respsamprate, 'String', respsamprate);
    end
    
    if ~isempty(filteroption)
        set(handles.filteroption, 'value', filteroption+1);
    end
    
    global rawpredDS predpredinitialFreq predpredendFreq
    if ~isempty(rawpredDS)
        [path, name, ext] = fileparts(rawpredDS{handles.index}.stimfiles); 
        set(handles.filename, 'String',[name ext]);
        [path, name, ext] = fileparts(rawpredDS{handles.index}.respfiles); 
        set(handles.rfilename, 'String',[name ext]);
        set(handles.lowfreq, 'String', predpredinitialFreq);
        set(handles.upfreq, 'String', predendFreq);
    else
        predinitialFreq = 250;
        predendFreq = 8000;
        set(handles.lowfreq, 'String', predpredinitialFreq);
        set(handles.upfreq, 'String', predendFreq);
    end
    
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
function varargout = fwidthHz_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    global fwidthHz predampsamprate
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
    fwidthHz = NewVal;
    predampsamprate = 10* fwidthHz;

    set(handles.fwidthms, 'String', 1000/(2*pi*fwidthHz));
    set(handles.predampsamprate, 'String', predampsamprate);
    
    

% --------------------------------------------------------------------
function varargout = fwidthms_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
     global fwidthHz;
     if ~isempty(fwidthHz)
         set(handles.fwidthms, 'String', 1000/(2*pi*fwidthHz));
         return;
     else
         newStr = get(h, 'String');
         NewVal = str2double(newStr);

        % Check that the entered value falls within the allowable range
        if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
            % Set default value for resp sample rate
            errordlg('Please enter valid frequency bandwith (positive number only).',...
                'Variable Error', 'modal')
            return; 
        end
     
        fwidthHz = 1000/(2*pi*NewVal);
        set(handles.fwidthHz, 'String', fwidthHz);
     end

% --------------------------------------------------------------------
function varargout = predampsamprate_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
     global predampsamprate fwidthHz
     if ~isempty(fwidthHz)
         fwhz = fwidthHz;
     else
         fwhz = 250;
     end
     newStr = get(h, 'String');
     NewVal = str2double(newStr);

     % Check that the entered value falls within the allowable range
     if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
         % Set default value for resp sample rate
         errordlg('Please enter valid frequency bandwith (positive number only).',...
                'Variable Error', 'modal')
         return; 
     elseif NewVal < 2*pi*sqrt(2)*fwhz
         ttt=warndlg('You are undersampling your data.', 'amp_samp_rate warning','modal');
         uiwait(ttt);
     end
     
     predampsamprate = NewVal;
     
% --------------------------------------------------------------------
function varargout = respsamprate_Callback(h, eventdata, handles, varargin)

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
   % 1. Check if we have valid parameter
   %
   global respsamprate
   if isempty(respsamprate)
       errordlg('Please enter your response data sampling rate.',...
           'Data Input Error', 'modal')
       return
   end
   global fwidthHz filteroption predampsamprate
   if isempty(fwidthHz)  | isempty(filteroption) | isempty(predampsamprate)
       errordlg('Please enter valid parameters first.',...
            'Data Input Error', 'modal')
        return
    end

    global rawpredDS
    if isempty(rawpredDS)
        errordlg('Please select wave files first.',...
            'Data input error', 'modal')
        return
    end
    
    % 2. Ask where you want to put your intermediate results
    %
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
                   'Warning Message', 'Overwrite the directory', 'Create new directory',...
                   'Create new directory');
               switch anw
                   case 'Overwrite the directory'
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
                            tt =warndlg('Directory exists. Please create new one.', 'Warning', 'modal');
                            uiwait(tt);
                           end
                       end
                   otherwise
                       errordlg('Computation is cancelled.', 'modal')
                       return;
               end
            end
        end

        outputPath = outdatafiledir;
    else   % not empty outputPath
        outdatafiledir = outputPath;
        tt = dir(outdatafiledir);
        if length(tt) >2    % more file there

            anw = questdlg({'The directory is not empty. What do you want to do?'},...
                'Warning Message', 'Overwrite the directory', 'Create new directory',...
                'Create new directory');
            switch anw
                case 'Overwrite the directory'
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
                            tt =warndlg('Directory exists. Please create new one.', 'Warning', 'modal');
                            uiwait(tt);
                        end
                    end
                otherwise
                    errordlg('Computation is cancelled.', 'modal')
                    return;
            end
        end
        outputPath = outdatafiledir;
    end
    
    set(handles.outdirshow, 'String', outputPath);
    % 3. Call ComplexSpetrum to compute specgram for input
    %     and assign them to global variable predDS
    %
    numfiles = length(rawpredDS);
    global predDS NBAND num_trials
    global stimsamprate predpredinitialFreq predendFreq sDim
    sDim='1-D';
        
    %outSpectrum = cell(numfiles, 1);
    %stimsamprate = cell(numfiles,1);
    currentPath = pwd;
    DBNOISE = 80;
    % do calculation
    tempWait = waitbar(0, 'Calculating spectrogram, please wait...');
    for ii=1:numfiles
      
        waitbar(ii/numfiles, tempWait);
        
        % 3.1. Take care of Stimulus file
        [input, fs] = wavread(rawpredDS{ii}.stimfiles);
        [path,name,ext,ver] = fileparts(rawpredDS{ii}.stimfiles);
        % Spectrogram calculation
        %
       [yy, xx] = ComplexSpectrum(input, floor(fs/predampsamprate),...
               floor(1/(2*pi*fwidthHz)*6*fs),fs);
        
        stimsamprate = fs;
        % Working on more options -JXZ- Aug. 2003
        if filteroption == 1 % linear-linear 
            tmp = abs(yy);
        else  % just temporary
            %tmp = log(abs(yy) +1);
            tmp = max(0,20*log10(abs(yy)./max(max(abs(yy))))+DBNOISE);
        end
       if isempty(predpredinitialFreq) & isempty(predendFreq)
           predpredinitialFreq = 1;
           predendFreq = fs/2;
       else
           if predendFreq > fs/2
               ttt=warndlg('Max frequency limit is stimsamprate/2.', 'high frequency warning','modal');
               uiwait(ttt);
               predendFreq = fs/2;
           end
       end
        freq_range = find(xx>=predinitialFreq & xx<=predendFreq);
        
%         if rem(size(tmp,1), 2)    % Odd
%             cutoff = (size(tmp,1) +1)/2;
%         else
%             cutoff = size(tmp,1)/2+1;
%         end

        outSpectrum = tmp(freq_range,:);
        fo = xx(freq_range);
 
        save(fullfile(outputPath,[name,'_Stim_',num2str(ii),'.mat']), 'outSpectrum');
        % Assign values to global variable predDS
        predDS{ii}.stimfiles = fullfile(outputPath,[name,'_Stim_',num2str(ii),'.mat']);
        predDS{ii}.nlen = round(length(input)*1000/fs) +1;
        
        %save(fullfile(outputPath,['specgram_Stimlabel_',num2str(ii),'.mat']), 'fo'); 
 
        % 3.2. Take care of Response file
        %    If you have multiple trial data, calculate psth first
        %    Then resample it using new amp_samp_rate
        
        %rawResp = load(rawpredDS{ii}.respfiles);
        % Modified by Junli, 2003 to read new spike arrivial time file
        % 
        [rawResp, trials] = read_spikeTime_2cell(rawpredDS{ii}.respfiles, round(length(input)*1000/fs) +1);
        [path,name,ext,ver] = fileparts(rawpredDS{ii}.respfiles);
        %rawResp = read_spikeTime(rawpredDS{ii}.respfiles, size(outSpectrum, 2));
                
        
        % save to the file for each data pair
        save(fullfile(outputPath,[name,'_Spike_time_',num2str(ii),'.mat']), 'rawResp');
 
        % Assign values to global variable predDS
        predDS{ii}.respfiles = fullfile(outputPath,[name,'_Spike_time_',num2str(ii),'.mat']);
        predDS{ii}.ntrials = trials;
        
    end
    close(tempWait)
    
    NBAND = size(outSpectrum, 1);
    save(fullfile(outputPath,'specgram_parameters.mat'),...
        'rawpredDS', 'fwidthHz', 'filteroption','stimsamprate','predampsamprate');
    set(handles.outdirshow, 'visible', 'on', 'String', outputPath);
    set(handles.outdirBrowser, 'visible', 'on');
    msgbox('Done computing and saving.', 'modal');
    
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
    clear global fstep fwidthHz filteroption outputPath 
    clear global respsamprate stimsamprate psth_smooth
    clear global predDS NBAND num_trials predendFreq predinitialFreq predampsamprate sDim
   
    set(handles.fwidthHz, 'String', ' ');
    set(handles.fwidthms, 'String', ' ');
    set(handles.predampsamprate, 'String', ' ');
    set(handles.respsamprate, 'String', ' ');
    set(handles.samprate, 'String', ' ');
    set(handles.smoothconst, 'String', ' ');
    set(handles.upfreq, 'String', ' ');
    set(handles.lowfreq, 'String',' ');
    set(handles.fstep, 'String', ' ');
    set(handles.filename, 'String', ' ');
    set(handles.rfilename, 'String', ' ');
    
    %v = get(handles.dimension, 'String');
    set(handles.filteroption, 'value',1);
    set(handles.stim, 'visible', 'off');
    child = get(handles.stim, 'Children');
    set(child, 'Visible', 'off')
    set(handles.psth, 'visible', 'off');
    child = get(handles.psth, 'Children');
    set(child, 'Visible', 'off')
    
    set(handles.prev, 'visible', 'off');
    set(handles.next, 'visible', 'off');
    set(handles.smoothpsth,'visible','off');
    set(handles.mstag, 'visible', 'off');
    set(handles.smoothconst,'visible','off');
    set(handles.outdirshow, 'visible', 'off', 'string', '');
    set(handles.outdirBrowser, 'visible', 'off')

    % Save all handles to figures
    guidata(h,handles);
% --------------------------------------------------------------------
function varargout = samprate_Callback(h, eventdata, handles, varargin)
    
    
% --------------------------------------------------------------------  
function displayspecgram(handles)
    
    %  Prepare for axes
    %
    global stimsamprate predDS psth_smooth predampsamprate
    if isempty(psth_smooth)
        psth_smooth = 30;
    end
    global outputPath fstep respsamprate
    
    forward = Check_And_Load(predDS{handles.index}.stimfiles);
    rawResp = Check_And_Load(predDS{handles.index}.respfiles);
    
    spiketrain = zeros(predDS{handles.index}.ntrials,predDS{handles.index}.nlen);
    for trial_ind =1:predDS{handles.index}.ntrials
        
        spiketrain(trial_ind, rawResp{trial_ind}) = ones(1,length(rawResp{trial_ind}));
    end
    newpsth = resample(spiketrain', predampsamprate, respsamprate);
    
    newpsth = newpsth'; % make sure new response data is trials x T.
    newpsth(find(newpsth < 0)) = 0;
    
    [xm ym] = size(forward);
    
    maxforward = max(max(forward));
    minforward = min(min(forward));
    absforward = max(abs(minforward),abs(maxforward));
    
    % Now display input stimulus file
    global predampsamprate
    nb = xm;
    nt = min(ym, size(newpsth,2));
    
    xlabelrange = ceil((1:nt)*1000/predampsamprate);
   
    global predinitialFreq predendFreq predampsamprate
        
    fstep = (predendFreq - predinitialFreq)/nb;
    faxis = predinitialFreq:fstep:predendFreq;
    axes(handles.stim);
    %imagesc(xlabelrange, faxis/1000, 3*forward(:,1:nt), [-absforward absforward]);
    imagesc(xlabelrange, faxis/1000, forward(:,1:nt));
    axis xy;
    ylabel('Frequency (kHz)')
    title('Stimulus')
    
    set(handles.fstep, 'String', xm);
    set(handles.samprate, 'String', stimsamprate);
    set(handles.smoothconst, 'String', psth_smooth);
    set(handles.stim, 'visible', 'on');
    set(handles.smoothpsth, 'visible','on');
    set(handles.smoothconst, 'visible','on');
    set(handles.mstag, 'visible', 'on');
    set(handles.upfreq, 'String', predendFreq);
    set(handles.lowfreq, 'String', predinitialFreq);
    
    %imagesc(1:size(forward,2), fo, forward,[-absforward absforward])
    %imagesc(1:size(forward,2), fo(1:fix(size(fo)/2)), 10*forward);
    %,[-absforward absforward]);
  

    % Display raw psth
    if size(newpsth) > 1
        spikeavg = mean(newpsth);
        %plot(xlabelrange, mean(newpsth)*predampsamprate);
        %plot(xlabelrange, mean(newpsth));
    else
        spikeavg = newpsth;
    end
    
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
    plot(xlabelrange,svagsm(halfwinsize+cutsize:length(svagsm)-halfwinsize)*predampsamprate);
    axis xy;
    ylabel('PSTH (Spikes/sec)')
    %xlabel(sprintf('Time ( in %2.2f ms)', 1000/predampsamprate))
    xlabel('Time (in ms)')
    p_axis = axis;
    p_axis(2) = ceil(nt*1000/predampsamprate);
    p_axis(4) = max(svagsm)*predampsamprate;
    axis(p_axis);
    
    set(handles.prev, 'visible', 'off');
    set(handles.next, 'visible', 'off');
    global rawpredDS
    if ~isempty(rawpredDS)
        [path, name, ext] = fileparts(rawpredDS{handles.index}.stimfiles); 
        set(handles.filename, 'String',[name ext]);
        [path, name, ext] = fileparts(rawpredDS{handles.index}.respfiles); 
        set(handles.rfilename, 'String',[name ext]);
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
        '   Short-time Fourier Transform (Help window)                               '
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
function upfreq_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor','white');



function upfreq_Callback(hObject, eventdata, handles)
newStr_respsamprate = get(hObject, 'String');
NewVal = str2double(newStr_respsamprate);

% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid high freq (positive number only).',...
        'Variable Error', 'modal')
    return; 
end

% Assign global variable to predampsamprate
global predendFreq;
if isempty(predendFreq)
    predendFreq = NewVal;
else
    predendFreq = NewVal;
    compsave_Callback(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function lowfreq_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor','white');




function lowfreq_Callback(hObject, eventdata, handles)

newStr_respsamprate = get(hObject, 'String');
NewVal = str2double(newStr_respsamprate);

% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid low frequency (positive number only).',...
        'Variable Error', 'modal')
    return; 
end

% Assign global variable to predampsamprate
global predinitialFreq;
if isempty(predinitialFreq)
    predinitialFreq = NewVal;
else
    predinitialFreq = NewVal;
    compsave_Callback(hObject, eventdata, handles);
end

% ====================================================================
%   END of songwave_specgram.m
% ====================================================================


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
% hObject    handle to datasetSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
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


