% --------------------------------------------------------------------
function varargout = displayraw_preprocessed_GUI(varargin)
% --------------------------------------------------------------------
%
% DISPLAYRAW_PREPROCESSED_GUI Application M-file for displayraw_preprocessed_GUI.fig
%    FIG = DISPLAYRAW_PREPROCESSED_GUI launch displayraw_preprocessed_GUI GUI.
%    DISPLAYRAW_PREPROCESSED_GUI('callback_name', ...) invoke the named callback.
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

% Created by JXZ, 2002.
%
if nargin == 0  % LAUNCH GUI

    fig = openfig(mfilename,'reuse');

    % for resize property
    set(fig, 'resize', 'on');
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
           hUIControls],'units','normalized',...
          'fontname', 'Times New Roman','fontunits','normalized');

    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    handles.frameIndex = 0;
    handles.trialIndex = 0;
    handles.displaywpsth = 0;
    guidata(fig, handles);
       
    % load the first data file for default display window
    global DS outputPath
    
%     if isempty(outputPath)
%         currentPath = pwd;
%         prompt={['Please Enter where you want to put your intermediate results:']};
%         def = {fullfile(currentPath, 'Output')};
%         dlgTitle='Path for intermediate results';
%         lineNo=1;
% 
%         % picture feature
%         AddOpts.Resize='on';
%         AddOpts.WindowStyle='normal';
%         AddOpts.Interpreter='tex';
%         datadir =inputdlg2(prompt,dlgTitle,lineNo,def,AddOpts);
% 
%         % Check if user input valid directory
%         if isempty(datadir)
%             errordlg('Please enter valid directory.','Input Error','modal')
%             return
%         end
%     
%         outdatafiledir = datadir{1};
%         if not(exist(outdatafiledir,'dir'))
%             disp('Directory not found. Creating new directory.');
%             [p, n, e] = fileparts(outdatafiledir);
%             if not(exist(p, 'dir'))
%                 errordlg('Even upper directory not found. existing...');
%                 return
%             end
%             cd (p)
%             mkdir(n)
%         else   % existing directory
%             tt = dir(outdatafiledir);
%             if length(tt) >2    % more file there
% 
%                 anw = questdlg({'The directory is not empty. What do you want to do?'},...
%                     'Warning Message', 'Empty the directory', 'Create new directory',...
%                     'Create new directory');
%                 switch anw
%                     case 'Empty the directory'
%                         current_dir = pwd;
%                         cd(outdatafiledir)
%                         delete *_*.mat;
%                         cd(current_dir);
% 
%                     case 'Create new directory'
% 
%                         newDir = 0;
%                         while newDir == 0
%                             currentPath = pwd;
%                             prompt={['Please Enter where you want to put your intermediate results:']};
%                             def = {fullfile(currentPath, 'Output')};
%                             dlgTitle='Path for intermediate results';
%                             lineNo=1;
% 
%                             % picture feature
%                             AddOpts.Resize='on';
%                             AddOpts.WindowStyle='normal';
%                             AddOpts.Interpreter='tex';
%                             datadir =inputdlg2(prompt,dlgTitle,lineNo,def,AddOpts);
% 
%                             % Check if user input valid directory
%                             if isempty(datadir)
%                                 errordlg('Please enter valid directory.','Input Error','modal')
%                                 return
%                             end
% 
%                             outdatafiledir = datadir{1};
% 
%                             if not(exist(outdatafiledir,'dir'))
%                                 disp('Directory not found. Creating new directory.');
%                                 [p, n, e] = fileparts(outdatafiledir);
%                                 if not(exist(p, 'dir'))
%                                     errordlg('Even upper directory not found. existing...');
%                                     return
%                                 end
%                                 cd (p)
%                                 mkdir(n)
%                                 newDir == 1;
%                             else
%                                 tt =warndlg('Directory is exsiting. Please create new one.', 'Warning', 'modal');
%                                 uiwait(tt);
%                             end
%                         end
%                 end
%             end
%         end
%         outputPath = outdatafiledir;
%     else   % not empty outputPath
%         outdatafiledir = outputPath;
%         tt = dir(outdatafiledir);
%         if length(tt) >2    % more file there
% 
%             anw = questdlg({'The directory is not empty. What do you want to do?'},...
%                 'Warning Message', 'Empty the directory', 'Create new directory',...
%                 'Create new directory');
%             switch anw
%                 case 'Empty the directory'
%                     current_dir = pwd;
%                     cd(outdatafiledir)
%                     delete *_*.mat;
%                     cd(current_dir);
% 
%                 case 'Create new directory'
%                     newDir = 0;
%                     while newDir == 0
%                         currentPath = pwd;
%                         prompt={['Please Enter where you want to put your intermediate results:']};
%                         def = {fullfile(currentPath, 'Output')};
%                         dlgTitle='Path for intermediate results';
%                         lineNo=1;
% 
%                         % picture feature
%                         AddOpts.Resize='on';
%                         AddOpts.WindowStyle='normal';
%                         AddOpts.Interpreter='tex';
%                         datadir =inputdlg2(prompt,dlgTitle,lineNo,def,AddOpts);
% 
%                         % Check if user input valid directory
%                         if isempty(datadir)
%                             errordlg('Please enter valid directory.','Input Error','modal')
%                             return
%                         end
% 
%                         outdatafiledir = datadir{1};
% 
%                         if not(exist(outdatafiledir,'dir'))
%                             disp('Directory not found. Creating new directory.');
%                             [p, n, e] = fileparts(outdatafiledir);
%                             if not(exist(p, 'dir'))
%                                 errordlg('Even upper directory not found. existing...');
%                                 return
%                             end
%                             cd (p)
%                             mkdir(n)
%                             newDir = 1;
%                         else
%                             tt =warndlg('Directory is exsiting. Please create new one.', 'Warning', 'modal');
%                             uiwait(tt);
%                         end
%                     end
%             end
%         end
%         outputPath = outdatafiledir;
%     end

 
    checknload(DS, 1, handles);
    
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
function varargout = next_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    global DS;

    % Get the index pointer and the files 
    index = handles.Index;
    
    % update index
    i = index + 1;
    % If the index is less then one then set it equal to the index of the
    % last element in the Addresses array
    if i > length(DS) 
        i = 1;
    end

    
    % Update the index pointer to reflect the new index
    handles.Index = i;
    % Also reset frameIndex
    handles.frameIndex = 0;
    guidata(h,handles)

    checknload(DS,i,handles);



% --------------------------------------------------------------------
function varargout = close_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    delete(handles.figure1);


% --------------------------------------------------------------------
function varargout = prev_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    global DS;

    % Get the index pointer and the files 
    index = handles.Index;
    
    % update index
    i = index - 1;
    % If the index is less then one then set it equal to the index of the
    % last element in the Addresses array
    if i < 1
        i = length(DS);
    end
    
    
    % Update the index pointer to reflect the new index
    handles.Index = i;
    % Also reset frameIndex
    handles.frameIndex = 0;
    guidata(h,handles)
    
    checknload(DS, i, handles);
    


% --------------------------------------------------------------------
function pass = checknload(DS, n, handles)
% --------------------------------------------------------------------
    
    global ampsamprate respsamprate
    % Check whether we have data files selected.
    if length(DS) == 0 
        errordlg('There are no data files. Please select data files first.',...
                'Data Files Error', 'modal')
        return;
    elseif n > length(DS) | (n < 0)
        errordlg('Index out of DS range.', 'File index wrong', 'modal')
        return;
    end 
    handles.Index = n;
    
    %load stimulus file
    [path,name,ext,ver] = fileparts(DS{n}.stimfiles);
    switch ext
    case {'.dat', '.txt', ''}
        handles.stimfiles = load (DS{n}.stimfiles);
        
    case {'.mat'}
        stimMat = load(DS{n}.stimfiles);   
        
        % Validate the MAT-file
         flds = fieldnames(stimMat);
         if (length(flds) == 1) 
             handles.stimfiles = getfield(stimMat, char(flds{1}));
             
         end
    otherwise
         errordlg(['Only ASCII and MAT binary fileformats work',...
                ' in this version.'], 'File Type Error', 'modal')
         return;
    end
       
    % Get the appropriate data for the index in selected
    set(handles.show_stimfile_text,'String',[name ext]);
    
    %load response file
    [path,name,ext,ver] = fileparts(DS{n}.respfiles);
    switch ext
    case {'.dat', '.txt', ''}
             
         handles.respfiles = load(DS{n}.respfiles);
         
         
    case {'.mat'}
        respMat = load(DS{n}.respfiles);   
        
        % Validate the MAT-file
        flds = fieldnames(respMat);
        
        % Check if response data is in spike arrival time 
        % or already preprocessed.
        if (length(flds) == 1)
            rawResp = getfield(respMat, char(flds{1}));
            if iscell(rawResp)
                spiketrain = zeros(DS{n}.ntrials,DS{n}.nlen);
                for trial_ind =1:DS{n}.ntrials

                    spiketrain(trial_ind, rawResp{trial_ind}) = ones(1,length(rawResp{trial_ind}));
                end
                newpsth = resample(spiketrain', ampsamprate, respsamprate);

                newpsth = newpsth'; % make sure new response data is trials x T.
                newpsth(find(newpsth < 0)) = 0;
            else
                newpsth = rawResp;
            end
            handles.respfiles = newpsth;
            
        end
    otherwise
         errordlg(['Only ASCII and MAT binary fileformats work',...
                ' in this version.'], 'File Type Error', 'modal')
         return;
    end
       
    % Get the appropriate data for the index in selected
    set(handles.show_respfile_text,'String',[name ext]);
    guidata(handles.figure1, handles)
    displayall(handles);
    
% --------------------------------------------------------------------
function displayall(handles)
% --------------------------------------------------------------------
    
    % Now display input stimulus file
    global ampsamprate psth_smooth DS
    if isempty(psth_smooth)
        psth_smooth = 21;
    end
    global sDim 
    set(handles.spatialdomain, 'String', sDim);
    
    set(handles.smoothconst,'String',psth_smooth);
    
    nb = size(handles.stimfiles,1);
    nt = min(size(handles.respfiles,2), size(handles.stimfiles, 2));
  
    xlabelrange = (1:nt)*ceil(1000/ampsamprate);
    
   
    if strcmp(sDim, '1-D')
        % Make 2d frame button invisible
        set(handles.prev5frame, 'Visible', 'off')
        set(handles.next5frame, 'Visible', 'off')
        
        forward = handles.stimfiles(:,1:nt);
        maxforward = max(max(forward));
        minforward = min(min(forward));
        absforward = max(abs(minforward),abs(maxforward));
       
        
        global initialFreq endFreq preprocessOption
        
        if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
            flabel = logspace(log10(initialFreq), log10(endFreq), nb);
            axes(handles.stim_axes);
            pcolor(xlabelrange, flabel/1000, forward); shading interp;
            caxis([-absforward absforward]);
        else
        
            fstep = (endFreq - initialFreq)/nb;
            faxis = initialFreq:fstep:endFreq;
            axes(handles.stim_axes);
            imagesc(xlabelrange, faxis/1000, forward);%, [-absforward absforward]);
            axis xy;
        
        end
        ylabel('Frequency (kHz)')
        title('Stimulus')
        
    else   % Display for 2-D data
        tt = handles.stimfiles;
        binsize= ceil((1/ampsamprate)*1000);
        splitX = sqrt(nb);
        if splitX - floor(splitX) ~= 0
           errordlg(['Data Error: Please check your input data by clicking ',...
               '"Get Files" Button in the main window: The first data file need ',...
               'to be stimuli and the second data file need to be its corresponding',...
               ' response file. If you made a mistake, please type "clear all" ',...
               ' or hit "reset" button first and then choose input data again.'],...
               'Input Data Error', 'modal');
            return;
        end
        
        stimstim = reshape(tt(:,1:nt), splitX, splitX, nt);
        numList = 10;
        whole=[];
        
        %         for i=handles.frameIndex:handles.frameIndex+numList
        %             if i+1 > nt
        %                 handles.frameIndex = 0;
        %                 guidata(handles.figure1, handles)
        %                 break;
        %             end
        %             whole = [whole stimstim(:,:,i+1)];
        %         end
        
        
        frameRange = handles.frameIndex+1:handles.frameIndex+numList;
        
        if length(find(frameRange>nt)) == 0
            whole = reshape(stimstim(:,:,frameRange), splitX, splitX*length(frameRange));
        else
            whole = reshape(stimstim(:,:,frameRange(find(frameRange<=nt))), splitX,...
                splitX*length(frameRange(find(frameRange<=nt))));
        end
        axes(handles.stim_axes);
        ttH = imagesc(whole);
        axis image
        set(get(ttH,'Parent'),'YTickLabel',[]);
        set(get(ttH,'Parent'),'YTick',[]);
        set(get(ttH,'Parent'),'XTickLabel',[]);
        set(get(ttH,'Parent'),'XTick',[]);
        
        % 5/12/03 - JXZ
        colormap(redblue)
        title([num2str(handles.frameIndex+1) '  to '...
                num2str(handles.frameIndex +numList) ' Video Frames'])
    end
    
    ntrials = DS{handles.Index}.ntrials;
    set(handles.ntrials, 'String', ntrials);
    
    
    if ntrials > 1
        numtrials = 10;
        if numtrials > ntrials-handles.trialIndex 
            numtrials = ntrials-handles.trialIndex;
        end
        
        zeroline = zeros(1,nt);
        
        for ii = 1:numtrials
            positionval = [0.106 0.3585+(ii-1)*(0.280612/numtrials)...
                    0.592014 0.280612/numtrials-0.001];
            %subplot('position', [0.106 0.34585+(ii-1)*(0.280612/numtrials)...
            % 0.592014 0.280612/numtrials-0.001]); 
            subplot('position', positionval); 
            if strcmp(sDim, '1-D')
                bar(xlabelrange,handles.respfiles(ii+handles.trialIndex,1:nt),0, 'k');
                hold on; plot(zeroline, 'k');
                hold off;
                p_axis = axis;
                p_axis(2) = nt*ceil(1000/ampsamprate);
                p_axis(3) = 0;
                p_axis(4) = 1.2;
                axis(p_axis)
                axis off;
            else  % Display 2-D data
                %labelrange = handles.frameIndex*binsize+1:(handles.frameIndex+6)*binsize;
                set(handles.displaywhole, 'Visible', 'On');
                if handles.displaywpsth == 1
                    labelrange = 1:nt;
                else
                    labelrange = frameRange(find(frameRange<=nt));
                end
                
                bar(xlabelrange(labelrange),handles.respfiles(ii+handles.trialIndex,labelrange),0, 'k');
                hold on; 
                plot(xlabelrange(labelrange),zeroline(labelrange),'k');
                hold off;
                p_axis = axis;
                p_axis(1) = xlabelrange(labelrange(1));
                p_axis(2) = xlabelrange(labelrange(end));
                p_axis(3) = 0;
                p_axis(4) = 1.2;
                axis(p_axis)
                %xlim([0, nt])
                axis off;
                
            end  % END of sDim trial display
        end
    else
        %axes(handles.spiketrain_axes)
        zeroline = zeros(1,nt);
        positionval = [0.106 0.3685 0.592014 0.280612-0.001];
        %labelrange = handles.frameIndex*binsize+1:(handles.frameIndex+6)*binsize;
        set(handles.displaywhole, 'Visible', 'On');
        if handles.displaywpsth == 1
            labelrange = 1:nt;
        else
            labelrange = frameRange(find(frameRange<=nt));
        end
        
        subplot('position', positionval);
        %plot(xlabelrange(labelrange),handles.respfiles(:,labelrange));
        bar(xlabelrange(labelrange), handles.respfiles(:,labelrange), 0,'k');
        hold on;
        plot(xlabelrange(labelrange), zeroline(labelrange),'k');
        p_axis = axis;
        p_axis(1) = xlabelrange(labelrange(1));
        p_axis(2) = xlabelrange(labelrange(end));
        axis(p_axis)
        axis off;
        ylabel('Spikes')
        hold off
    end

    axes(handles.psth_axes);
    if ntrials > 1
        meanPSTH = mean(handles.respfiles);
    else
        meanPSTH = handles.respfiles;
    end
    
     % Junli: 9/29/04
    % Check if psth_smoothwindow is even or odd and to set cutsize
    if mod(psth_smooth, 2) == 0
        cutsize = 0;
    else
        cutsize = 1;
    end
   
    % Junli: 7/7/2005
    % Also plot psth with time-vary mean rate and with const mean rate
    
    % 1. First calculate time-vary mean rate by calling cal_AVG
    global NBAND outputPath 
    avg_done_check = fullfile(outputPath,'stim_avg.mat');
    if ~exist(avg_done_check)
        cal_AVG(DS, NBAND,1);    % this seems innefficient - how about having this somewhere 
    end
    
    loadpsth = load(fullfile(outputPath, 'stim_avg.mat'), 'Avg_psth', 'constmeanrate'); 
    vary_rate = loadpsth.Avg_psth;
    const_rate = loadpsth.constmeanrate;
    
    %Show smoothed version of psth at window = 15
    halfwinsize = floor(psth_smooth/2);
    wind1 = hanning(psth_smooth)/sum(hanning(psth_smooth)); 
    svagsm=conv(meanPSTH(1:nt),wind1);
    smoothpsth = svagsm(halfwinsize+cutsize:length(svagsm)-halfwinsize);
    
    if strcmp(sDim, '1-D')
        plot(xlabelrange,smoothpsth*ampsamprate,'b'); hold on;
        plot(xlabelrange,(smoothpsth-vary_rate(handles.Index,1:size(smoothpsth,2)))*ampsamprate, 'r'); hold on;
        plot(xlabelrange, (smoothpsth-const_rate)*ampsamprate, 'g'); 
        p_axis = axis;
        p_axis(2) = nt*ceil(1000/ampsamprate);
        axis(p_axis)
    else
        plot(xlabelrange(labelrange),smoothpsth(labelrange)*ampsamprate,'b'); hold on;
        plot(xlabelrange(labelrange),(smoothpsth(labelrange)-vary_rate(handles.Index,1:size(smoothpsth(labelrange),2)))*ampsamprate, 'r'); hold on;
        plot(xlabelrange(labelrange), (smoothpsth(labelrange)-const_rate)*ampsamprate, 'g'); 
        p_axis = axis;        
        p_axis(1) = xlabelrange(labelrange(1));
        p_axis(2) = xlabelrange(labelrange(end));
        
        axis(p_axis)
    end
    ylabel('Smoothed PSTH (spikes/s)')
    %xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
    xlabel('Time (ms)')
     hold off;
    %xlim([0,nt])
    legend('PSTH', 'Time-varying removed', 'Mean removed')
    

% --------------------------------------------------------------------
function varargout = spatialdomain_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
   
    newStr = get(h, 'String');
    if strcmp(newStr, '1-D') | strcmp(newStr, '2-D')
        global sDim
        sDim = newStr;
    else
        errordlg('Please enter valid spatial domain( now 1-D/2-D only).',...
          'Variable Error', 'modal')
        return;
    end  

    global DS;
    checknload(DS,handles.Index,handles);

    % save handles
    guidata(h, handles);
    


% --------------------------------------------------------------------
function varargout = prev5frame_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global DS;

% Get the index pointer and the files 
frameindex = handles.frameIndex;

% update index
i = frameindex - 10;
% If the index is less then one then set it equal to the index of the
% last element in the Addresses array
if i < 0
    i = 0;
end
% Update the index pointer to reflect the new index
handles.frameIndex = i;

% Reset displaywhole option if next5frame has been clicked.
handles.displaywpsth = 0;
guidata(h,handles)

checknload(DS, handles.Index, handles);


% --------------------------------------------------------------------
function varargout = next5frame_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global DS ampsamprate;
nt = min(size(handles.respfiles,2), size(handles.stimfiles, 2));
% Update the index pointer to reflect the new index
range = handles.frameIndex+10;
if range < nt
    handles.frameIndex = range;
else
    handles.frameIndex = 0;
end
% Reset displaywhole option if next5frame has been clicked.
handles.displaywpsth = 0;
guidata(h,handles)

checknload(DS, handles.Index, handles);


% --- Executes during object creation, after setting all properties.
function smoothconst_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor','white');

% --------------------------------------------------------
function smoothconst_Callback(hObject, eventdata, handles)
% --------------------------------------------------------
global psth_smooth

newStr = get(hObject, 'String');
NewVal = str2double(newStr);

% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter postive number.',...
        'Variable Error', 'modal')
    return; 
    
end

psth_smooth = NewVal;
displayall(handles);
     
% --------------------------------------------------------
function smoothpsth_Callback(hObject, eventdata, handles)
% --------------------------------------------------------
helpdlg({'Smoothing PSTH: ',...
      ' ',...
      ' This flag is used for smoothing psth when displaying psth.',...
      ' It is in ms. You can type new number in the text box to modify it.',...
      ' '},...
      'Smooth PSTH Help');


% --- Executes on button press in firsttrials.
function firsttrials_Callback(hObject, eventdata, handles)
    
    % Get the index pointer and the files 
    index = handles.trialIndex;
    
    % update index
    i = index - 10;
    % If the index is less then one then set it equal to the index of the
    % last element in the Addresses array
    if i < 1
        i = 0;
    end
    
    %checknload(DS, i, handles);
    % Update the index pointer to reflect the new index
    handles.trialIndex = i;
    guidata(hObject,handles)
    
    displayall(handles);
    
% --- Executes on button press in secondtrials.
function secondtrials_Callback(hObject, eventdata, handles)
    global DS
    
    % Get the index pointer and the files 
    index = handles.trialIndex;
    
    % update index
    i = index + 10;
    % If the index is less then one then set it equal to the index of the
    % last element in the Addresses array
    if i >= DS{handles.Index}.ntrials
        i = 0;
    end
    
    %checknload(DS, i, handles);
    % Update the index pointer to reflect the new index
    handles.trialIndex = i;
    guidata(hObject,handles)
    
    displayall(handles);


% --- Executes during object creation, after setting all properties.
function ntrilas_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor','white');

function ntrilas_Callback(hObject, eventdata, handles)


% --- Executes on button press in displaywhole.
function displaywhole_Callback(hObject, eventdata, handles)
handles.displaywpsth = 1;
guidata(hObject, handles);
displayall(handles);
% --------------------------------------------------------------------
function varargout = help_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    ttlStr = get(handles.figure1, 'Name');
    hlpStr = [...
            '                                                               '
            '   DISPLAY INPUT WINDOW                                        '
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
            '     2. Plot of spike train (ntrials X TIME)                   '
            '        Note: We only display 10 trials here. If you have more '
            '        than 10 trials, you can click next 10 trials to see the'
            '        rest trials.                                           '
            '     3. Plot of psth (TIME)                                    '
            '        PSTH is computed as average across all the trials and  '
            '        is shown in blue line. The green line shows the psth   '
            '        removed constant mean rate. The red lines shows the    '
            '        psth removed time-varying mean rate.                  '
            '                                                               '
            '     Spatial Domain: shows the dimension size of stimulus.     '
            '     Stim files: shows the filename of the stimulus.           '
            '     Resp files: shows the filename of the response.           '
            '     Prev_files: click this button to see previous data files. '
            '     Next_files: see next data files visually.                 '
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
            ' Updated by Junli, Sept 2005.                                  '
            '                                                               '];

    myFig = handles.figure1;
    helpwin(hlpStr, ttlStr);

% --------------------------------------------------------------------
%  END of displayraw_preprocessed_GUI
% --------------------------------------------------------------------
