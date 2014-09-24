% --------------------------------------------------------------------
function varargout = displayrawstim_GUI(varargin)
% --------------------------------------------------------------------
%
% DISPLAYRAWSTIM_GUI Application M-file for displayrawstim_GUI.fig
%    FIG = DISPLAYRAWSTIM_GUI launch displayrawstim_GUI GUI.
%    DISPLAYRAWSTIM_GUI('callback_name', ...) invoke the named callback.
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

    fig = openfig(mfilename);

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
    % load the first data file for default display window
    global rawData
    if rawData == 0
        tt=errordlg('This option is only for raw stimulus.',...
            'Raw Stimulus only','modal');
        uiwait(tt);
        delete(handles.figure1);
        return;
    end
    
    handles.frameIndex = 1;
    handles.trialIndex = 0;
    handles.displaywpsth = 0;
    handles.stopflag = 0;
    handles.Index = 1;    % Data Set index
    handles.numList = 6;
    set(handles.DataSet, 'string', num2str(handles.Index));
    guidata(fig, handles);
       
    global rawDS ampsamprate respsamprate 
    %set(handles.datasetSlider, 'min', 1);
    if length(rawDS) ==1
        set(handles.datasetSlider, 'min', 0, 'max', length(rawDS), 'value',0,'sliderstep', [1 1]);
    elseif length(rawDS) > 1
        set(handles.datasetSlider, 'Min', 1, 'Max', length(rawDS),'value',1,...
            'sliderstep', [1/(length(rawDS)-1) 1/(length(rawDS)-1)]);
    else
        warndlg('No input selected yet. Exiting...')
        return;
    end
    
    if isempty(ampsamprate) || isempty(respsamprate)
        ampsamprate = 1000;
        respsamprate = 1000;
    end
    set(handles.frameSet, 'String', 1);
    checknload(rawDS, 1, handles);
    
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
    global rawDS;

    % Get the index pointer and the files 
    index = handles.Index;
    
    % update index
    i = index + 1;
    % If the index is less then one then set it equal to the index of the
    % last element in the Addresses array
    if i > length(rawDS) 
        i = 1;
    end

    
    % Update the index pointer to reflect the new index
    handles.Index = i;
    % Also reset frameIndex
    handles.frameIndex = 0;
    guidata(h,handles)

    checknload(rawDS,i,handles);



% --------------------------------------------------------------------
function varargout = close_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    delete(handles.figure1);


% --------------------------------------------------------------------
function varargout = prev_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    global rawDS;

    % Get the index pointer and the files 
    index = handles.Index;
    
    % update index
    i = index - 1;
    % If the index is less then one then set it equal to the index of the
    % last element in the Addresses array
    if i < 1
        i = length(rawDS);
    end
    
    
    % Update the index pointer to reflect the new index
    handles.Index = i;
    % Also reset frameIndex
    handles.frameIndex = 0;
    guidata(h,handles)
    
    checknload(rawDS, i, handles);
    


% --------------------------------------------------------------------
function pass = checknload(rawDS, n, handles)
% --------------------------------------------------------------------
    
    % Get data's extension from the filename
    [path,name,ext,ver] = fileparts(rawDS{n}.stimfiles);
    set(handles.show_stimfile_text, 'String', [name, ext]);

    % Response data
%     [rpath, rname, rext, rver] = fileparts(rawDS{n}.respfiles);
%     set(handles.show_respfile_text, 'String', [rname, rext]);

    % 1. Check if the data is raw data and preprocessed data.
    global rawData
    if rawData == 0    % is preprocessed data

    else          % is raw data
        % Check if the data file is song
        switch ext
            case {'.dat', '.txt', '.mat'}  % 2-D raw data (movie)
                stimstim = load(rawDS{n}.stimfiles);
                
                % Check if stimstim is struct
                if isstruct(stimstim)
                    flds = fieldnames(stimstim);
                    if (length(flds) == 1)
                        stimstim = getfield(stimstim, char(flds{1}));

                    end
                end
                
                % Check if 2-D movie data or 1-D song data
                global sDim NBAND
                if length(size(stimstim))>2  % 2-D movie data
                    [xsize, ysize, framesize] = size(stimstim);
                    sDim = '2-D';
                    NBAND = xsize*ysize;
                    set(handles.spatialdomain, 'String', [sDim, '(Movie)']);
                else      % 1-D song data
                    warndlg({'The data is identified as auditory data and but make sure the data ',...
                        'is raw data with .wav extension.'}, '!!WARNING!!','modal')
                    sDim = '1-D';
                    [NBAND framesize] = size(stimstim);
                    set(handles.spatialdomain, 'String', [sDim, '(Sound)']);
                end
                    
                % Initial starting point and ending point
                handles.startFrame = 1;
                handles.endFrame = framesize;
                handles.Data = stimstim;
                numList = handles.numList;
                if framesize ==1
                    set(handles.frameSetSlider, 'min', 0, 'max', 1, 'value',0,'sliderstep', [1 1]);
                else
                    set(handles.frameSetSlider, 'Min', 1, 'Max', framesize,'value',1,...
                        'sliderstep', [numList/(framesize-1) numList/(framesize-1)]);
                end
                guidata(handles.figure1, handles);

                displayraw(handles);

            case {'.wav'}

                [Data,Freq,Bits] = wavread(rawDS{n}.stimfiles);
                Data = Data./max(max(abs(Data)));
                DataSize = length(Data);

                handles.Data = Data;
                global sDim stimsamprate
                sDim = '1-D';
                set(handles.spatialdomain, 'String', [sDim, '(Sound)']);
                stimsamprate = Freq;

                % Initial starting point and ending point
                handles.startFrame = 1;
                handles.endFrame = DataSize;
                guidata(handles.figure1, handles);
                

                displayraw(handles);   % Display raw stimulis file
        end

    end
    set(handles.framestart, 'String', num2str(handles.startFrame));
    set(handles.frameend, 'String', num2str(handles.endFrame));
    % Display response data: including spike train and psth
    %displayall(handles);

% --------------------------------------------------------------------
function displayraw(handles)
% --------------------------------------------------------------------

% Draw x-axes indicator
axes(handles.indicator);
Data = handles.Data;
DataSize = length(Data);
ttH=plot(ones(1,DataSize));  box on;
patch([0, DataSize, DataSize, 0], [ 0, 0, 2,2], 'y');
axis tight;
set(get(ttH,'Parent'),'YTickLabel',[]);
set(get(ttH,'Parent'),'YTick',[]);
ylabel({'  data  '; 'selector'})

global stimsamprate sDim
if ~isempty(stimsamprate)
    if sDim == '2-D'
        timebin = round(1000/stimsamprate);  % in ms
    end
    %         xlabel(['Time (in ', num2str(timebin), ' ms)']);
    set(handles.stimsamprateD, 'String', stimsamprate);
    ss = get(get(ttH,'Parent'),'XTick');
    ss = num2str(round(ss'*1000/stimsamprate));
    set(get(ttH,'Parent'),'XTickLabel', ss);
    xlabel('Time (in ms)');
end
% Show starting position and ending position

% draw stimulus based on sDim
if sDim == '1-D'
    set(handles.stimsamprateD, 'String', num2str(stimsamprate));

    axes(handles.stim_axes);
    selectRange = handles.startFrame:handles.endFrame;
    clear ttH
    ttH = plot(Data(selectRange), 'r-');
    axis tight
    ss = get(get(ttH,'Parent'),'XTick');
    if ~isempty(stimsamprate)
        ss = num2str(round(ss'*1000/stimsamprate));
        set(get(ttH,'Parent'),'XTickLabel', ss);
        xlabel('Time (in ms)');
    else
        xlabel('Points');
    end
    box on; grid on;

%     % Save for play /stop button
%     player = audioplayer(Data, stimsamprate);
%     handles.DataSize = DataSize;
%     handles.player = player;
    guidata(handles.figure1, handles);
else

    [xsize, ysize, DataSize] = size(handles.Data);
    % Set up maximum numbers of frames can be shown in the window
    numList = handles.numList;  % Right now only 6 frames for one screen
    % Set up frameSet slider

    frameRange = handles.frameIndex:handles.frameIndex+numList;
    %frameRange = handles.startFrame:handles.endFrame;
    if length(find(frameRange>DataSize)) == 0
        whole = reshape(Data(:,:,frameRange), xsize, ysize*length(frameRange));
    else
        whole = reshape(Data(:,:,frameRange(find(frameRange<=DataSize))),...
            xsize, ysize*length(frameRange(find(frameRange<=DataSize))));
    end

    axes(handles.stim_axes);
    ttH = imagesc(whole);
    axis image
    set(get(ttH,'Parent'),'YTickLabel',[]);
    set(get(ttH,'Parent'),'YTick',[]);
    set(get(ttH,'Parent'),'XTickLabel',[]);
    set(get(ttH,'Parent'),'XTick',[]);
    colormap(gray)
    if isempty(stimsamprate)
        title([num2str(frameRange(1)) '  to '...
            num2str(frameRange(end)) ' Video Frames'])
    else
        title([num2str(frameRange(1)) '  to '...
            num2str(frameRange(end)) ' Video Frames (per\_frame\_dur = ',...
            num2str(timebin), ' ms )'])
    end

end

% --------------------------------------------------------------------
function displayrawbyframe(handles)
% --------------------------------------------------------------------

% Show starting position and ending position
set(handles.framestart, 'String', num2str(handles.startFrame));
set(handles.frameend, 'String', num2str(handles.endFrame));
Data = handles.Data;
global stimsamprate sDim
if ~isempty(stimsamprate)
    if sDim == '2-D'
        timebin = round(1000/stimsamprate);  % in ms
        xlabel(['Time (in ', num2str(timebin), ' ms)']);
        set(handles.stimsamprateD, 'String', stimsamprate);
    else
        ss = get(get(ttH,'Parent'),'XTick');
        ss = num2str(round(ss'*1000/stimsamprate));
        set(get(ttH,'Parent'),'XTickLabel', ss);
        xlabel('Time (in ms)');
    end
end
% draw stimulus based on sDim
if sDim == '1-D'
    set(handles.stimsamprateD, 'String', num2str(stimsamprate));

    axes(handles.stim_axes);
    selectRange = handles.startFrame:handles.endFrame;
    plot(Data(selectRange), 'r-');
    axis tight; box on; grid on;

    guidata(handles.figure1, handles);
else

    [xsize, ysize, DataSize] = size(handles.Data);
    % Set up maximum numbers of frames can be shown in the window
    numList = handles.numList;  % Right now only 6 frames for one screen
    % Set up frameSet slider

    %frameRange = handles.frameIndex+1:handles.frameIndex+numList;
    frameRange = handles.startFrame:handles.endFrame;
    %frameRange = handles.startFrame:handles.endFrame;
    if length(find(frameRange>DataSize)) == 0
        whole = reshape(Data(:,:,frameRange), xsize, ysize*length(frameRange));
    else
        whole = reshape(Data(:,:,frameRange(find(frameRange<=DataSize))),...
            xsize, ysize*length(frameRange(find(frameRange<=DataSize))));
    end

    axes(handles.stim_axes);
    ttH = imagesc(whole);
    %axis image
    set(get(ttH,'Parent'),'YTickLabel',[]);
    set(get(ttH,'Parent'),'YTick',[]);
    set(get(ttH,'Parent'),'XTickLabel',[]);
    set(get(ttH,'Parent'),'XTick',[]);
    colormap(gray)
    global stimsamprate
    if isempty(stimsamprate)
        title([num2str(frameRange(1)) '  to '...
            num2str(frameRange(end)) ' Video Frames'])
    else
        title([num2str(frameRange(1)) '  to '...
            num2str(frameRange(end)) ' Video Frames (per\_frame\_dur = ',...
            num2str(timebin), ' ms )'])
    end
   

end

% --------------------------------------------------------------------
function displayall(handles)
% --------------------------------------------------------------------
% Now display input stimulus file
global ampsamprate psth_smooth rawDS
if isempty(psth_smooth)
    psth_smooth = 21;
end
global sDim rawDS
set(handles.spatialdomain, 'String', sDim);
set(handles.smoothconst,'String',psth_smooth);


% Response data

n = handles.Index;
[rpath, rname, rext, rver] = fileparts(rawDS{n}.respfiles);
set(handles.show_respfile_text, 'String', [rname, rext]);

switch rext
    case {'.dat', '.txt', ''}
        handles.respfiles = load(rawDS{n}.respfiles);

    case {'.mat'}
        respMat = load(rawDS{n}.respfiles);

        % Validate the MAT-file
        flds = fieldnames(respMat);

        % Check if response data is in spike arrival time
        % or already preprocessed.
        if (length(flds) == 1)
            rawResp = getfield(respMat, char(flds{1}));
            if iscell(rawResp)
                spiketrain = zeros(rawDS{n}.ntrials,rawDS{n}.nlen);
                for trial_ind =1:rawDS{n}.ntrials

                    spiketrain(trial_ind, rawResp{trial_ind}) = ones(1,length(rawResp{trial_ind}));
                end
                newpsth = resample(spiketrain', ampsamprate, respsamprate);

                newpsth = newpsth'; % make sure new response data is trials x T.
                newpsth(find(newpsth < 0)) = 0;
            else
                newpsth = rawResp;
            end
            [xs, ys] = size(newpsth);
            if xs > ys
                newpsth = newpsth';
            end
            rgoodidx = find(~isnan(newpsth));
            newpsth = newpsth(rgoodidx);
            handles.respfiles = newpsth;

        end
    otherwise
        errordlg(['Only ASCII and MAT binary fileformats work',...
            ' in this version.'], 'File Type Error', 'modal')
        return;
end

spikeSize = size(handles.respfiles);
rawDS{handles.Index}.ntrials = min(spikeSize);
ntrials = rawDS{handles.Index}.ntrials;
set(handles.ntrials, 'String', ntrials);

nt = min(length(handles.Data), length(newpsth));
xlabelrange = (1:nt)*ceil(1000/ampsamprate);
numList = handles.numList;

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

    frameRange = handles.frameIndex+1:handles.frameIndex+numList;
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

rgoodidx = find(~isnan(meanPSTH));
meanPSTH = meanPSTH(rgoodidx);
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

global rawData
if rawData ==0
    global NBAND outputPath
    avg_done_check = fullfile(outputPath,'stim_avg.mat');
    if ~exist(avg_done_check)
        cal_AVG(rawDS, NBAND,1);    % this seems innefficient - how about having this somewhere
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
else
    %Show smoothed version of psth at window = 15
    halfwinsize = floor(psth_smooth/2);
    wind1 = hanning(psth_smooth)/sum(hanning(psth_smooth));
    svagsm=conv(meanPSTH(1:nt),wind1);
    smoothpsth = svagsm(halfwinsize+cutsize:length(svagsm)-halfwinsize);
    if strcmp(sDim, '1-D')
        plot(xlabelrange,smoothpsth*ampsamprate,'b'); 
       
        p_axis = axis;
        p_axis(2) = nt*ceil(1000/ampsamprate);
        axis(p_axis)
    else
        plot(xlabelrange(labelrange),smoothpsth(labelrange)*ampsamprate,'b'); 
        
        p_axis = axis;
        p_axis(1) = xlabelrange(labelrange(1));
        p_axis(2) = xlabelrange(labelrange(end));

        axis(p_axis)
    end
    ylabel('Smoothed PSTH (spikes/s)')
    %xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
    xlabel('Time (ms)')
    
end

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

global rawDS;
checknload(rawDS,handles.Index,handles);

% save handles
guidata(h, handles);



% --------------------------------------------------------------------
function varargout = prev5frame_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global rawDS;

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

checknload(rawDS, handles.Index, handles);


% --------------------------------------------------------------------
function varargout = next5frame_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
global rawDS ampsamprate;
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

checknload(rawDS, handles.Index, handles);


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
    
    %checknload(rawDS, i, handles);
    % Update the index pointer to reflect the new index
    handles.trialIndex = i;
    guidata(hObject,handles)
    
    displayall(handles);
    
% --- Executes on button press in secondtrials.
function secondtrials_Callback(hObject, eventdata, handles)
    global rawDS
    
    % Get the index pointer and the files 
    index = handles.trialIndex;
    
    % update index
    i = index + 10;
    % If the index is less then one then set it equal to the index of the
    % last element in the Addresses array
    if i >= rawDS{handles.Index}.ntrials
        i = 0;
    end
    
    %checknload(rawDS, i, handles);
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
            '   STRFPAK: Display raw stimulus window                        '
            '                                                               '
            '  Right Panel:                                                 '
            '        Spatial dimension: STRFPAK identifies the selected data'
            '        dimensionality. If the data is raw data:               '
            '             If length(size(stim)) == 2                        '
            '                 Spatial Dimension = 1-D                       '
            '             else if length(size(stim)) > 2                    '
            '                 Spatial Dimension = 2-D                       '
            '        Data Set:  the text box shows the number of data sets  '
            '             which are displaying in the left panel. The user  '
            '             can either click the slider or type in the text   '
            '             box to display the desired stimulus.              '
            '        Frame Set: this option is only for 2-D. 2-D data are   '
            '             collection of the movie frames. The current window'
            '             only shows six frames one time. The user can      '
            '             either click the slider or type in the text box   '
            '             to display the desired 6 frames.                  '
            '        Stim Filename: shows the filename of the stimulus which'
            '             is displaying on the left panel.                  '
            '        Stim SampRate: shows the sampling rate of the stimulus '
            '             If STRFPAK can identify, the box shows the value. '
            '             Otherwise, it shows nothing. The user can type in '
            '             if the user knows the stim sampling rate. For     '
            '             example, the stim samp rate for demo movie data is'
            '             60 Hz and for auditory sound is 32000 Hz.         '
            '   Left Panel:                                                 '
            '                                                               '
            '  Start position:  shows the starting point of the display. By '
            '        default, it is 1. It is changed based on Select Data   '
            '        button and frame Set.                                  '
            '  End position: shows the ending point of the display. By      '
            '        default, it is total length (in point) of the stimulus.'
            '        It is changed based on Select Data button and frame Set'
            '                                                               '
            '   Select Data: When this button is clicked, it uses the mouse '
            '        click to select two position from the yellow data      '
            '        selector patch.  So the user need click on the yellow  '
            '        area two times to get both start and end positions.    '
            '                                                               '
            '   Data Selector: the yellow area is used for choosing the new '
            '        start and end positions.  If the stim samp rate is     '
            '        specified, the label of the data selector shows in time'
            '        Otherwise, it only shows the points.                   '
            '   Stim Axes: display the amplitude of the stimulus for 1-D    '
            '        data or movie frames for 2-D.                          '
            '   Play:  play the sound for 1-D data or play the movie frame  '
            '        for 2-D data.                                          '
            '   Stop: stop playing the sound or click the movie frames to   '
            '        stop movie playing.                                    '
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
            ' Updated by Junli, May 2006.                                   '
            '                                                               '];

    myFig = handles.figure1;
    helpwin(hlpStr, ttlStr);

% --------------------------------------------------------------------
%  END of displayrawstim_GUI
% --------------------------------------------------------------------


% --- Executes on slider movement.
function datasetSlider_Callback(hObject, eventdata, handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
     newSlider = round(get(hObject, 'Value'));
     if newSlider <= 0 || round(newSlider) == 0
         set(handles.DataSet, 'String', 1);
     else
        set(handles.DataSet, 'String', round(newSlider));
        handles.Index = round(newSlider);
     end
     guidata(handles.figure1, handles);
     global rawDS
     checknload(rawDS, handles.Index, handles);

% --- Executes during object creation, after setting all properties.
function datasetSlider_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function DataSet_Callback(hObject, eventdata, handles)
     global rawDS
     maxlength = length(rawDS);
     newIndex = str2double(get(hObject,'String'));
     if newIndex <=0 || newIndex > maxlength
         tt=errordlg(['DataSet is out of range. It need to be between 1 and ', num2str(maxlength)]);
         uiwait(tt);
         set(hObject, 'String', num2str(handles.Index));
     else
         handles.Index = newIndex;
         %set(hObject, 'String', num2str(handles.Index));
         % Update slider value based on dataSet value
         set(handles.datasetSlider, 'String', newIndex);
     end
     guidata(handles.figure1, handles);
     checknload(rawDS, handles.Index, handles);
     
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


% --- Executes on slider movement.
function trialSlider_Callback(hObject, eventdata, handles)
% hObject    handle to trialSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function trialSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trialSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
global exitLoop
global sDim
if sDim == '1-D'
     % Save for play /stop button
    global stimsamprate
    data = handles.Data;
    selecteddata = data(handles.startFrame:handles.endFrame, :);
    handles.player = audioplayer(selecteddata, stimsamprate);
    guidata(handles.figure1, handles);
    play(handles.player);
elseif sDim == '2-D'
    set(handles.stop, 'visible', 'off');
    set(handles.stopmov_flag, 'visible', 'on');
    mov = handles.Data;
    axes(handles.stim_axes);
    %set(handles.stim_axes, 'ButtonDownFcn', 'break');
    exitLoop=0;
    for ii=handles.startFrame:handles.endFrame
        hh=imagesc(mov(:,:,ii)); 
        axis off;
        axis image;
        title(ii);
        set(hh,'ButtonDownFcn','setExitLoop'); 
        pause(0.05);
        if exitLoop == 1
            break;
        end
    end
    set(handles.stop, 'visible', 'on');
    set(handles.stopmov_flag, 'visible', 'off');
end

function setExitLoop()
global exitLoop
exitLoop=1;
   
    
% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
global sDim
if sDim == '1-D'
    stop(handles.player);
else
   warndlg('Please play the movie first.', 'color', 'r');
end



function framestart_Callback(hObject, eventdata, handles)

Data = handles.Data;
DataSize = length(Data);
newvalue = str2double(get(hObject,'String'));
if newvalue <=0 || newvalue > DataSize
    warndlg('Starting point is out of range.');
    newvalue = 1;
elseif newvalue >= handles.endFrame
    warndlg('Starting point must be less than ending point.');
    newvalue = handles.endFrame-1;
end
handles.startFrame = newvalue;
set(handles.framestart, 'String', num2str(handles.startFrame));

% draw starting frame in indicator
% Draw x-axes indicator
axes(handles.indicator);
ttH=plot(ones(1,DataSize));  box on;
patch([0, DataSize, DataSize, 0], [ 0, 0, 2,2], 'y');
axis tight;
set(get(ttH,'Parent'),'YTickLabel',[]);
set(get(ttH,'Parent'),'YTick',[]);
ylabel({'  data  '; 'selector'})

global stimsamprate sDim
if ~isempty(stimsamprate)
    if sDim == '2-D'
        timebin = round(1000/stimsamprate);  % in ms
        xlabel(['Time (in ', num2str(timebin), ' ms)']);
        set(handles.stimsamprateD, 'String', stimsamprate);
    else
        ss = get(get(ttH,'Parent'),'XTick');
        ss = num2str(round(ss'*1000/stimsamprate));
        set(get(ttH,'Parent'),'XTickLabel', ss);
        xlabel('Time (in ms)');
    end
end
hold on

% draw starting position in indicator
plot([newvalue newvalue], [0, 2], 'r-');
text((newvalue-DataSize/100),(2-2/4),num2str(newvalue),'FontWeight','bold');
plot([handles.endFrame handles.endFrame], [0,2],'r-');
text((handles.endFrame-5*DataSize/100),(2-2/4),num2str(handles.endFrame),'FontWeight','bold');
hold off
guidata(handles.figure1, handles);
displayrawbyframe(handles);
    

% --- Executes during object creation, after setting all properties.
function framestart_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frameend_Callback(hObject, eventdata, handles)
DataSize = length(handles.Data);
newvalue = str2double(get(hObject,'String'));
if newvalue <=0 || newvalue > length(handles.Data)
    warndlg('Ending point is out of range.');
    newvalue = length(handles.Data);
elseif newvalue <= handles.startFrame
    warndlg('Ending point must be greater than ending point.');
    newvalue = handles.startFrame-1;
end
handles.endFrame = newvalue;
set(handles.frameend, 'String', num2str(handles.endFrame));
guidata(handles.figure1, handles);
% draw starting frame in indicator
% Draw x-axes indicator
axes(handles.indicator);
ttH=plot(ones(1,DataSize));  box on;
patch([0, DataSize, DataSize, 0], [ 0, 0, 2,2], 'y');
axis tight;
set(get(ttH,'Parent'),'YTickLabel',[]);
set(get(ttH,'Parent'),'YTick',[]);
ylabel({'  data  '; 'selector'})

global stimsamprate sDim
if ~isempty(stimsamprate)
    if sDim == '2-D'
        timebin = round(1000/stimsamprate);  % in ms
        xlabel(['Time (in ', num2str(timebin), ' ms)']);
        set(handles.stimsamprateD, 'String', stimsamprate);
    else
        ss = get(get(ttH,'Parent'),'XTick');
        ss = num2str(round(ss'*1000/stimsamprate));
        set(get(ttH,'Parent'),'XTickLabel', ss);
        xlabel('Time (in ms)');
    end
end
hold on

% draw starting position in indicator
plot([handles.startFrame handles.startFrame], [0,2], 'r-');
text((handles.startFrame-DataSize/100),(2-2/4),num2str(handles.startFrame),'FontWeight','bold');
plot([newvalue newvalue], [0, 2], 'r-');
text((newvalue-5*DataSize/100),(2-2/4),num2str(newvalue),'FontWeight','bold');
hold off
displayrawbyframe(handles);
     

% --- Executes during object creation, after setting all properties.
function frameend_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stimsamprateD_Callback(hObject, eventdata, handles)

   newStrValue = get(hObject, 'String');
   NewVal = str2double(newStrValue);

   % Check if the entered value falls within the allowable range
   if (NewVal< 0) | ~isnumeric(NewVal)
       % Set default value for step
       errordlg('Not valid value for stim sampling rate.');
       return;
   elseif isnan(NewVal)
       clear global stimsamprate;
       set(handles.stimsamprateD, 'String', '');
       displayraw(handles);
       return;
   end

   % Assign the global variable 
   global stimsamprate rawDS
   stimsamprate = NewVal;
   set(handles.stimsamprateD, 'String', num2str(stimsamprate));
   %checknload(rawDS, handles.Index,handles);
   displayraw(handles);


% --- Executes during object creation, after setting all properties.
function stimsamprateD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimsamprateD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rawdataDisplay.
function rawdataDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to rawdataDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rawdataDisplay


% --- Executes on button press in ppdataDisplay.
function ppdataDisplay_Callback(hObject, eventdata, handles)

% --- Executes on button press in selectdata.
function selectdata_Callback(hObject, eventdata, handles)
indicatorHandles = handles.indicator;
axes(indicatorHandles);
Data = handles.Data;
DataSize = length(Data);
%[xsize, ysize,DataSize] = size(Data); 

ttH=plot(ones(1,DataSize));  box on;
patch([0, DataSize, DataSize, 0], [ 0, 0, 2,2], 'y');
axis tight;
set(get(ttH,'Parent'),'YTickLabel',[]);
set(get(ttH,'Parent'),'YTick',[]);
ylabel({'  data  '; 'selector'})

global stimsamprate sDim
if ~isempty(stimsamprate)
    if sDim == '2-D'
        timebin = round(1000/stimsamprate);  % in ms
    end
    %         xlabel(['Time (in ', num2str(timebin), ' ms)']);
    set(handles.stimsamprateD, 'String', stimsamprate);
    ss = get(get(ttH,'Parent'),'XTick');
    ss = num2str(round(ss'*1000/stimsamprate));
    set(get(ttH,'Parent'),'XTickLabel', ss);
    xlabel('Time (in ms)');
end
hold on;

exitclick = 0;
Count = 1;
xmin = 1; 
LowerLim = xmin;
xmax = DataSize;
ymin = 0;
ymax = 2;
HandleText = text((xmax-xmin)/10,(ymin + (ymax-ymin)/4),...
    'Out of Range! Please Re-Select!','tag','outofrange');
set(HandleText,'color',[0.5 0 0.6],'visible','off','FontWeight','bold');

% -----------------------------------------------------------------
% Select Data by Selecting an Interval via Mouse Click on the Graph
% -----------------------------------------------------------------
while (exitclick ~= 1);
    MouseP = ginput(1);
    % First check whether it is right handles
    if gca == indicatorHandles
        xCoord = round(MouseP(1));

        if (xCoord <= LowerLim);
            set(HandleText,'visible','on');
        elseif ((xCoord >= xmax) & (Count == 1));
            set(HandleText,'visible','on');
        elseif (xCoord >= xmax);
            xCoord = xmax;
            set(HandleText,'visible','off');
            Select(Count) = xCoord;
            plot([xCoord,xCoord],[ymin,ymax],'r-');
            if (Count == 1);
                text((xCoord + (xmax-xmin)/100),(ymax - (ymax-ymin)/4),num2str(xCoord),'FontWeight','bold');
            else
                text((xCoord - 5*(xmax-xmin)/100),(ymax - (ymax-ymin)/4),num2str(xCoord),'FontWeight','bold');
            end;
            exitclick = 1;
        elseif (Count >= 2);
            set(HandleText,'visible','off');
            Select(Count) = xCoord;
            plot([xCoord,xCoord],[ymin,ymax],'r-');
            if (Count == 1);
                text((xCoord + (xmax-xmin)/100),(ymax - (ymax-ymin)/4),num2str(xCoord),'FontWeight','bold');
            else
                text((xCoord - 5*(xmax-xmin)/100),(ymax - (ymax-ymin)/4),num2str(xCoord),'FontWeight','bold');
            end;
            exitclick = 1;
        else
            set(HandleText,'visible','off');
            Select(Count) = xCoord;
            plot([xCoord,xCoord],[ymin,ymax],'r-');
            if (Count == 1);
                text((xCoord + (xmax-xmin)/100),(ymax - (ymax-ymin)/4),num2str(xCoord),'FontWeight','bold');
            else
                text((xCoord - 5*(xmax-xmin)/100),(ymax - (ymax-ymin)/4),num2str(xCoord),'FontWeight','bold');
            end;
            Count = Count + 1;
            %LowerLim = xCoord;
        end;
    else  % incorrect handles
        set(HandleText,'visible','on');
    end
end;
hold off
handles.startFrame = min(Select);
handles.endFrame = max(Select);
set(handles.framestart, 'String', num2str(handles.startFrame));
set(handles.frameend, 'String', num2str(handles.endFrame));
guidata(handles.figure1, handles);

axes(handles.stim_axes);
selectRange = handles.startFrame:handles.endFrame;
global sDim
if sDim == '1-D'
    plot(selectRange,Data(selectRange), 'r-');
    axis tight; box on; grid on;
elseif sDim == '2-D'
    [xsize, ysize, DataSize] = size(Data);
    frameRange = handles.startFrame:handles.endFrame;
    %frameRange = handles.startFrame:handles.endFrame;
    if length(find(frameRange>DataSize)) == 0
        whole = reshape(Data(:,:,frameRange), xsize, ysize*length(frameRange));
    else
        whole = reshape(Data(:,:,frameRange(find(frameRange<=DataSize))),...
            xsize, ysize*length(frameRange(find(frameRange<=DataSize))));
    end

    ttH = imagesc(whole);
    %axis image
    set(get(ttH,'Parent'),'YTickLabel',[]);
    set(get(ttH,'Parent'),'YTick',[]);
    set(get(ttH,'Parent'),'XTickLabel',[]);
    set(get(ttH,'Parent'),'XTick',[]);
    colormap(gray)
    if isempty(stimsamprate)
        title([num2str(frameRange(1)) '  to '...
            num2str(frameRange(end)) ' Video Frames'])
    else
        title([num2str(frameRange(1)) '  to '...
            num2str(frameRange(end)) ' Video Frames (per\_frame\_dur = ',...
            num2str(timebin), ' ms )'])
    end
   
end


function frameSet_Callback(hObject, eventdata, handles)
   
global sDim
if strcmp(sDim, '2-D')
    
   maxlength = length(handles.Data);
   newIndex = str2double(get(hObject,'String'));
   if newIndex <=0 || newIndex >= maxlength
       warndlg(['FrameSet is out of range. It need to be between 1 and ', num2str(maxlength)]);
       return;
       
   else
       handles.frameIndex = newIndex;
       %set(hObject, 'String', num2str(handles.Index));
       % Update slider value based on dataSet value
       set(handles.frameSetSlider, 'value', newIndex);
   end
   set(handles.frameSet, 'String', handles.frameIndex);
   % Show starting position and ending position
   handles.startFrame = handles.frameIndex;
   handles.endFrame = handles.frameIndex+handles.numList;
  set(handles.framestart, 'String', handles.frameIndex);
  set(handles.frameend, 'String', handles.frameIndex+handles.numList);
   guidata(handles.figure1, handles);
   displayraw(handles);
   %displayall(handles);
end

% --- Executes during object creation, after setting all properties.
function frameSet_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function frameSetSlider_Callback(hObject, eventdata, handles)
global sDim
if strcmp(sDim, '2-D')
    maxlength = length(handles.Data);
    newSlider = round(get(hObject, 'Value'));
    if newSlider <= 0
        set(handles.frameSet, 'String', 1);
    elseif newSlider >= maxlength
        set(handles.frameSet', 'String', maxlength);
    else
        set(handles.frameSet, 'String', round(newSlider));
        handles.frameIndex = round(newSlider);
    end
    handles.startFrame = handles.frameIndex;
    handles.endFrame = handles.frameIndex+handles.numList;
    set(handles.framestart, 'String', handles.frameIndex);
    set(handles.frameend, 'String', handles.frameIndex+6);
    guidata(handles.figure1, handles);
    displayraw(handles);
    %displayall(handles);
end

% --- Executes during object creation, after setting all properties.
function frameSetSlider_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function stimsamprateUnit_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

