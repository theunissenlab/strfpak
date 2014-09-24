 % --------------------------------------------------------------------
function varargout = displayinput2_GUI(varargin)
% --------------------------------------------------------------------
%
% DISPLAYINPUT_GUI Application M-file for displayinput_GUI.fig
%    FIG = DISPLAYINPUT_GUI launch displayinput_GUI GUI.
%    DISPLAYINPUT_GUI('callback_name', ...) invoke the named callback.
% 
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

    % for resize property
    set(fig, 'resize', 'on');
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
           hUIControls],'units','normalized','fontunits','normalized');

    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    handles.frameIndex = 0;
    handles.trialindex = 1;
    guidata(fig, handles);
       
    % load the first data file for default display window
    global DS
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

    checknload(DS,i,handles);

    % Update the index pointer to reflect the new index
    handles.Index = i;
    guidata(h,handles)



% --------------------------------------------------------------------
function varargout = help_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    ttlStr = get(handles.figure1, 'Name');
    hlpStr = [...
            '                                                               '
            'The DISPLAY INPUT WINDOW shows 3 graphes of the data set:      '
            '                                                               '
            '     1. Plot of stimuli file:                                  '
            '             X-axis is TIME (in seconds);                      '
            '             Y-axis is Spatial domain (will be frequency in Hz '
            '             if spatio-domain is 1-Dimension and wrapped space '
            '             axes if 2-Dimension).                             '
            '     2. Plot of spike train (ntrials X TIME)                   '
            '        Note: We only display 10 trials here. If we have more  '
            '        than 10 trials, we display the first 10 trials.        '
            '     3. Plot of psth (TIME)                                    '
            '        Note: If we have one trail spike train, we take its    '
            '        psth as itself. Otherwise, we compute avg of ntrials   '
            '        as their psth.                                         '
            '                                                               '
            'You can choose prev and next to display different data sets.   '
            '                                                               '
            'You can also see which data set you are displaying from stim   '
            'file text field and response file text field.                  '
            '                                                               '
            'Zoom in/out properties:                                        '
            '   MATLAB has built-in function for zooming the figure. If you '
            '   want to zoom in/out the figure, you go to VIEW menu and     '
            '   then choose FIGURE TOOLBAR. Then you choose ZOOM IN/OUT     '
            '   button and select the part of the figure for zooming.       ' 
            '                                                               '];

    myFig = handles.figure1;
    helpwin(hlpStr, ttlStr);



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
    
    checknload(DS, i, handles);
    % Update the index pointer to reflect the new index
    handles.Index = i;
    guidata(h,handles)


% --------------------------------------------------------------------
function pass = checknload(DS, n, handles)
% --------------------------------------------------------------------
    
    % Check whether we have data files selected.
    if length(DS) == 0 
        errordlg('There are no data files. Please select data files first.',...
                'Data Files Error', 'modal')
        return;
    elseif n > length(DS) | (n < 0)
        errordlg('Index out of DS range.', 'File index wrong', 'modal')
        return;
    end 
    
    %load stimulus file
    [path,name,ext,ver] = fileparts(DS{n}.stimfiles);
    switch ext
    case {'.dat', '.txt', ''}
        handles.stimfiles = load (DS{n}.stimfiles);
        guidata(handles.figure1, handles)
    case {'.mat'}
        stimMat = load(DS{n}.stimfiles);   
        
        % Validate the MAT-file
         flds = fieldnames(stimMat);
         if (length(flds) == 1) 
             handles.stimfiles = getfield(stimMat, char(flds{1}));
             guidata(handles.figure1, handles)
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
         guidata(handles.figure1, handles)
    case {'.mat'}
         respMat = load(DS{n}.respfiles);   
           
         % Validate the MAT-file
         flds = fieldnames(respMat);
         if (length(flds) == 1)
             handles.respfiles = getfield(respMat, char(flds{1}));
             guidata(handles.figure1, handles)
         end
    otherwise
         errordlg(['Only ASCII and MAT binary fileformats work',...
                ' in this version.'], 'File Type Error', 'modal')
         return;
    end
       
    % Get the appropriate data for the index in selected
    set(handles.show_respfile_text,'String',[name ext]);

    % Now display input stimulus file
    global ampsamprate
    nb = size(handles.stimfiles, 1);
    nt = min(size(handles.respfiles,2),min(DS{n}.nlen,size(handles.stimfiles, 2)));
  
    xlabelrange = (1:nt)*ceil(1000/ampsamprate);
    
    global sDim 
    set(handles.spatialdomain, 'String', sDim);
   
    if strcmp(sDim, '1-D')
        % Make 2d frame button invisible
        set(handles.prev5frame, 'Visible', 'off')
        set(handles.next5frame, 'Visible', 'off')
        
        global initialFreq endFreq
        
        fstep = (endFreq - initialFreq)/nb;
        faxis = initialFreq:fstep:endFreq;
        axes(handles.stim_axes);
        imagesc(xlabelrange, faxis/1000, handles.stimfiles(:,1:nt));
        axis xy;
        %xlim([0,nt])
        ylabel('Frequency (kHz)')
        title('Stimulus')
        
    else
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
        
        stimstim = reshape(handles.stimfiles, splitX, splitX, nt);
        numList = 5;
        whole=[];
        
        for i=handles.frameIndex:handles.frameIndex+numList
            if i+1 > nt
                handles.frameIndex = 0;
                guidata(handles.figure1, handles)
                break;
            end
            whole = [whole stimstim(:,:,i+1)];
        end
        axes(handles.stim_axes);
        ttH = imagesc(whole);
        axis image
        set(get(ttH,'Parent'),'YTickLabel',[]);
        set(get(ttH,'Parent'),'YTick',[]);
        set(get(ttH,'Parent'),'XTickLabel',[]);
        set(get(ttH,'Parent'),'XTick',[]);
        
        % 5/12/03 - JXZ
        %colormap(redblue)
        colormap(gray)
        title([num2str(handles.frameIndex+1) '  to '...
                num2str(handles.frameIndex +6) ' Video Frames'])
    end
    
    
    ntrials = DS{n}.ntrials;
    if ntrials > 1
        numtrials = 10;
        if numtrials > ntrials
            numtrials = ntrials;
        end
        
        for ii = 1:numtrials
            subplot('position', [0.106 0.34585+(ii-1)*(0.280612/numtrials)...
               0.592014 0.280612/numtrials-0.001]); 
            plot(xlabelrange,handles.respfiles(ii,:), 'k');
            p_axis = axis;
            p_axis(2) = nt*ceil(1000/ampsamprate);
            axis(p_axis)
            %xlim([0, nt])
            axis off;
         
        end
    else
        axes(handles.spiketrain_axes)
        plot(xlabelrange,handles.respfiles);
        p_axis = axis;
        p_axis(2) = nt*ceil(1000/ampsamprate);
        axis(p_axis)
        
        ylabel('Spikes')
    end

    axes(handles.psth_axes);
    if ntrials > 1
        meanPSTH = mean(handles.respfiles);
    else
        meanPSTH = handles.respfiles;
    end
    
    %Show smoothed version of psth at window = 15
    wind1 = hanning(15)/sum(hanning(15));
    
    svagsm=conv(meanPSTH(1:nt),wind1);
    plot(xlabelrange,svagsm(8:length(svagsm)-7));
    p_axis = axis;
    p_axis(2) = nt*ceil(1000/ampsamprate);
    axis(p_axis)
    
    ylabel('Smoothed PSTH (spikes/s)')
    xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
    %xlim([0,nt])
    
    % Set the index pointer to 1 and save handles
    handles.Index = n;
    guidata(handles.figure1, handles)
    

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
    i = frameindex - 5;
    % If the index is less then one then set it equal to the index of the
    % last element in the Addresses array
    if i < 0
        i = 0;
    end
    % Update the index pointer to reflect the new index
    handles.frameIndex = i;
    guidata(h,handles)
    
    checknload(DS, handles.Index, handles);
    

% --------------------------------------------------------------------
function varargout = next5frame_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
    global DS;
    
    % Update the index pointer to reflect the new index
    handles.frameIndex = handles.frameIndex + 5;
    guidata(h,handles)
    
    checknload(DS, handles.Index, handles);
    
% --------------------------------------------------------------------
%  END of displayinput_GUI
% --------------------------------------------------------------------


% --------------------------------------------------------------------
function varargout = playstim_Callback(h, eventdata, handles, varargin)
    
    global rawDS
    % Get the index pointer and the files 
    index = handles.Index;
    
    [sig, fs] = wavread(rawDS{index}.stimfiles);
    soundsc(sig,fs);


% --------------------------------------------------------------------
function varargout = playspike_Callback(h, eventdata, handles, varargin)
    global DS ampsamprate
    % Get the index pointer and the files 
    index = handles.Index;
    tindex = handles.trialindex;
    
    spike = load(DS{index}.respfiles);
    %soundsc(spike(tindex,:), ampsamprate);
    soundsc(mean(spike), ampsamprate);


% --- Executes during object creation, after setting all properties.
function spiketrial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spiketrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function spiketrial_Callback(hObject, eventdata, handles)
    global DS
    newStr = get(hObject, 'String');
    
    NewVal = str2double(newStr);

    % Check that the entered value falls within the allowable range
    
    if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
        % Set default value for resp sample rate
        errordlg('Please enter valid trial number (positive number only).',...
            'Variable Error', 'modal')
        return; 
    elseif NewVal > DS{handles.Index}.ntrials
        NewVal = DS{handles.Index}.ntrials;
        set(h, 'String', NewVal);
    end

   handles.trialindex = NewVal;
   guidata(hObject, handles);
    
