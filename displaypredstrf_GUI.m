function varargout = displaypredstrf_GUI(varargin)
% DISPLAYPREDSTRF_GUI Application M-file for displaypredstrf_GUI.fig
%    FIG = DISPLAYPREDSTRF_GUI launch displaypredstrf_GUI GUI.
%    DISPLAYPREDSTRF_GUI('callback_name', ...) invoke the named callback.
%
%             STRFPAK: STRF Estimation Software
% Copyright �?003. The Regents of the University of California (Regents).
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
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
        hUIControls],'fontname', 'Times New Roman','units',...
        'normalized','fontunits','normalized');

    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);

    % Assign data index and control index
    handles.dataIndex = 1;
    set(handles.DataSet, 'string', num2str(handles.dataIndex));
    handles.Index = 1;
    handles.frameIndex = 1;
    handles.displaywpsth = 0;

    % 7/14/2005
    handles.stdindex = 1;
    handles.tolindex = 1;

    % Save handles to fig
    guidata(fig, handles);

    % Check if you have done prediction. If yes, go to result directory
    global Tol_val sDim predDS predinitialFreq predendFreq DS TimeLag
    global initialFreq endFreq ampsamprate Std_val preprocessOption


    set(handles.frameSet, 'String', 1);


    global outputPath
    if isempty(outputPath)|isempty(Tol_val)|isempty(sDim)|isempty(preprocessOption)
        % |isempty(predDS) %|isempty(predinitialFreq) |isempty(predendFreq)
        anws = questdlg('Have you done prediction of PSTH?',...
            'Calculation Question','Yes', 'No', 'Yes');
        switch anws
            case 'Yes'
                prompt={['Enter the path of prediction result files:']};
                defaultdir = fullfile(pwd, 'Output');
                def={defaultdir};
                dlgTitle='Directory for output';
                lineNo=1;
                % picture feature
                AddOpts.Resize='on';
                AddOpts.WindowStyle='normal';
                AddOpts.Interpreter='tex';
                datadir =inputdlg2(prompt,dlgTitle,lineNo,def,AddOpts);
                if isempty(datadir)
                    errordlg('You dont enter the path. Bye.')
                    return
                end
                outputdir = datadir{1};

                if not(exist(outputdir,'dir'))
                    errordlg('Directory not found... exiting.');
                    return
                end
                outputPath = outputdir;
                tempGV = fullfile(outputPath, 'GlobalVariables.mat');
                if not(exist(tempGV, 'file'))
                    errordlg('Wrong path for output files')
                    return
                end
                load(tempGV, 'Tol_val','Std_val','sDim','DS','TimeLag',...
                    'preprocessOption','initialFreq','endFreq','ampsamprate')

            case 'No'
                msgbox('You need do prediction first', 'Warning', 'modal');
                return
        end
    end
    if isempty(predDS)
        tempPV = fullfile(outputPath, 'predVariables.mat');
        if ~exist(tempPV, 'file')
            errordlg('There is no prediction result for this input set.')
            return
        end
        load(fullfile(outputPath, 'predVariables.mat'),'predDS',...
            'predinitialFreq','predendFreq')
    end

    pfilename = fullfile(outputPath, sprintf('predResult_EstSpike_Tol%d.mat',...
        handles.Index));
    if not(exist(pfilename))
        errordlg('No validation yet, please do validation first','!!No Validation!!','modal');
        return;
    end
    handles.estSpike = Check_And_Load(pfilename);

    pfilename = fullfile(outputPath, 'predResult_avgSpike1.mat');
    handles.avgSpike1 = Check_And_Load(pfilename);

    pfilename = fullfile(outputPath, 'predResult_avgSpike2.mat');
    handles.avgSpike2 = Check_And_Load(pfilename);

    spike_psth = Check_And_Load(fullfile(outputPath,'spike_psth_count.mat'));
    ntrials_proper=spike_psth(end);
    spike_psth=spike_psth(1:end-1);
    %ntrials_index=find(spike_psth>=ntrials_proper);
    ntrials_index = 1:length(predDS);
    handles.ntrialsIndex = ntrials_index;
    clear spike_psth;

    predlen = length(handles.ntrialsIndex);
    if predlen ==1
        set(handles.datasetSlider, 'min', 0, 'max', predlen, 'value',0,'sliderstep', [1 1]);
    elseif predlen > 1
        set(handles.datasetSlider, 'Min', 1, 'Max', predlen,'value',1,...
            'sliderstep', [1/(predlen-1) 1/(predlen-1)]);
    else
        warndlg('No input selected yet. Exiting...')
        return;
    end
    framesize = length(handles.estSpike{handles.dataIndex}{handles.stdindex});
    numList = 6;
    if framesize ==1
        set(handles.frameSetSlider, 'min', 0, 'max', 1, 'value',0,'sliderstep', [1 1]);
    else
        set(handles.frameSetSlider, 'Min', 1, 'Max', framesize,'value',1,...
            'sliderstep', [numList/(framesize-1) numList/(framesize-1)]);
    end

    tollen = length(Tol_val);
    stdlen = length(Std_val);

    if tollen ==1
        set(handles.tolvalSlider, 'min', 0, 'max', tollen, 'value',0,'sliderstep', [1 1]);
    elseif tollen > 1
        set(handles.tolvalSlider, 'Min', 1, 'Max', tollen,'value',1,...
            'sliderstep', [1/(tollen-1) 1/(tollen-1)]);
    end
    if stdlen ==1
        set(handles.stdvalSlider, 'min', 0, 'max', stdlen, 'value',0,'sliderstep', [1 1]);
    elseif stdlen > 1
        set(handles.stdvalSlider, 'Min', 1, 'Max', stdlen,'value',1,...
            'sliderstep', [1/(stdlen-1) 1/(stdlen-1)]);
    end
    % Save handles to fig
    guidata(fig, handles);

    set(handles.tolvalshow, 'visible','on');
    set(handles.stdvalshow, 'visible', 'on');
    set(handles.tolvalSlider, 'visible', 'on');
    set(handles.stdvalSlider, 'visible', 'on');
    set(handles.tolvalshow, 'string', Tol_val(handles.tolindex));
    set(handles.tolvalSlider,'value', handles.tolindex);
    set(handles.stdvalshow, 'string', Std_val(handles.stdindex));
    set(handles.stdvalSlider,'value', handles.stdindex);

    % show the current path
    set(handles.outdirshow, 'string', outputPath);
    set(handles.outdirshow, 'visible', 'on');

    %     handles.tol_val = Tol_val;
    %     handles.std_val = Std_val;
    %     set(handles.predds,'struct', predDS);
    %
    %     % Initially display tol val and pred files
    %     tol = handles.tol_val(handles.Index);
    %     set(handles.tolval_show, 'String', tol);
    %     stds = handles.std_val(handles.stdIndex);
    %     set(handles.stdval_show, 'String', stds);

    [path, name, ext] = fileparts(predDS{1}.stimfiles);
    set(handles.stimfile_show, 'String',[name ext]);
    [path, rame, rxt] = fileparts(predDS{1}.respfiles);
    set(handles.respfile_show, 'String',[rame rxt]);

    global smoothwindow
    newStr_smoothwindow = get(handles.smoothwindow, 'String');
    smoothwindow = str2double(newStr_smoothwindow);
    display_predstrf(handles);

    % save handles new data
    guidata(fig, handles);

    if nargout > 0
        varargout{1} = fig;
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    %try
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    %catch
    %    disp(lasterr);
    %end

end

% ====================================================================
function varargout = prev_button_Callback(h, eventdata, handles, varargin)
% ====================================================================
% Get the index pointer and the files
index = handles.Index;

% update index
i = index - 1;
% If the index is less then one then set it equal to the index of the
% last element in the Addresses array
if i < 1
    i = length(handles.tol_val);
end
handles.Index = i;


tol = handles.tol_val(handles.Index);
set(handles.tolval_show, 'String', tol);

% Assign the predict result to handles
global outputPath
pfilename = fullfile(outputPath, sprintf('predResult_EstSpike_Tol%d.mat',...
    handles.Index));
handles.estSpike = Check_And_Load(pfilename);

%     pfilename = fullfile(outputPath, sprintf('predResult_avgSpike1_Tol%d.mat',...
%                  handles.Index));
%     handles.avgSpike1 = Check_And_Load(pfilename);
%
%     pfilename = fullfile(outputPath, sprintf('predResult_avgSpike2_Tol%d.mat',...
%                  handles.Index));
%     handles.avgSpike2 = Check_And_Load(pfilename);
%
%     pfilename = fullfile(outputPath, sprintf('predResult_meanrate_Tol%d.mat',...
%                  handles.Index));
%     handles.meanrate = Check_And_Load(pfilename);

guidata(h, handles);
display_predstrf(handles);


% ====================================================================
function varargout = next_button_Callback(h, eventdata, handles, varargin)
% ====================================================================
% Get the index pointer and the files
index = handles.Index;

% update index
i = index + 1;
% If the index is less then one then set it equal to the index of the
% last element in the Addresses array
if i > length(handles.tol_val)
    i = 1;
end
handles.Index = i;

tol = handles.tol_val(handles.Index);
set(handles.tolval_show, 'String', tol);


% Assign the predict result to handles
global outputPath
pfilename = fullfile(outputPath, sprintf('predResult_EstSpike_Tol%d.mat',...
    handles.Index));
handles.estSpike = Check_And_Load(pfilename);

guidata(h, handles);
display_predstrf(handles);


% ====================================================================
function varargout = close_button_Callback(h, eventdata, handles, varargin)
% ====================================================================
delete(handles.figure1);



% ====================================================================
%  core function displaypredstrf to plot the graph
% ====================================================================
function display_predstrf(handles)
global sDim predendFreq predinitialFreq ampsamprate smoothwindow predDS
global preprocessOption

% get data index and predResult
trialindex = handles.dataIndex;
index = handles.ntrialsIndex(trialindex);
%tolindex = handles.Index;
tolindex=1;
estSpike = handles.estSpike;
%estSpike = estSpike_temp{handles.stdIndex};
avgSpike1 = handles.avgSpike1;
avgSpike2 = handles.avgSpike2;

%save('debugging.mat','estSpike','avgSpike1','avgSpike2');
handles.stimfile = Check_And_Load(predDS{index}.stimfiles);
handles.respfile = Check_And_Load(predDS{index}.respfiles);
[nb, nt] = size(handles.stimfile);


xlabelrange = ceil((1:nt)*1000/ampsamprate);
% Display predicted STRF

if strcmp(sDim, '1-D')

    % Display stim 1-D
    taxis =ceil((1:nt)*1000/ampsamprate);
    if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
        flabel = logspace(log10(predinitialFreq), log10(predendFreq), nb);
        axes(handles.stim_axes);

        pcolor(taxis, flabel, handles.stimfile); shading interp;

    else

        fstep = (predendFreq - predinitialFreq)/nb;
        faxis = predinitialFreq:fstep:predendFreq;
        axes(handles.stim_axes);
        imagesc(taxis, faxis, handles.stimfile);
        axis xy;
    end
    s_axis = axis;
    s_axis_y = s_axis(2);
    title('Stimulus file for prediction');
    ylabel('Frequency (Hz)')
    %xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
    global stimsamprate sDim
    if ~isempty(stimsamprate)
        if sDim == '2-D'
            timebin = round(1000/stimsamprate);  % in ms
            xlabel(['Time (in ', num2str(timebin), ' ms)']);

        else

            xlabel('Time (in ms)');
        end
    else
        xlabel('Frames');
    end

    %xlim([0,nt])

elseif strcmp(sDim, '2-D')
    % 2-D display
    global preprocessOption
    if strcmp(preprocessOption, 'Fourier Power transform')
        phasecount=1;
        tbincount = nt;
        kcount = 1;
        spacebincount = nb;
        chancount=spacebincount/phasecount;
        Xmax=sqrt(chancount*2);
        [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);

        tsf=zeros(Xmax*Xmax,tbincount,kcount);
        tsf(cfilt,:,:)=squeeze(sum(reshape(handles.stimfile(:,1:nt),chancount,phasecount,...
            tbincount,kcount),2));

        tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
        stimstim=reshape(tsf,Xmax,Xmax,tbincount);
        splitX = Xmax;
    else
        nt = min(size(handles.respfile,2), nt);
        splitX = floor(sqrt(nb));
        stimstim = reshape(handles.stimfile(:,1:nt), splitX, splitX, nt);
    end
    numList = 10;
    whole=[];
    frameRange = handles.frameIndex:handles.frameIndex+numList;

    if length(find(frameRange>nt)) == 0
        whole = reshape(stimstim(:,:,frameRange), splitX, splitX*length(frameRange));
    else
        whole = reshape(stimstim(:,:,frameRange(find(frameRange<=nt))), splitX,...
            splitX*length(frameRange(find(frameRange<=nt))));
    end

    axes(handles.stim_axes);
    ttH = imagesc(whole);

    colormap(redblue)
    axis image
    set(get(ttH,'Parent'),'YTickLabel',[]);
    set(get(ttH,'Parent'),'YTick',[]);
    set(get(ttH,'Parent'),'XTickLabel',[]);
    set(get(ttH,'Parent'),'XTick',[]);
    title([num2str(handles.frameIndex) '  to '...
        num2str(handles.frameIndex +numList) ' Video Frames'])
    %s_axis_y = 6*splitX;
elseif strcmp(sDim, '0-D')
    % Display stim 1-D
    taxis =ceil((1:nt)*1000/ampsamprate);

    axes(handles.stim_axes);
    plot(taxis, handles.stimfile);
    axis tight;

    title('Stimulus file for prediction');


end

% Display avg of the part of neuron response
axes(handles.mixstrf_axes);
%tt = mean([avgSpike1{index}  avgSpike2{index}], 2);

% Smooth two half psth by window size = smoothwindow
if smoothwindow == 0
    smoothwindow = 1;
end
nshift = floor(smoothwindow/2);
wind1 = hanning(2*nshift +1)/sum(hanning(nshift*2+1));
rawPSTH = (avgSpike1{trialindex} + avgSpike2{trialindex})/2;
svagsm=conv(rawPSTH, wind1);
len_psth = length(rawPSTH);
smoothpsth = svagsm(nshift+1:end-nshift);

% JXZ: 9/13/2005
%  Adding displaying portion of psth corresponding to 2-D stim display
if strcmp(sDim, '1-D') | strcmp(sDim, '0-D')

    plot(ceil((1:len_psth)*1000/ampsamprate),...
        smoothpsth*ampsamprate, 'b');
    hold on;
elseif strcmp(sDim, '2-D')
    set(handles.displaywhole, 'Visible', 'On');
    if handles.displaywpsth == 1
        labelrange = 1:nt;
    else
        labelrange = frameRange(find(frameRange<=nt));
    end
    smoothpsth = svagsm(nshift+1:end-nshift);
    plot(xlabelrange(labelrange),smoothpsth(labelrange)*ampsamprate, 'b');
    hold on;
end

% Now display predicted psth
%svagsm=conv(estSpike{index}+handles.meanrate,wind1);
global  allow_negative_rates 

% Junli 6/30/2005
%  Check which psth-option is used.
global timevary_PSTH outputPath predDS
[predstim_avg, avgpsth] = cal_AVG(predDS, nb);
if timevary_PSTH == 0
    %loadpsth = load(fullfile(outputPath, 'stim_avg.mat'), 'constmeanrate');
    meanrate = avgpsth;
    model_prediction = (estSpike{trialindex}{handles.stdindex}+meanrate)*ampsamprate;
    if ~allow_negative_rates
        model_prediction(find(model_prediction<0))=0.0;
    end

    curtitle = 'Prediction model: Mean rate + STRF prediction';
else
    %loadpsth = load(fullfile(outputPath, 'stim_avg.mat'), 'Avg_psth');
    mean_rate = avgpsth;
    curtitle = 'Prediction model: Time varying rate + STRF prediction';
    model_prediction = (estSpike{trialindex}{handles.stdindex}+mean_rate(index,1:len_psth)')*ampsamprate;
    if isempty(allow_negative_rates)
        allow_negative_rates = 0;
    end
    if ~allow_negative_rates
        model_prediction(find(model_prediction<0))=0.0;
    end

end

if strcmp(sDim, '1-D')
    plot(ceil((1:len_psth)*1000/ampsamprate), model_prediction, 'r');
    p_axis = axis;
    p_axis(2) = s_axis_y;
    axis(p_axis);
elseif strcmp(sDim, '2-D')   % Display for 2-D option
    %labelrange = frameRange(find(frameRange<=nt));
    plot(xlabelrange(labelrange), model_prediction(labelrange), 'r');
    p_axis = axis;
    p_axis(1) = xlabelrange(labelrange(1));
    p_axis(2) = xlabelrange(labelrange(end));
    axis(p_axis);
elseif strcmp(sDim, '0-D')
    plot(ceil((1:len_psth)*1000/ampsamprate), model_prediction, 'r');
    axis tight;
end


legend('Response', 'Prediction');
title(curtitle);

%xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
global stimsamprate sDim
if ~isempty(stimsamprate)
    if sDim == '2-D'
        timebin = round(1000/stimsamprate);  % in ms
        xlabel(['Time (in ', num2str(timebin), ' ms)']);

    else

        xlabel('Time (in ms)');
    end
else
    xlabel('Frames');
end

ylabel('Spike Rate (spikes/s)');
hold off;

% Calculate raw r between raw psth and predicted psth
raw_r = diag(corrcoef(smoothpsth*ampsamprate,model_prediction), 1);
set(handles.pred_r, 'String', raw_r);


% ====================================================================
function varargout = back_datafile_Callback(h, eventdata, handles, varargin)
% ====================================================================
% Get the index pointer and the files
index = handles.dataIndex;
ntrials_Index = handles.ntrialsIndex;

% update index
i = index - 1;
% If the index is less then one then set it equal to the index of the
% last element in the Addresses array
if i < 1
    i = 1;
end
handles.dataIndex = i;
global predDS
[path, name, ext] = fileparts(predDS{ntrials_Index(i)}.stimfiles);
set(handles.stimfile_show, 'String',[name ext]);
[path, rame, rxt] = fileparts(predDS{ntrials_Index(i)}.respfiles);
set(handles.respfile_show, 'String',[rame rxt]);
% Also reset frameIndex
handles.frameIndex = 0;
guidata(h, handles);

display_predstrf(handles);


% ====================================================================
function varargout = forward_datafile_Callback(h, eventdata, handles, varargin)
% ====================================================================
% Get the index pointer and the files
index = handles.dataIndex;
ntrials_Index = handles.ntrialsIndex;

% update index
i = index + 1;
% If the index is less then one then set it equal to the index of the
% last element in the Addresses array
if i > length(ntrials_Index)
    i = length(ntrials_Index);
end
handles.dataIndex = i;
global predDS
[path, name, ext] = fileparts(predDS{ntrials_Index(i)}.stimfiles);
set(handles.stimfile_show, 'String',[name ext]);
[path, rame, rxt] = fileparts(predDS{ntrials_Index(i)}.respfiles);
set(handles.respfile_show, 'String',[rame rxt]);

% Also reset frameIndex
handles.frameIndex = 1;
guidata(h, handles);

display_predstrf(handles);


% --------------------------------------------------------------------
function varargout = prev5frames_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
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
display_predstrf(handles);

% --------------------------------------------------------------------
function varargout = next5frames_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Update the index pointer to reflect the new index
global ampsamprate
nt = length(handles.estSpike{handles.dataIndex}{handles.stdIndex});

% Update the index pointer to reflect the new index
range = handles.frameIndex+10;
if range < nt
    handles.frameIndex = range;
else
    handles.frameIndex = 1;
end

% Reset displaywhole option if next5frame has been clicked.
handles.displaywpsth = 0;

guidata(h,handles)
display_predstrf(handles);


% --------------------------------------------------------------------
function varargout = smoothwindow_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
newStr_smoothwindow = get(h, 'String');
NewVal = str2double(newStr_smoothwindow);

% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 0) | isnan(NewVal)
    % Set default value for resp sample rate
    errordlg('Please enter valid sampling rate (positive number only).',...
        'Variable Error', 'modal')
    return;
end


% Assign global variable to ampsamprate
global smoothwindow;
smoothwindow = NewVal;

display_predstrf(handles);


% --- Executes on button press in prestdval.
function prestdval_Callback(hObject, eventdata, handles)
% hObject    handle to prestdval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the index pointer and the files
index = handles.stdIndex;

% update index
i = index - 1;
% If the index is less then one then set it equal to the index of the
% last element in the Addresses array
if i < 1
    i = length(handles.std_val);
end
handles.stdIndex = i;


stds = handles.std_val(handles.stdIndex);
set(handles.stdval_show, 'String', stds);

guidata(hObject, handles);
display_predstrf(handles);

% --- Executes on button press in nextstdval.
function nextstdval_Callback(hObject, eventdata, handles)
% hObject    handle to nextstdval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the index pointer and the files
index = handles.stdIndex;

% update index
i = index + 1;
% If the index is less then one then set it equal to the index of the
% last element in the Addresses array
if i > length(handles.std_val)
    i = 1;
end
handles.stdIndex = i;

stds = handles.std_val(handles.stdIndex);
set(handles.stdval_show, 'String', stds);

guidata(hObject, handles);
display_predstrf(handles);





% --- Executes on button press in displaywhole.
function displaywhole_Callback(hObject, eventdata, handles)
handles.displaywpsth = 1;
guidata(hObject, handles);
display_predstrf(handles);


% --- Executes on button press in smoothpsth.
function smoothpsth_Callback(hObject, eventdata, handles)
helpdlg({'Smoothing PSTH: ',...
    ' ',...
    ' This flag is used for smoothing psth when displaying psth.',...
    ' It is in ms. You can type new number in the text box to modify it.',...
    ' '},...
    'Smooth PSTH Help');

% ====================================================================
function varargout = help_button_Callback(h, eventdata, handles, varargin)
% ====================================================================
ttlStr = get(handles.figure1, 'Name');
hlpStr = [...
    '                                                                            '
    '    STRFPAK: Display Predicted Results Window                               '
    '                                                                            '
    ' The "Display predicted reults" window show prediction and input data on one'
    ' figure. The left panel displays stimulus plots on the top and predcition   '
    ' with neuron response (PSTH) on the bottom. The right panel shows all the   '
    ' related information about the figures on the left.                         '
    '                                                                            '
    '  On Left Panel:                                                            '
    '         1. The plot of stimulus file for prediction: If it is 1-D, this    '
    '            is spectrogram. If it is 2-D, it is a series of frames.         '
    '         2. The plot of predicted neuron response with actual neuron        '
    '            response. The title of the plot is based on how you set Model   '
    '            parameter from "Calculation Parameter".                         '
    '         3. Smooth_PSTH field provide the window size for smoothing neuron  '
    '            response. You can type different values in the text field and   '
    '            see the difference of response line.                            '
    '         4. Display whole psth: If input data is 2-D in spatio domain, this '
    '            button will show up. If it is clicked, the plot of prediction   '
    '            with response will show the whole time duration instead of just '
    '            10 frames long.                                                 '
    '                                                                            '
    '  On Right Panel:                                                           '
    '         Data Set: This field is for more than one predicted files chosen.  '
    '         Data Set Slider: The user can use mouse click to display the next  '
    '            results from the next predicted files if available.             '
    '         Pred Stim filename: The text box shows the file name.              '
    '         Pred Resp filesname: The text box show sthe file name.             '
    '         Frame Set: If the data are 2-D data, this option shows the next    '
    '             frame set if more than 6 frames.                               '
    '         Tol Val box: shows the tolerance values used for the STRF.         '
    '         Tol val slider: use mouse can choose next or prev tol val          '
    '                  and the corresponding STRF.                               '
    '         Spareness box: shows the sparseness values used for the            '
    '                  STRF.                                                     '
    '         Spareness slider: use mouse can choose next or prev                '
    '               sparseness value and the corresponding STRF.                 '
    '         Predicted r: shows the corr coeff between the prediction and       '
    '              neural response for the stimulus shown above.                 '
    '         Help: This window will pop up. For detailed doc, please refer      '
    '              to the user manual from http://strfpak.berkeley.edu.          '
    '         Close: close all the figures.                                      '
    '                                                                            '
    ' Updated by Junli, May 2006.                                                '
    '                                                                            '];
myFig = handles.figure1;
helpwin(hlpStr, ttlStr);



% --- Executes on slider movement.
function tolvalSlider_Callback(hObject, eventdata, handles)
global Tol_val
newSlider = round(get(hObject, 'Value'));
if newSlider <= 0
    set(handles.tolvalshow, 'String', Tol_val(handles.tolindex));
else
    newSlider = round(newSlider);
    set(handles.tolvalshow, 'String', Tol_val(newSlider));
    handles.tolindex = newSlider;
end
global outputPath
pfilename = fullfile(outputPath, sprintf('predResult_EstSpike_Tol%d.mat',...
    handles.tolindex));
handles.estSpike = Check_And_Load(pfilename);

guidata(handles.figure1, handles);
display_predstrf(handles);


% --- Executes during object creation, after setting all properties.
function tolvalSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolvalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function tolvalshow_Callback(hObject, eventdata, handles)
global Tol_val
maxlength = length(Tol_val);
newIndex = str2double(get(hObject,'String'));
if newIndex <=0 || newIndex > maxlength
    errordlg(['Tol_val_index is out of range. It need to be between 1 and ', num2str(maxlength)]);
    return;
else
    newIndex = round(newIndex);
    handles.tolindex = newIndex;
    %set(hObject, 'String', num2str(handles.Index));
    % Update slider value based on dataSet value
    set(handles.tolvalSlider, 'String', newIndex);
end
global outputPath
pfilename = fullfile(outputPath, sprintf('predResult_EstSpike_Tol%d.mat',...
    handles.tolindex));
handles.estSpike = Check_And_Load(pfilename);
guidata(handles.figure1, handles);
display_predstrf(handles);


% --- Executes during object creation, after setting all properties.
function tolvalshow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolvalshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function stdvalSlider_Callback(hObject, eventdata, handles)
global Std_val
newSlider = round(get(hObject, 'Value'));
if newSlider <= 0
    set(handles.stdvalshow, 'String', Std_val(handles.stdindex));
else
    newSlider = round(newSlider);
    set(handles.stdvalshow, 'String', Std_val(newSlider));
    handles.stdindex = newSlider;
end
guidata(handles.figure1, handles);
display_predstrf(handles);

% --- Executes during object creation, after setting all properties.
function stdvalSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stdvalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function stdvalshow_Callback(hObject, eventdata, handles)
global Std_val
maxlength = length(Std_val);
newIndex = str2double(get(hObject,'String'));
if newIndex <=0 || newIndex > maxlength
    errordlg(['Std_val_index is out of range. It need to be between 1 and ', num2str(maxlength)]);
    return;
else
    newIndex = round(newIndex);
    handles.stdindex = newIndex;
    %set(hObject, 'String', num2str(handles.Index));
    % Update slider value based on dataSet value
    set(handles.stdvalSlider, 'String', newIndex);
end
guidata(handles.figure1, handles);
display_predstrf(handles);


% --- Executes during object creation, after setting all properties.
function stdvalshow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stdvalshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outdirshow_Callback(hObject, eventdata, handles)
global outputPath
set(handles.outdirshow, 'String', outputPath);


% --- Executes during object creation, after setting all properties.
function outdirshow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function datasetSlider_Callback(hObject, eventdata, handles)
% global predDS
% maxlength = length(predDS);
maxlength = length(handles.ntrialsIndex);
newSlider = round(get(hObject, 'Value'));
if newSlider <= 0 | newSlider > maxlength | round(newSlider) == 0

    set(handles.DataSet, 'String', handles.dataIndex);
else
    set(handles.DataSet, 'String', round(newSlider));
    handles.dataIndex = round(newSlider);
end
guidata(handles.figure1, handles);
global predDS
[path, name, ext] = fileparts(predDS{handles.dataIndex}.stimfiles);
set(handles.stimfile_show, 'String',[name ext]);
[path, rame, rxt] = fileparts(predDS{handles.dataIndex}.respfiles);
set(handles.respfile_show, 'String',[rame rxt]);
% Also reset frameIndex
handles.frameIndex = 1;
display_predstrf(handles);


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
global predDS
maxlength = length(predDS);
newIndex = str2double(get(hObject,'String'));
if newIndex <=0 || newIndex > maxlength
    warndlg(['DataSet is out of range. It need to be between 1 and ', num2str(maxlength)]);
    set(hObject, 'String', num2str(handles.Index));
else
    handles.Index = newIndex;
    %set(hObject, 'String', num2str(handles.Index));
    % Update slider value based on dataSet value
    set(handles.datasetSlider, 'String', newIndex);
end
guidata(handles.figure1, handles);
[path, name, ext] = fileparts(predDS{handles.dataIndex}.stimfiles);
set(handles.stimfile_show, 'String',[name ext]);
[path, rame, rxt] = fileparts(predDS{handles.dataIndex}.respfiles);
set(handles.respfile_show, 'String',[rame rxt]);
% Also reset frameIndex
handles.frameIndex = 1;
display_predstrf(handles);

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

function frameSet_Callback(hObject, eventdata, handles)
maxlength = length(handles.estSpike{handles.dataIndex}{handles.stdindex});
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
guidata(handles.figure1, handles);
handles.displaywpsth = 0;
display_predstrf(handles);


% --- Executes during object creation, after setting all properties.
function frameSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function frameSetSlider_Callback(hObject, eventdata, handles)
maxlength = length(handles.estSpike{handles.dataIndex}{handles.stdindex});
newSlider = round(get(hObject, 'Value'));
if newSlider <= 0
    set(handles.frameSet, 'String', 1);
elseif newSlider >= maxlength
    set(handles.frameSet', 'String', maxlength);
else
    set(handles.frameSet, 'String', round(newSlider));
    handles.frameIndex = round(newSlider);
end
handles.displaywpsth = 0;
guidata(handles.figure1, handles);
display_predstrf(handles);



% --- Executes during object creation, after setting all properties.
function frameSetSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameSetSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% ====================================================================
%  END of DISPLAYPREDSTRF_GUI.m
% ====================================================================


% --- Executes on button press in Show_STRF.
function Show_STRF_Callback(hObject, eventdata, handles)
% hObject    handle to Show_STRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tolmax = handles.tolindex;
stdmax = handles.stdindex;
global outputPath DS Tol_val Std_val sDim preprocessOption predendFreq predinitialFreq ampsamprate
% then load predicted STRF and display the best STRF
bestStrf = load(fullfile(outputPath,...
    ['strfResult_Tol',num2str(tolmax),'.mat']));

% find the best filter
% If only one data set is chosen for STRF estimation, we have to choose
% STRF instead of STRFJN and could not do smoothing on the STRF since STRFJNstd is zero.
if length(DS) == 1
    forward = bestStrf.STRF_Cell;
    [nb, nt] = size(forward);
    handles.nJN = 1;
else
    forwardJN = bestStrf.STRFJN_Cell;
    forwardJNstd = bestStrf.STRFJNstd_Cell;
    [nb, nt,nJN] = size(forwardJN);

    forward = fast_filter_filter(bestStrf.STRF_Cell,mean(bestStrf.STRFJNstd_Cell,3),Std_val(stdmax));
    %
    %     forwardJN_s = zeros(nb,nt,nJN);
    %     for iJN =1:nJN
    %         %  find the filtered JN STRF
    %         forwardJN_s(:,:,iJN) = fast_filter_filter(mean(forwardJN(:,:,[1:(iJN-1) (iJN+1):nJN]),3), forwardJNstd(:,:,iJN), Std_val(stdmax));
    %     end
    %
    %     % Take mean of Jackknifed STRFs
    %     forward = squeeze(mean(forwardJN_s,3));
end

% Now display it
maxforward = max(max(forward));
minforward = min(min(forward));
absforward = max(abs(minforward),abs(maxforward));

t = -(nt-1)/2:(nt-1)/2;
f = 1:nb;
figure;
if strcmp(sDim, '1-D')

    if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
        flabel = logspace(log10(predinitialFreq), log10(predendFreq), nb);
        pcolor(t*ceil(1000/ampsamprate), flabel/1000, 3*forward); shading interp;
        caxis([-absforward absforward]);
    else
        fstep = (predendFreq - predinitialFreq)/nb;
        faxis = predinitialFreq:fstep:predendFreq;
        imagesc(t*ceil(1000/ampsamprate),faxis/1000,forward,[-absforward absforward]);
        axis xy;
    end
    xlabel('Time (ms)');
    ylabel('Frequency (kHz)');
    title(['STRF Tol ' num2str(Tol_val(tolmax)) ' Std ' num2str(Std_val(stdmax)) '.']);
else
    disp(['Sorry, this option doesn''t exist yet.'  char(10) ...
        'Please email the STRFPAK team if you would like it.']);
    %     set(handles.resp_axes, 'visible', 'off');
    %     % 2-D display
    %     splitX = floor(sqrt(nb));
    %     STA = reshape(forward, splitX, splitX, nt);
    %     amax = max(abs(STA(:)));
    %     amin = -amax;
    %
    %     numList = handles.numList;
    %     if numList > nt
    %         numList = nt;
    %     end
    %     global ampsamprate
    %     binsize = ceil((1/ampsamprate)*1000);
    %     titlestr = 'Latency(ms):';
    %
    %     % display total 30 images
    %     for fr = handles.imageIndex*numList+1:numList+handles.imageIndex*numList-(nt-1)/2
    %         if fr > nt
    %             break
    %         end
    %         % display them in appropriate axes
    %         axesStr = strcat('subaxes', num2str(fr-handles.imageIndex * numList-1));
    %         axes(getfield(handles, axesStr));
    %
    %         if amin~=amax,
    %             ttH = imagesc(STA(:,:,fr+(nt-1)/2),[amin,amax]);
    %
    %         else
    %             ttH = imagesc(zeros(size(STA,1),size(STA,2)));
    %         end
    %
    %         if size(STA,1)>1 & size(STA,2)>1,
    %             axis image
    %             set(get(ttH,'Parent'),'YTickLabel',[]);
    %             set(get(ttH,'Parent'),'XTickLabel',[]);
    %             %axis off
    %         end
    %
    %         if (fr-handles.imageIndex * numList -1 ) ~= 0
    %             curTitle = num2str((fr-1)*binsize);
    %         else
    %             curTitle = 'Latency(ms)';
    %         end
    %         title(curTitle);
    %         colormap(redblue);
    %     end
end

