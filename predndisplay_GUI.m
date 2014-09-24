function varargout = predndisplay_GUI(varargin)
%
% PREDNDISPLAY_GUI Application M-file for predndisplay_GUI.fig
%    FIG = PREDNDISPLAY_GUI launch predndisplay_GUI GUI.
%    PREDNDISPLAY_GUI('callback_name', ...) invoke the named callback.
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
%
% Modified by JXZ, 7/14/2005
%   Change plot for different std values


if nargin == 0  % LAUNCH GUI

    fig = openfig(mfilename, 'reuse');
    %set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

    % for resize property
    set(fig, 'resize', 'on');
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
           hUIControls],'units','normalized','fontunits','normalized', 'fontname', 'times');

    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    guidata(fig, handles);
   
    % Check whether if we have done validation
    global Tol_val sDim TimeLag DS Std_val smoothVect
    global StimDS predinitialFreq predendFreq
    global outputPath 
    if isempty(outputPath) | isempty(Tol_val) | isempty(Std_val)
         anws = questdlg('Have you done validation of STRF?',...
                'Calculation Question','Yes', 'No', 'Yes');
         switch anws
         case 'Yes'
                prompt={['Enter the path of validation result files:']};
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
               
               % Rest part for tol_values
               tempGV = fullfile(outputPath, 'Global_Variables.mat');
               if not(exist(tempGV, 'file'))
                    errordlg('Wrong path for output files')
                    return
               end
               load(tempGV, 'Tol_val', 'sDim', 'Std_val')

%                % Load prediction and validate data set and their related parameters
%                load(fullfile(outputPath, 'predVariables.mat'),'predDS',...
%                        'predinitialFreq','predendFreq')

           otherwise
               msgbox('You need do validation first', 'Warning', 'modal');
               return
           end
    end
   
    handles.fidx = 1;
    set(handles.DataSet, 'string', num2str(handles.fidx));
    
    predlen = length(StimDS);
    if predlen ==1
        set(handles.datasetSlider, 'min', 0, 'max', predlen, 'value',0,'sliderstep', [1 1]);
    elseif predlen > 1
        set(handles.datasetSlider, 'Min', 1, 'Max', predlen,'value',1,...
            'sliderstep', [1/(predlen-1) 1/(predlen-1)]);
    else
        warndlg('No input selected yet. Exiting...')
        return;
    end
    
    [path, name, ext] = fileparts(StimDS{handles.fidx}.stimfiles); 
    set(handles.stimfile_show, 'String',[name ext]);
   
    global novalidation
    if novalidation ==1
        % show tolval and stdval tags
        set(handles.tolvaltext, 'visible', 'on');
        set(handles.tolvalshow, 'visible', 'on');
        set(handles.tolvalSlider, 'visible', 'on');
        set(handles.stdvaltext, 'visible', 'on');
        set(handles.stdvalshow, 'visible', 'on');
        set(handles.stdvalSlider, 'visible', 'on');

        % Assign tol_val and tol_val index to handles
        handles.stdIndex = 1;
        handles.Index = 1;
        handles.optionflg = 1;
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

        % display tol_val on the screen
        set(handles.tolvalshow, 'String', Tol_val(handles.Index));

        % display tol_val on the screen
        set(handles.stdvalshow, 'String', Std_val(handles.stdIndex));
        
        get_strf(handles);
    else
        % 
        answ = questdlg('Do you want to predict using best strf only?','Warning Message',...
            'Use best strf', 'Use all strfs', 'Use best strf');
        switch answ
            case 'Use best strf'
                get_beststrf(handles);
            case 'Use all strfs'
                % show tolval and stdval tags
                set(handles.tolvaltext, 'visible', 'on');
                set(handles.tolvalshow, 'visible', 'on');
                set(handles.tolvalSlider, 'visible', 'on');
                set(handles.stdvaltext, 'visible', 'on');
                set(handles.stdvalshow, 'visible', 'on');
                set(handles.stdvalSlider, 'visible', 'on');

                % Assign tol_val and tol_val index to handles
                handles.stdIndex = 1;
                handles.Index = 1;
                handles.optionflg = 1;
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

                % display tol_val on the screen
                set(handles.tolvalshow, 'String', Tol_val(handles.Index));

                % display tol_val on the screen
                set(handles.stdvalshow, 'String', Std_val(handles.stdIndex));

                get_strf(handles);
        end
    end
    % save figure's content
    %guidata(fig, handles);

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
function varargout = close_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
delete(handles.figure1);

% --------------------------------------------------------------------
function varargout = help_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
ttlStr = get(handles.figure1, 'Name');
hlpStr = [...
    '                                                                            '
    '    STRFPAK: PREDICT & DISPLAY RESULTS                                      '
    '                                                                            '
    '  The window is also called the Neural Prediction Challenge window. Since it'
    ' predicts on the new stimulus and displays prediction results only. If the  '
    ' validation has been done early, it can predict using the best strf only or '
    ' the different strfs chosen by smoothness and sparseness parameters.        '
    '                                                                            '
    ' The top right panel display the strf used for this prediction.             '
    ' The top left figure shows the new stimulus and the bottom left figure      '
    ' displays the neuron prediction. The bottom right panel shows information.  '
    '                                                                            '
    ' Updated by Junli, May. 2006.                                               '
    '                                                                            '];
myFig = handles.figure1;
helpwin(hlpStr, ttlStr);



% --- Executes on button press in compsave.
function compsave_Callback(hObject, eventdata, handles)
% calculate  avg of new  stim
global StimDS
nband = handles.Nband;
ndata_files = length(StimDS);
stim_avg = zeros(nband, 1);
count_avg = 0;
lin_flag = 1;
sil_window = 0;
end_window = 0;

tempWait = waitbar(0,...
        'Predicting is in progress, please wait ...');
% ========================================================
% calculate the output over all the data file 
% ========================================================
for n = 1:ndata_files
     
    % load stimulus files
    stim_env = Check_And_Load(StimDS{n}.stimfiles);
    
    nt = size(stim_env, 2);
    % take logrithm of data based on lin_flag 
    if lin_flag == 0
        stim_env = log10(stim_env + 1.0);
    end
    
    % calculate stim avg
    %
    % Before Do calculation, we want to check if we got the correct input
    tempXsize = size(stim_env,1);
    if tempXsize ~= nband
        errordlg(['The selected new stimulus has different spatial dimension size',...
                ' from estimation data. Please reselect the new stimulus. '],...
            'Prediction Data Error', 'modal');
        close(tempWait);
        
        return;
    end
    
    stim_avg = stim_avg + sum(stim_env(1:nband, :), 2);
    count_avg = count_avg +(nt + 2*sil_window);  
    
    % clear workspace
    clear stim_env
    
end
stim_avg = stim_avg/count_avg;

for fidx = 1:ndata_files
    waitbar(n/ndata_files, tempWait);
    % load stimulus files
    stim_env = Check_And_Load(StimDS{fidx}.stimfiles);
    nlen = size(stim_env,2);
    forwardJN_s = handles.beststrf;
    nband = min(nband, size(forwardJN_s, 1));
   
    % subtract mean of stim from each stim
    stimval = zeros(nband, nlen);
    for tgood = 1:nlen
        stimval(:, tgood) = stim_env(1:nband, tgood) - stim_avg(1:nband);
    end

    % NEW algorithm
    est_spike{fidx} = zeros(nlen+end_window,1);
    clear tempResult
    
    for ib = 1:nband
        tempResult = conv(stimval(ib, :), forwardJN_s(ib, :));
        chopsize = size(forwardJN_s, 2);

        % Chop to our wanted length
        chopsize = (chopsize +1) /2;
        est_spike{fidx} = est_spike{fidx} +...
            tempResult(1,chopsize:size(tempResult, 2) - chopsize +1)';
    end     % END of ib

    % Clear workspace
    clear stimval;
    clear stim_env;

end
close(tempWait);
handles.predicted_estspike = est_spike;

guidata(handles.figure1, handles);

% --- Executes on slider movement.
function tolvalSlider_Callback(hObject, eventdata, handles)
global Tol_val
newSlider = round(get(hObject, 'Value'));
if newSlider <= 0
    set(handles.tolvalshow, 'String', Tol_val(handles.tolindex));
else
    newSlider = round(newSlider);
    set(handles.tolvalshow, 'String', Tol_val(newSlider));
    handles.Index = newSlider;
end
guidata(handles.figure1, handles);
get_strf(handles);
set(handles.predstim, 'visible', 'off');
child = get(handles.predstim, 'Children');
set(child, 'Visible', 'off')
set(handles.predpsth, 'visible', 'off');
child = get(handles.predpsth, 'Children');
set(child, 'Visible', 'off')


% --- Executes during object creation, after setting all properties.
function tolvalSlider_CreateFcn(hObject, eventdata, handles)

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
    handles.Index = newIndex;
    %set(hObject, 'String', num2str(handles.Index));
    % Update slider value based on dataSet value
    set(handles.tolvalSlider, 'String', newIndex);
end
guidata(handles.figure1, handles);

get_strf(handles);
set(handles.predstim, 'visible', 'off');
child = get(handles.predstim, 'Children');
set(child, 'Visible', 'off')
set(handles.predpsth, 'visible', 'off');
child = get(handles.predpsth, 'Children');
set(child, 'Visible', 'off')

% --- Executes during object creation, after setting all properties.
function tolvalshow_CreateFcn(hObject, eventdata, handles)
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
    handles.stdIndex = newSlider;
end
guidata(handles.figure1, handles);
get_strf(handles);
set(handles.predstim, 'visible', 'off');
child = get(handles.predstim, 'Children');
set(child, 'Visible', 'off')
set(handles.predpsth, 'visible', 'off');
child = get(handles.predpsth, 'Children');
set(child, 'Visible', 'off')

% --------------------------------------------------------------------
function get_beststrf(handles)
% --------------------------------------------------------------------
global Tol_val sDim outputPath
global Std_val DS

inforesult = load(fullfile(outputPath,'display_INFO_result.mat'));
ccresult = load(fullfile(outputPath,'info_r_result.mat'));

% JXZ 7/14/2005
infoV = zeros(length(Tol_val), length(Std_val));
for ii = 1:length(Tol_val)
    for jj = 1:length(Std_val)
        infoV(ii,jj) = inforesult.infopre{ii}{jj};
    end
end

% Get the largest predicted info value and index
infoMax = max(max(infoV));

% Get the position of the best infomax
[tolmax, stdmax] = find(infoV==infoMax);

% In order to remove redundent infoMax
tolmax = tolmax(1);
stdmax = stdmax(1);

% 5/03/2004: give a warning message if the best filter is with largest
% tol value or smallest tol value
% if tolmax == 1 | tolmax == length(Tol_val)
%     ttt = warndlg(['The best filter is found at boundaries of tolerance values.',...
%         ' You can change the tolerance ranges to redo the calculation.'],...
%         ' Tol Values Warning','modal');
%     uiwait(ttt);
% 
% end

% then load predicted STRF and display the best STRF
bestStrf = load(fullfile(outputPath,...
    ['strfResult_Tol',num2str(tolmax),'.mat']));

% find the best filter
% If only one data set is chosen for STRF estimation, we have to choose
% STRF instead of STRFJN and could not do smoothing on the STRF since STRFJNstd is zero.
if length(DS) == 1
    forward = bestStrf.STRF_Cell;
    [nb, nt] = size(forward);
  
else
    forwardJN = bestStrf.STRFJN_Cell;
    forwardJNstd = bestStrf.STRFJNstd_Cell;
    [nb, nt,nJN] = size(forwardJN);

    forwardJN_s = zeros(nb,nt,nJN);
    for iJN =1:nJN
        %  find the filtered JN STRF
        forwardJN_s(:,:,iJN) = fast_filter_filter(mean(forwardJN(:,:,[1:(iJN-1) (iJN+1):nJN]),3),...
            forwardJNstd(:,:,iJN), Std_val(stdmax));
    end

    % Take mean of Jackknifed STRFs
    forward = squeeze(mean(forwardJN_s,3));
end
handles.beststrf = forward;
handles.Nband = nb;
t = -(nt-1)/2:(nt-1)/2;
global initialFreq endFreq
if isempty(initialFreq) | isempty(endFreq)
    f = 1:nb;
    flabel = logspace(log10(1), log10(nb), nb);
    xtext = 'Channels';
else
    f = (linspace(initialFreq, endFreq, nb))/1000;
    flabel = logspace(log10(initialFreq), log10(endFreq), nb)/1000;
    xtext = 'Frequency (kHz)';
end
maxforward = max(max(forward));
minforward = min(min(forward));
absforward = max(abs(minforward),abs(maxforward));
axes(handles.predict_axes);
imagesc(t,f, forward, [-absforward absforward]); axis xy;
xlabel('Time'); ylabel(xtext)
title('Best STRF')
guidata(handles.figure1, handles);



% --- Executes on button press in displayresult.
function displayresult_Callback(hObject, eventdata, handles)
global StimDS
if ~isfield(handles, 'predicted_estspike')
    errordlg('No prediction results yet.', 'Pls predict first', 'modal');
    return;
end 

stim_env = Check_And_Load(StimDS{handles.fidx}.stimfiles);
set(handles.predstim, 'visible', 'on');
axes(handles.predstim);
imagesc(stim_env); axis xy;
axis tight;
title('Stim')

set(handles.predpsth, 'visible', 'on');
axes(handles.predpsth);
global ampsamprate
psthprediction = handles.predicted_estspike{handles.fidx}*ampsamprate;
plot(psthprediction); axis xy;
axis tight;
xlabel('Time');
title('prediction')

% Save to the file for NPC
global outputPath
save(fullfile(outputPath, 'NPC_prediction.mat'), 'psthprediction');

% --- Executes on slider movement.
function datasetSlider_Callback(hObject, eventdata, handles)
newSlider = round(get(hObject, 'Value'));
if newSlider <= 0
    set(handles.DataSet, 'String', 1);
else
    set(handles.DataSet, 'String', round(newSlider));
    handles.fidx = round(newSlider);
end
guidata(handles.figure1, handles);
global StimDS
[path, name, ext] = fileparts(StimDS{handles.fidx}.stimfiles);
set(handles.stimfile_show, 'String',[name ext]);

displayresult_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function datasetSlider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function DataSet_Callback(hObject, eventdata, handles)
global StimDS
maxlength = length(StimDS);
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

[path, name, ext] = fileparts(StimDS{handles.fidx}.stimfiles);
set(handles.stimfile_show, 'String',[name ext]);

displayresult_Callback(hObject, eventdata, handles);

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


% --------------------------------------------------------------------
function get_strf(handles)
% --------------------------------------------------------------------
global  sDim outputPath DS Std_val

% then load predicted STRF and display the best STRF
bestStrf = load(fullfile(outputPath,...
    ['strfResult_Tol',num2str(handles.Index),'.mat']));

% find the best filter
% If only one data set is chosen for STRF estimation, we have to choose
% STRF instead of STRFJN and could not do smoothing on the STRF since STRFJNstd is zero.
if length(DS) == 1
    forward = bestStrf.STRF_Cell;
    [nb, nt] = size(forward);
  
else
    forwardJN = bestStrf.STRFJN_Cell;
    forwardJNstd = bestStrf.STRFJNstd_Cell;
    [nb, nt,nJN] = size(forwardJN);

    forwardJN_s = zeros(nb,nt,nJN);
    for iJN =1:nJN
        %  find the filtered JN STRF
        forwardJN_s(:,:,iJN) = fast_filter_filter(mean(forwardJN(:,:,[1:(iJN-1) (iJN+1):nJN]),3),...
            forwardJNstd(:,:,iJN), Std_val(handles.stdIndex));
    end

    % Take mean of Jackknifed STRFs
    forward = squeeze(mean(forwardJN_s,3));
end
handles.beststrf = forward;
handles.Nband = nb;
guidata(handles.figure1, handles);
t = -(nt-1)/2:(nt-1)/2;
global initialFreq endFreq
if isempty(initialFreq) | isempty(endFreq)
    f = 1:nb;
    flabel = logspace(log10(1), log10(nb), nb);
    xtext = 'Channels';
else
    f = (linspace(initialFreq, endFreq, nb))/1000;
    flabel = logspace(log10(initialFreq), log10(endFreq), nb)/1000;
    xtext = 'Frequency (kHz)';
end
maxforward = max(max(forward));
minforward = min(min(forward));
absforward = max(abs(minforward),abs(maxforward));
axes(handles.predict_axes);
imagesc(t,f, forward, [-absforward absforward]); axis xy;
xlabel('Time'); ylabel(xtext)
title('STRF')




% --- Executes during object creation, after setting all properties.
function stdvalSlider_CreateFcn(hObject, eventdata, handles)
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
    handles.stdIndex = newIndex;
    %set(hObject, 'String', num2str(handles.Index));
    % Update slider value based on dataSet value
    set(handles.stdvalSlider, 'String', newIndex);
end
guidata(handles.figure1, handles);
get_strf(handles);
set(handles.predstim, 'visible', 'off');
child = get(handles.predstim, 'Children');
set(child, 'Visible', 'off')
set(handles.predpsth, 'visible', 'off');
child = get(handles.predpsth, 'Children');
set(child, 'Visible', 'off')

% --- Executes during object creation, after setting all properties.
function stdvalshow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
% END of PREDNDISPLAY_GUI Window
% --------------------------------------------------------------------
