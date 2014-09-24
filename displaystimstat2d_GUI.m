function varargout = displaystimstat2d_GUI(varargin)
% DISPLAYSTIMSTAT2D_GUI Application M-file for displaystimstat_GUI.fig
%    FIG = DISPLAYSTIMSTAT2D_GUI launch displaystimstat_GUI GUI.
%    DISPLAYSTIMSTAT2D_GUI('callback_name', ...) invoke the named callback.
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
% Last Modified by GUIDE v2.5 06-Apr-2006 10:36:38
%
if nargin == 0  % LAUNCH GUI

    fig = openfig(mfilename);

    % Use system color scheme for figure:
    %set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

    % for resize property
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
        hUIControls],'fontname', 'Times New Roman','units','normalized','fontunits','normalized');

    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    guidata(fig, handles);

    %set(handles.stimauto, 'Visible', 'off')
    set(handles.stimspikecross, 'Visible', 'off')
    set(handles.corr_axes, 'Visible', 'off')
    set(handles.normalcorr_axes, 'Visible', 'off')

    % turn off total numbers of small axes for 2-d image
    handles.numList = 40;
    for ii = 1: handles.numList
        axesStr = strcat('subaxes', num2str(ii-1));
        set(getfield(handles, axesStr), 'Visible', 'off');
    end

    handles.imageIndex = 0;
    handles.nt = 0;
    guidata(fig, handles);

    % Load stim AutoCorrelation matrix
    global outputPath sDim ampsamprate DS TimeLag Tol_val NBAND
    if isempty(outputPath) | isempty(sDim)
        anws = questdlg('Have you done calculating stimulus second-order statistics?',...
            'Calculation Question','Yes', 'No', 'Yes');
        switch anws
            case 'Yes'
                prompt={['Enter the path of output data files:']};
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
                tempGV = fullfile(outputPath, 'Global_Variables.mat');
                if not(exist(tempGV, 'file'))
                    errordlg('Wrong path for output files')
                    return
                end
                load(tempGV,'sDim','ampsamprate','DS',...
                    'initialFreq', 'endFreq', 'TimeLag', 'Tol_val','NBAND')

            case 'No'
                msgbox('You need do calculation first', 'Warning', 'modal');
                return
        end
    end
    crossfile = fullfile(outputPath, 'StimResp_crosscorr.mat');
    if not(exist(crossfile))
        errordlg('No strf results yet, please do calculation first.','Results missing','modal');
        return;
    end
    handles.crossfile = crossfile;
    guidata(fig, handles);
    plot_stimspike(handles, crossfile);


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
function varargout = displaystimstat_option_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Get the display option to display
v = get(handles.displaystimstat_option, 'value');
switch v
    % Stim Spike Cross Matrix Display
    case 1
        %set(handles.stimauto, 'Visible', 'off')
        set(handles.stimspikecross, 'Visible', 'off')
        set(handles.corr_axes, 'Visible', 'off')
        set(handles.normalcorr_axes, 'Visible', 'off')

        axes(handles.corr_axes);
        handles.cross_only = 1;

        % Load stim spike crossCorr
        global outputPath
        crossfile = fullfile(outputPath, 'StimResp_crosscorr.mat');
        plot_stimspike(handles, crossfile);

        % Stim AutoCorr Display
    case 2

        %set(handles.stimauto, 'Visible', 'off')
        set(handles.stimspikecross, 'Visible', 'off')
        child = get(handles.stimspikecross, 'Children');
        set(child, 'Visible', 'Off')

        set(handles.corr_axes, 'Visible', 'off')
        child = get(handles.corr_axes,'Children');
        set(child, 'Visible', 'Off')

        set(handles.normalcorr_axes, 'Visible', 'off')
        child = get(handles.normalcorr_axes,'Children');
        set(child, 'Visible', 'Off')

        % display working signal
        set(handles.figure1,'Pointer', 'watch');
        OldBC=get(handles.displaystimstat_option, 'BackgroundColor');
        set(handles.displaystimstat_option, 'BackgroundColor', [1.0 0.5 0.5]);

        % Load stim AutoCorrelation matrix
        %axes(handles.corr_axes);

        global outputPath
        autoCorrfile = fullfile(outputPath, 'Stim_autocorr.mat');
        plot_stimauto( autoCorrfile);

        set(handles.displaystimstat_option, 'BackgroundColor', OldBC);
        set(handles.figure1, 'Pointer', 'Arrow');

        % Display both cross-correlation in separate window
        %     case 3
        %         global sDim ampsamprate
        %         if strcmp(sDim, '1-D')
        %             errordlg('This option is only for 2-D')
        %             return
        %         end
        %
        %         % Load stim spike crossCorr
        %         global outputPath
        %         crossfile = fullfile(outputPath, 'StimResp_crosscorr.mat');
        %         stim_spike = Check_And_Load(crossfile);
        %         [nb, nt] = size(stim_spike);
        %
        %         stim_spike = fliplr(stim_spike);
        %
        %         %smoothed stim_spike
        %         movx = sqrt(nb);
        %         w = hanning(nt);
        %
        %         for ib=1:nb
        %             stim_spike_nor(ib,:) = stim_spike(ib,:).*w';
        %         end
        %
        %         % display
        %         binsize = ceil((1/ampsamprate)*1000);
        %         tkern=cat(3,stim_spike,stim_spike_nor);
        %         titles={'Raw stimulus-spike cross correlation',...
        %            'Smoothed stimulus-spike cross correlation'};
        %         % display results
        %         figure;
        %         showkern(tkern,'space',[movx movx],titles, 0, binsize);

    case 3
        global setSep
        if setSep == 0
            global initialFreq endFreq fwidthHz ampsamprate NBAND
            if isempty(fwidthHz)
                fwidthHz = 125;
            end

            % Right now, we only display Modulation Spectrum for nonseparable
            %
            cal_ModulationSpectrum_plot(fwidthHz,ampsamprate,...
                fix((endFreq-initialFreq)/NBAND));
        else
            set(handles.displaystimstat_option, 'value', 1);
            errordlg(['We only display modulation spectrum for nonseparable',...
                ' calculation now.'], 'Modulation Spectrum Error', 'modal');
            return;
        end
end


% --------------------------------------------------------------------
function varargout = helpbutton_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------



% --------------------------------------------------------------------
function varargout = close_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
%delete(handles.figure1)
close all

% --------------------------------------------------------------------
function plot_stimspike(handles, crossfile)
% --------------------------------------------------------------------
% Loads and process the stim - Spike cross-correlation
stim_spike = Check_And_Load(crossfile);

disp('Done loading stim_spike cross correlatin.');

[nb, nt] = size(stim_spike);
w = hanning(nt);
stim_spike_nor = fliplr(stim_spike);
stim_spike = stim_spike_nor;
for ib=1:nb
    stim_spike_nor(ib,:)=stim_spike_nor(ib,:).*w';
end

global sDim ampsamprate
[nb, nt] = size(stim_spike);

handles.nt = nt;
guidata(handles.figure1, handles);

t = -(nt -1)/2:(nt -1)/2;


minval = min(min(stim_spike));
maxval = max(max(stim_spike));

% 2-D display
if strcmp(sDim, '2-D')

    % two-dimension display
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
        tsf(cfilt,:,:)=squeeze(sum(reshape(stim_spike_nor,chancount,phasecount,...
            tbincount,kcount),2));

        tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
        STA=reshape(tsf,Xmax,Xmax,tbincount);
    else
        % 2-D display
        splitX = floor(sqrt(nb));
        STA = reshape(stim_spike_nor, splitX, splitX, nt);
        
    end
    
    amax = max(abs(STA(:)));
    amin = -amax;

    numList = nt;
    if isempty(ampsamprate)
        binsize = 1;
        titlestr = 'Latency(frames)=';
    else
        binsize = ceil((1/ampsamprate)*1000);
        titlestr = 'Latency(msec)=';
    end


    % display total 12 images
    for fr = handles.imageIndex * numList +1:numList + handles.imageIndex * numList-(nt-1)/2

        if fr > nt |fr -handles.imageIndex*numList > handles.numList
            break
        end

        axesStr = strcat('subaxes', num2str(fr-handles.imageIndex * numList -1 ));
        axes(getfield(handles, axesStr));

        if amin~=amax,
            ttH = imagesc(STA(:,:,fr+(nt-1)/2),[amin,amax]);
        else
            ttH = imagesc(zeros(size(STA,1),size(STA,2)));
        end

        if size(STA,1)>1 & size(STA,2)>1,
            axis image
            set(get(ttH,'Parent'),'YTickLabel',[]);
            %set(get(ttH,'Parent'),'YTick',[]);
            set(get(ttH,'Parent'),'XTickLabel',[]);
            %set(get(ttH,'Parent'),'XTick',[]);
            %axis off
        end

        %curTitle = strcat(titlestr, num2str((fr-1-(nt-1)/2)*binsize));
        if (fr-handles.imageIndex * numList -1 ) ~= 0
            curTitle = num2str((fr-1)*binsize);
        else
            curTitle = titlestr;
        end
        title(curTitle);
        colormap(redblue);

    end

else
    % Make 2d frame button invisible
    set(handles.prev, 'Visible', 'off')
    set(handles.next, 'Visible', 'off')

    % 1-D display
    global initialFreq endFreq preprocessOption
    if isempty(initialFreq) | isempty(endFreq)
        f = 1:nb;
        label_yaxis = 'Frequency Band';
    else
        f = initialFreq:endFreq;
        label_yaxis = 'Frequency (Hz)';
    end

    if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
        flabel = logspace(log10(initialFreq), log10(endFreq), nb);
        axes(handles.corr_axes);
        pcolor(t*ceil(1000/ampsamprate), flabel, stim_spike); shading interp;
        %axes(handles.normalcorr_axes);
        ylabel(label_yaxis)
        title('Raw Cross Correlation');

        axes(handles.normalcorr_axes);
        pcolor(t*ceil(1000/ampsamprate), flabel, stim_spike_nor); shading interp;
        ylabel(label_yaxis)
        xlabel('Time (ms)');
        title('Smoothed Cross Correlation');
    else  % for spectrogram

        axes(handles.corr_axes);
        imagesc(t*ceil(1000/ampsamprate),f,stim_spike)
        axis xy;
        ylabel(label_yaxis)
        title('Raw Cross Correlation');

        axes(handles.normalcorr_axes);
        imagesc(t*ceil(1000/ampsamprate),f,stim_spike_nor)
        axis xy;
        ylabel(label_yaxis)
        xlabel('Time (ms)');
        title('Smoothed Cross Correlation');
    end

end


% --------------------------------------------------------------------
function plot_stimauto(autofile)
% --------------------------------------------------------------------
% Loads and process the stim - Spike cross-correlation
stim = Check_And_Load(autofile);
disp('Done loading stim auto correlatin.');

global setSep ampsamprate
if setSep == 0    % display stimuli auto-correlation for nonseparable case

    ncorr=size(stim,1);
    nt=size(stim,2);
    nt2=(nt-1)/2;
    t=-nt2:nt2;
    minval = min(min(stim));
    maxval = max(max(stim));

    % Set layout position parameter
    axisval(1)=t(1);
    axisval(2)=t(nt);
    axisval(3)=minval;
    axisval(4)=maxval;

    nb = (-1 + sqrt(1+8*ncorr))/2;

    % ASK user to input how much of auto-corr matrix he wants to display
    prompt={sprintf(['To display the whole auto-correlation matrix will take longer time\n',...
        'and get low-quality figures. For better view and performance,\n ',...
        'please enter the number of elements you want to display\n ',...
        '(The default value is total number of spatial parameters,\n',...
        'please enter a number less than or equal to the default value.)',...
        '                              \n'])};

    def = {num2str(nb)};
    dlgTitle='Number of auto-correlation display';
    lineNo=1;
    % picture feature
    AddOpts.Resize='on';
    AddOpts.WindowStyle='normal';
    AddOpts.Interpreter='tex';
    datadir =inputdlg2(prompt,dlgTitle,lineNo,def,AddOpts);

    % Check if user input is valid
    dis = str2double(datadir);
    if isempty(datadir)
        errordlg('You choose CANCEL, Bye.','Input Error','modal')
        return
    end
    if dis > nb
        errordlg('Please enter smaller number.','Input Error','modal')
        return
    end

    disp('Now displaying the large matrix.');
    figure('Name', sprintf('Upper %d  auto-correlation entries', dis*dis));
    clf

    RC = 0;
    xpp = 0;

    for RC = 0:dis

        for i=1:dis - RC
            subplot(dis, dis, i + RC * dis + RC)
            xpp = xpp +1;
            plot(t,stim(xpp,:))
            axis(axisval);
            axis off;

        end
        drawnow
        if dis < nb
            xpp = xpp + (nb -dis);
        end

    end
else    % for separable case

    figure('Name', 'Stimuli auto correlation matrix')
    lenT = size(stim{2}, 2);
    if isempty(ampsamprate)
        xlabelrange = (1:lenT);
        xtext = 'Time (in frames)';
    else
        xlabelrange = (1:lenT)*ceil(1000/ampsamprate);
        xtext = 'Time (ms)';
    end
    subplot(6,1,[1:5])
    imagesc(stim{1});
    axis xy; axis image;
    title('Second-order Spacial Correlation')
    ylabel('Spatial Channel')
    xlabel('Spatial Channel')
    subplot(6,1,6)
    plot(xlabelrange,stim{2});
    axis xy; axis tight;
    title('Temporal stimuli auto-correlation')
    xlabel(xtext)

end


disp('Done display the stimulus auto-correlation.');
% ===============================================================
% END of plot_stimauto
% ===============================================================

% --------------------------------------------------------------------
function varargout = next_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------

handles.imageIndex = handles.imageIndex +1;

if handles.imageIndex* handles.numList > handles.nt
    % update imageIndex
    handles.imageIndex = 0;
end
for ii = 1:handles.numList
    axesStr = strcat('subaxes', num2str(ii-1));
    subhandle = getfield(handles, axesStr);
    set(subhandle, 'Visible', 'off');
    child = get(subhandle, 'Children');
    set(child, 'Visible', 'off')
    titlethis = get(subhandle, 'Title');
    set(titlethis, 'Visible', 'off');
end
guidata(h, handles);
set(handles.displaystimstat_option, 'Value', 1);
plot_stimspike(handles, handles.crossfile);


% --------------------------------------------------------------------
function varargout = prev_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
handles.imageIndex = handles.imageIndex - 1;
if handles.imageIndex < 0
    handles.imageIndex = 0;
end
guidata(h, handles);
set(handles.displaystimstat_option, 'Value', 1);
plot_stimspike(handles, handles.crossfile);

% --------------------------------------------------------------------
% END of display stim stat window -JXZ
% --------------------------------------------------------------------


% --- Executes on button press in helpbutton.
function helpbutton_Callback(hObject, eventdata, handles)
% hObject    handle to helpbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ttlStr = get(handles.figure1, 'Name');
hlpStr = [...
    '                                                               '
    '             STRFPAK:  Display Stim Statistics Window          '
    '                                                               '
    ' Stim and response statistics include stim response cross      '
    ' correlation and stim auto correlation. This window can display'
    ' those statistics results.                                     '
    '                                                               '
    ' Display Option:                                               '
    '    Display stimulus spike cross correlation: This option is by'
    '    fault. Top panel shows raw cross correlation and the bottom'
    '    panel shows smoothed cross correlation. The smoothed cross '
    '    correlation is smoothed in time domain using hanning       '
    '    window.                                                    '
    '                                                               '
    '    Display stimulus auto correlation: When the user chooses   '
    '    this option, STRFPAK asks how many items the user want to  '
    '    display. The trick is to choose the smaller number because '
    '    displaying the whole auto matrix will take long long time. '
    '                                                               '
    '    Display modulation spectrum: The modulation spectrum is a  '
    '    representative group of sounds decomposed into its ripple  '
    '    components and the power density. For more information, see'
    '    the following paper:                                       '
    '                                                               '
    '    Singh NC and Theunissen FE, "Modulation spectra of natural '
    '    sounds and ethological theories of auditory processing",   '
    '    J Acoust Soc Am, 2003 Dec, 114 (6 Pt 1): 3394:411.         '
    '                                                               '
    '                                                               '
    '    Updated by Junli, May 2006.                                '];

myFig = handles.figure1;
helpwin(hlpStr, ttlStr);

