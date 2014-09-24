function varargout = displaybestfit_GUI(varargin)
% DISPLAYBESTFIT_GUI Application M-file for displaybestfit_GUI.fig
%    FIG = DISPLAYBESTFIT_GUI launch displaybestfit_GUI GUI.
%    DISPLAYBESTFIT_GUI('callback_name', ...) invoke the named callback.
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
% Modified by JXZ 7/14/2005
%  Add more parameter display 'Std_val'
%  Add display option choice
%

if nargin == 0  % LAUNCH GUI

    fig = openfig(mfilename);
    %set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

    % for resize property
    %set(fig, 'resize', 'on');
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
        hUIControls],'fontname', 'Times New Roman','units','normalized',...
        'fontunits','normalized');

    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);

    % Set imageIndex to handles and for displaying
    handles.imageIndex = 0;
    % limited number of images shown on one screen
    handles.numList = 30;
    handles.nt = 0;

    guidata(fig, handles);

    % display best filter and related results
    for ii = 1:handles.numList
        axesStr = strcat('subaxes', num2str(ii-1));
        subhandle = getfield(handles, axesStr);
        set(subhandle, 'Visible', 'off');
        child = get(subhandle, 'Children');
        set(child, 'Visible', 'off');
        titlethis = get(subhandle, 'Title');
        set(titlethis, 'Visible', 'off');

    end

    % Check whether if we have done validation
    global Tol_val sDim ampsamprate predinitialFreq predendFreq
    global outputPath Std_val
    if isempty(outputPath) | isempty(Tol_val) | isempty(sDim) | isempty(Std_val)

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

                tempGV = fullfile(outputPath, 'GlobalVariables.mat');
                if not(exist(tempGV, 'file'))
                    errordlg('Wrong path for output files')
                    return
                end
                load(tempGV, 'Tol_val', 'sDim','ampsamprate','Std_val')

                % For 1-D display, we need two more global variables
                if strcmp(sDim, '1-D')
                    load(fullfile(outputPath, 'predVariables.mat'),...
                        'predinitialFreq','predendFreq')
                end


            otherwise
                msgbox('You need do prediction first', 'Warning', 'modal');
                return
        end
    end

    validate_file = fullfile(outputPath, 'info_r_result.mat');
    if ~exist(validate_file)
        errordlg('There is no validation result available for this input set.');
        return;
    end
    guidata(fig, handles);
    % Display two buttons if sDim = 2-D
    if strcmp(sDim, '1-D')
        set(handles.prevlag, 'visible','off')
        set(handles.nextlag, 'visible', 'off')
    end

    displayBestFit(handles);

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
%delete(handles.figure1);
close all;

% --------------------------------------------------------------------
function displayBestFit(handles)
% --------------------------------------------------------------------
% Get global variables
global Tol_val sDim outputPath predinitialFreq predendFreq ampsamprate
global preprocessOption Std_val DS

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
if tolmax == 1 | tolmax == length(Tol_val)
    ttt = warndlg(['The best filter is found at boundaries of tolerance values.',...
        ' You can change the tolerance ranges to redo the calculation.'],...
        ' Tol Values Warning','modal');
    uiwait(ttt);

end

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

    handles.nJN = nJN;

    forwardJN_s = zeros(nb,nt,nJN);
    for iJN =1:nJN
        %  find the filtered JN STRF
        forwardJN_s(:,:,iJN) = fast_filter_filter(mean(forwardJN(:,:,[1:(iJN-1) (iJN+1):nJN]),3), forwardJNstd(:,:,iJN), Std_val(stdmax));
    end

    % Take mean of Jackknifed STRFs
    forward = squeeze(mean(forwardJN_s,3));
end

handles.nb= nb;
handles.nt = nt;
% Now display it
maxforward = max(max(forward));
minforward = min(min(forward));
absforward = max(abs(minforward),abs(maxforward));

t = -(nt-1)/2:(nt-1)/2;
f = 1:nb;

if strcmp(sDim, '1-D')

    axes(handles.resp_axes);
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
    title('STRF');
else
    set(handles.resp_axes, 'visible', 'off');
    % 2-D display
    splitX = floor(sqrt(nb));
    STA = reshape(forward, splitX, splitX, nt);
    amax = max(abs(STA(:)));
    amin = -amax;

    numList = handles.numList;
    if numList > nt
        numList = nt;
    end
    global ampsamprate
    binsize = ceil((1/ampsamprate)*1000);
    titlestr = 'Latency(ms):';

    % display total 30 images
    for fr = handles.imageIndex*numList+1:numList+handles.imageIndex*numList-(nt-1)/2
        if fr > nt
            break
        end
        % display them in appropriate axes
        axesStr = strcat('subaxes', num2str(fr-handles.imageIndex * numList-1));
        axes(getfield(handles, axesStr));

        if amin~=amax,
            ttH = imagesc(STA(:,:,fr+(nt-1)/2),[amin,amax]);

        else
            ttH = imagesc(zeros(size(STA,1),size(STA,2)));
        end

        if size(STA,1)>1 & size(STA,2)>1,
            axis image
            set(get(ttH,'Parent'),'YTickLabel',[]);
            set(get(ttH,'Parent'),'XTickLabel',[]);
            %axis off
        end

        if (fr-handles.imageIndex * numList -1 ) ~= 0
            curTitle = num2str((fr-1)*binsize);
        else
            curTitle = 'Latency(ms)';
        end
        title(curTitle);
        colormap(redblue);
    end
end
loaded = load(fullfile(outputPath,'best_strf'));


% Show info value and cc value
set(handles.infopred_show, 'String', loaded.infopre);
set(handles.info_show, 'String', loaded.info);
set(handles.cc_show, 'String', loaded.cc_ratio_max);
set(handles.tolval_show, 'String', Tol_val(tolmax));
set(handles.stdval_show, 'String', Std_val(stdmax));
set(handles.ccc_show, 'String', loaded.cc_ratio_const);


guidata(handles.figure1, handles);


% --------------------------------------------------------------------
function varargout = prevlag_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------

handles.imageIndex = handles.imageIndex - 1;
if handles.imageIndex < 0
    handles.imageIndex = 0;
end

guidata(h, handles);
for ii = 1:handles.numList
    axesStr = strcat('subaxes', num2str(ii-1));
    subhandle = getfield(handles, axesStr);
    set(subhandle, 'Visible', 'off');
    child = get(subhandle, 'Children');
    set(child, 'Visible', 'off');
    titlethis = get(subhandle, 'Title');
    set(titlethis, 'Visible', 'off');

end

displayBestFit(handles);


% --------------------------------------------------------------------
function varargout = nextlag_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------

handles.imageIndex = handles.imageIndex +1;
if handles.imageIndex*handles.numList > handles.nt
    % update imageIndex
    handles.imageIndex = 0;
end

guidata(h, handles);

for ii = 1:handles.numList
    axesStr = strcat('subaxes', num2str(ii-1));
    subhandle = getfield(handles, axesStr);
    set(subhandle, 'Visible', 'off');
    child = get(subhandle, 'Children');
    set(child, 'Visible', 'off');
    titlethis = get(subhandle, 'Title');
    set(titlethis, 'Visible', 'off');

end
displayBestFit(handles);




% --------------------------------------------------------------------
function varargout = ccc_show_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.ccc_show.
disp('ccc_show Callback not implemented yet.')


% --- Executes during object creation, after setting all properties.
function displayoption_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');


% --- Executes on selection change in displayoption.
function displayoption_Callback(hObject, eventdata, handles)
v = get(handles.displayoption, 'value');
switch v
    case 1 % Display best STRF only
        displayBestFit(handles);

    case 2  % Display all STRFs in one window

        if handles.nJN == 1
            errordlg(['Since only one data set is chosen for estimation, ',...
                'we can not do smoothing on the STRF.']);
            return;
        end

        global Std_val Tol_val preprocessOption sDim
        global outputPath initialFreq endFreq ampsamprate
        nstd = length(Std_val);
        ntol = length(Tol_val);
        inforesult = load(fullfile(outputPath,'display_INFO_result.mat'));

        nx = handles.nb;
        nt = handles.nt;
        nJN = handles.nJN;

        t = -(nt-1)/2:(nt-1)/2;
        if isempty(initialFreq) | isempty(endFreq)
            f = 1:nx;
        else
            f = linspace(initialFreq, endFreq, nx);
        end
        forwardJN_s = zeros(nx, nt, nJN);
        if strcmp(sDim, '1-D')
            figure('name', 'Display all filtered STRFs');
        end
        for itol=1:ntol
            for istd=1:nstd
                % Get infopre values for all filtered filters
                infoV(itol,istd) = inforesult.infopre{itol}{istd};
            end
        end
        for itol=1:ntol

            % Load all the results for each tol value
            clear strfresult forwardJN forwardJNstd forwardJN_s
            strfresult = load(fullfile(outputPath,['strfResult_Tol',num2str(itol),'.mat']));
            forwardJN = strfresult.STRFJN_Cell;
            forwardJNstd = strfresult.STRFJNstd_Cell;
            forwardcat = [];
            for istd=1:nstd
                stdfilt = Std_val(istd);
                titlestr{istd} = sprintf('std = %3.2f\n',Std_val(istd));



                % Filter the filters
                for iJN =1:nJN
                    %  find the filtered JN STRF
                    forwardJN_s(:,:,iJN) = fast_filter_filter(mean(forwardJN(:,:,[1:(iJN-1) (iJN+1):nJN]),3), forwardJNstd(:,:,iJN), stdfilt);
                end

                forward = squeeze(mean(forwardJN_s, 3));
                maxforward = max(max(forward));
                minforward = min(min(forward));
                absforward = max(abs(minforward),abs(maxforward));

                if strcmp(sDim, '1-D')
                    % Now plotting all filtered strfs for all tol values
                    h = subplot(nstd,ntol,itol + (istd-1)*ntol);
                    if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
                        flabel = logspace(log10(initialFreq), log10(endFreq), nx);
                        ttH =pcolor(t*ceil(1000/ampsamprate), f/1000, forward); shading interp;
                        caxis([-absforward absforward]);
                    else
                        ttH = imagesc(t*ceil(1000/ampsamprate),f/1000,forward,[-absforward absforward]);
                        axis xy;
                    end

                    xlim([-20 80]);
                    set(get(ttH,'Parent'),'YTickLabel',[]);
                    set(get(ttH,'Parent'),'XTickLabel',[]);


                    curTitle = ['I=' num2str(.01*round(100*infoV(itol,istd))) ' Tol ' num2str(Tol_val(itol)) ' Std ' num2str(stdfilt) ]; %sprintf('I=%3.2f', infoV(itol,istd));
                    title(curTitle);
                    if infoV(itol,istd) == max(max(infoV));
                        set(h,'LineWidth',2);
                    end
                    drawnow;
                else
                    forwardcat = cat(3, forwardcat, forward);
                end

            end  % END of istd

            % Now display all 2-D smoothed STRFs
            if strcmp(sDim, '2-D')
                % display
                splitX = floor(sqrt(nx));
                binsize = ceil((1/ampsamprate)*1000);
                figure('name', ['All smoothed strfs for tol val = ', num2str(Tol_val(itol))]);
                showkern(forwardcat,'space',[splitX splitX],titlestr, 0, binsize);

            end
        end    % END of itol
end  % END of switch

% --------------------------------------------------------------------
function varargout = help_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
ttlStr = get(handles.figure1, 'Name');
hlpStr = [...
    '                                                                            '
    '        DISPLAY BEST FILTER WINDOW                                          '
    '                                                                            '
    ' The best filter is chosen with the largest predicted information value.    '
    ' The left panel shows the best strf plot. The right panel shows the related '
    ' information about this strf.                                               '
    ' It includes:                                                               '
    '      1. Tol value: what tol value is used for this strf.                   '
    '      2. Std value: what std value is used for this strf.                   '
    '      3. Info value: information value of the prediction data.              '
    '      4. Info Pred value: predicted information value of the prediction.    '
    '      5. max_CC_ratio: the max ratio between the predicted CC and CC among  '
    '         the range of smoothing window size.                                '
    '      6. const_CC_ratio: the ratio between the predicted CC and CC at fixed '
    '          smoothing window size, e.g. 21 ms.                                '
    '      7. Help: This window will pop up. For detailed doc, please refer      '
    '              to the user manual from http://strfpak.berkeley.edu.          '
    '      8. Close: close all the figures.                                      '
    '                                                                            '
    ' This window also can display all the estimated STRFs in separate windows.  '
    ' If input data is 1-D in spatio domain, all the STRFs are ploted in one     '
    ' figure with their best predicted info values. If it is 2-D case, there will'
    ' be more figures with the smoothed STRFs if more than one Tol value is used.'
    '                                                                            '
    ' Updated by Junli, Sept. 2005.                                              '
    '                                                                            '];
myFig = handles.figure1;
helpwin(hlpStr, ttlStr);


% --------------------------------------------------------------------
% END of DISPLAYBESTFIT_GUI
% --------------------------------------------------------------------


% --- Executes on button press in all_strf_new_window.
function all_strf_new_window_Callback(hObject, eventdata, handles)
% hObject    handle to all_strf_new_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.nJN == 1
    errordlg(['Since only one data set is chosen for estimation, ',...
        'we can not do smoothing on the STRF.']);
    return;
end

global Std_val Tol_val preprocessOption sDim
global outputPath initialFreq endFreq ampsamprate
nstd = length(Std_val);
ntol = length(Tol_val);
inforesult = load(fullfile(outputPath,'display_INFO_result.mat'));

nx = handles.nb;
nt = handles.nt;
nJN = handles.nJN;

t = -(nt-1)/2:(nt-1)/2;
if isempty(initialFreq) | isempty(endFreq)
    f = 1:nx;
else
    f = linspace(initialFreq, endFreq, nx);
end
forwardJN_s = zeros(nx, nt, nJN);
if strcmp(sDim, '1-D')
    figure('name', 'Display all filtered STRFs');
end

for itol=1:ntol
    for istd=1:nstd
        % Get infopre values for all filtered filters
        infoV(itol,istd) = inforesult.infopre{itol}{istd};
    end
end
for itol=1:ntol

    % Load all the results for each tol value
    clear strfresult forwardJN forwardJNstd forwardJN_s
    strfresult = load(fullfile(outputPath,['strfResult_Tol',num2str(itol),'.mat']));
    forwardJN = strfresult.STRFJN_Cell;
    forwardJNstd = strfresult.STRFJNstd_Cell;
    forwardcat = [];
    for istd=1:nstd
        stdfilt = Std_val(istd);
        titlestr{istd} = sprintf('std = %3.2f\n',Std_val(istd));



        % Filter the filters
        for iJN =1:nJN
            %  find the filtered JN STRF
            forwardJN_s(:,:,iJN) = fast_filter_filter(mean(forwardJN(:,:,[1:(iJN-1) (iJN+1):nJN]),3), forwardJNstd(:,:,iJN), stdfilt);
        end

        forward = squeeze(mean(forwardJN_s, 3));
        maxforward = max(max(forward));
        minforward = min(min(forward));
        absforward = max(abs(minforward),abs(maxforward));

        if strcmp(sDim, '1-D')
            % Now plotting all filtered strfs for all tol values
            h = subplot(nstd,ntol,itol + (istd-1)*ntol);
            if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
                flabel = logspace(log10(initialFreq), log10(endFreq), nx);
                ttH =pcolor(t*ceil(1000/ampsamprate), f/1000, forward); shading interp;
                caxis([-absforward absforward]);
            else
                ttH = imagesc(t*ceil(1000/ampsamprate),f/1000,forward,[-absforward absforward]);
                axis xy;
            end

            xlim([-20 80]);
            set(get(ttH,'Parent'),'YTickLabel',[]);
            set(get(ttH,'Parent'),'XTickLabel',[]);


            curTitle = ['I=' num2str(.01*round(100*infoV(itol,istd))) ' Tol ' num2str(Tol_val(itol)) ' Std ' num2str(stdfilt) ]; %sprintf('I=%3.2f', infoV(itol,istd));
            title(curTitle);
            if infoV(itol,istd) == max(max(infoV));
                set(h,'LineWidth',2);
            end
            drawnow;
        else
            forwardcat = cat(3, forwardcat, forward);
        end

    end  % END of istd

    % Now display all 2-D smoothed STRFs
    if strcmp(sDim, '2-D')
        % display
        splitX = floor(sqrt(nx));
        binsize = ceil((1/ampsamprate)*1000);
        figure('name', ['All smoothed strfs for tol val = ', num2str(Tol_val(itol))]);
        showkern(forwardcat,'space',[splitX splitX],titlestr, 0, binsize);

    end
end    % END of itol

