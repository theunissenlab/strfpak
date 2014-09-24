function varargout = displayinfocc_GUI(varargin)
%
% DISPLAYINFOCC_GUI Application M-file for displayinfocc_GUI.fig
%    FIG = DISPLAYINFOCC_GUI launch displayinfocc_GUI GUI.
%    DISPLAYINFOCC_GUI('callback_name', ...) invoke the named callback.
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

    fig = openfig(mfilename);
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
    global predDS predinitialFreq predendFreq
    global outputPath
    if ~isempty(outputPath)
        if ~exist(fullfile(outputPath,'best_strf.mat'),'file')
            
                compsave_Callback;
            
        end
    end
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

                % Rest part for tol_values
                tempGV = fullfile(outputPath, 'GlobalVariables.mat');
                if not(exist(tempGV, 'file'))
                    errordlg('Wrong path for output files')
                    return
                end
                load(tempGV, 'Tol_val', 'sDim', 'Std_val')

                % Load prediction and validate data set and their related parameters
                load(fullfile(outputPath, 'predVariables.mat'),'predDS',...
                    'predinitialFreq','predendFreq')

            otherwise
                msgbox('You need do validation first', 'Warning', 'modal');
                return
        end
    end

    % Assign tol_val and tol_val index to handles
    handles.tol_val = Tol_val;
    handles.std_val = Std_val;
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
    tol = handles.tol_val(handles.Index);
    set(handles.tolvalshow, 'String', tol);

    % display tol_val on the screen
    stds= handles.std_val(handles.stdIndex);
    set(handles.stdvalshow, 'String', stds);

    % display default option
    %displayinfocc(handles);

    % save figure's content
    guidata(fig, handles);

    if nargout > 0
        varargout{1} = fig;
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    %try
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    %     catch
    %         disp(lasterr);
    %     end

end


% --------------------------------------------------------------------
function varargout = displayoption_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% get option value
v = get(handles.displayoption, 'value');
switch v
    case 1  % display Correlation coefficience
        errordlg('Please choose display option', 'Display Option', 'modal');
        return;
    case 2

        handles.optionflg = 2;
        guidata(h, handles);

    case 3 % display info values
        handles.optionflg = 3;
        guidata(h, handles);
    case 4
        handles.optionflg = 4;
        guidata(h, handles);
end


% Based on optional values to display predicted info and cc values
displayinfocc(handles);

% --------------------------------------------------------------------
function varargout = prev_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
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

displayinfocc(handles);
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = next_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
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

displayinfocc(handles);
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = close_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
delete(handles.figure1);

% --------------------------------------------------------------------
function displayinfocc(handles)
% --------------------------------------------------------------------
% Assign the predict result to handles
global smoothVect  ampsamprate
global outputPath % This should be user input - but in ms not in number of points
if isempty(smoothVect)
    validate_file = fullfile(outputPath, 'info_r_result.mat');
    if ~exist(validate_file)
        set(handles.displayoption, 'value', 1);
        errordlg('Please compute goodness of fit first.');
        return;
    end
    load(validate_file,'smoothVect');

end
handles.cc = load(fullfile(outputPath, 'display_CC_result.mat'));
handles.info = load(fullfile(outputPath, 'display_INFO_result.mat'));

% get index and optionflag
% But index value is not useful for display all tol values option
index = handles.Index;
optionflag = handles.optionflg;

% get filter width range from global variable width or hard coded it

if isempty(smoothVect)
    width = 5:8:29;  % rhis is hard coded in two places - BIG NO NO!
else
    width = smoothVect(1):smoothVect(2):smoothVect(3);
end
%width = smoothVect*(1000.0/ampsamprate);   % Width in ms

% Get proper axes as current
set(handles.smooth, 'Visible', 'off')
child = get(handles.smooth, 'Children');
set(child, 'Visible', 'Off')
% Get proper axes as current
set(handles.ccratiotolval, 'Visible', 'off')
child = get(handles.ccratiotolval, 'Children');
set(child, 'Visible', 'Off')
% Get proper axes as current
set(handles.tolfixed, 'Visible', 'off')
child = get(handles.tolfixed, 'Children');
set(child, 'Visible', 'Off')

switch optionflag
    case 2
        % Get proper axes as current
        set(handles.smooth, 'Visible', 'off')
        child = get(handles.smooth, 'Children');
        set(child, 'Visible', 'Off')
        % Get proper axes as current
        set(handles.ccratiotolval, 'Visible', 'off')
        child = get(handles.ccratiotolval, 'Children');
        set(child, 'Visible', 'Off')
        % Get proper axes as current
        set(handles.tolfixed, 'Visible', 'off')
        child = get(handles.tolfixed, 'Children');
        set(child, 'Visible', 'Off')

        % get cc data
        cc = handles.cc;

        % plot on the axes
        axes(handles.predict_axes);
        plot(width, cc.cc_two_halves_corrval{index}{handles.stdIndex});
        hold on;
        plot(width, cc.cc_two_halves_corrval_low{index}{handles.stdIndex},'r');
        plot(width, cc.cc_two_halves_corrval_high{index}{handles.stdIndex},'g');
        plot(width, cc.cc_spike_pre_corrval{index}{handles.stdIndex},':');
        plot(width, cc.cc_spike_pre_corrval_low{index}{handles.stdIndex},'r:');
        plot(width, cc.cc_spike_pre_corrval_high{index}{handles.stdIndex},'g:');
        hold off
        legend('r', 'r\_low','r\_high','r predicted', 'r\_predicted\_low', 'r\_predicted\_high');
        ylabel('Corr Coeff (r)');

        % Display overall CC
        %set(handles.overallCC, 'String', cc.raw_r{index}{handles.stdIndex});

        % Plot CC ratio = prediced/orignal
        axes(handles.orignal_axes);

        % 8/17/2005: Avoid dividing by zero warning
        halvescorrval = cc.cc_two_halves_corrval{index}{handles.stdIndex};
        precorrval = cc.cc_spike_pre_corrval{index}{handles.stdIndex};
        notzeroind = find(cc.cc_two_halves_corrval{index}{handles.stdIndex}~=0);
        prehalves = zeros(size(halvescorrval));
        prehalves(notzeroind) = precorrval(notzeroind)./halvescorrval(notzeroind);

        plot(width,prehalves, 'k');
        xlabel('Smoothing Filter Width (ms)');
        ylabel('CC Ratio');

        %          twohalvesconstant = cc.cc_two_halves_constant{index}{handles.stdIndex};
        %          if twohalvesconstant == 0
        %              fixedccratio = 0.0;
        %          else
        %              fixedccratio = cc.cc_spike_pre_constant{index}{handles.stdIndex} ./twohalvesconstant;
        %          end
        global psth_smoothconst
        titlestring = sprintf('MaxCCRatio = %g at t = %g\nFixedCCRatio = %g at t = %g',...
            cc.cc_ratio_max{index}{handles.stdIndex}, cc.tmax_ratio{index}{handles.stdIndex}(1),...
            cc.cc_ratio_const{index}{handles.stdIndex}, psth_smoothconst);
        title(titlestring);

    case 3
        % Get proper axes as current
        set(handles.smooth, 'Visible', 'off')
        child = get(handles.smooth, 'Children');
        set(child, 'Visible', 'Off')
        % Get proper axes as current
        set(handles.ccratiotolval, 'Visible', 'off')
        child = get(handles.ccratiotolval, 'Children');
        set(child, 'Visible', 'Off')
        % Get proper axes as current
        set(handles.tolfixed, 'Visible', 'off')
        child = get(handles.tolfixed, 'Children');
        set(child, 'Visible', 'Off')

        % get info data
        info = handles.info;

        % Display overall CC
        %set(handles.overallCC, 'String', info.raw_r{index}{handles.stdIndex});

        % plot precited coherence
        axes(handles.predict_axes);
        plot(info.fpxypre{index}{handles.stdIndex}, info.cxypre{index}{handles.stdIndex}, info.fpxypre{index}{handles.stdIndex},...
            info.cxyuppre{index}{handles.stdIndex}, ':', info.fpxypre{index}{handles.stdIndex},...
            info.cxydownpre{index}{handles.stdIndex}, ':');
        global sDim
        if strcmp(sDim, '1-D')
            xlabel('Frequency (Hz)');
        end

        ylabel('Predicted coherence');
        xlim([0 50])
        plot_title = sprintf('predInfo = %g (predInfoLow=%g predInfoUp=%g)',...
            info.infopre{index}{handles.stdIndex},...
            info.infodownpre{index}{handles.stdIndex}, info.infouppre{index}{handles.stdIndex});
        title(plot_title);

        % plot coherence
        axes(handles.orignal_axes);
        plot(info.fpxypre{index}{handles.stdIndex}, info.cxy{index}{handles.stdIndex}, info.fpxypre{index}{handles.stdIndex},...
            info.cxyup{index}{handles.stdIndex}, ':', info.fpxypre{index}{handles.stdIndex}, info.cxydown{index}{handles.stdIndex}, ':');

        if strcmp(sDim, '1-D')
            xlabel('Frequency (Hz)');
        end

        ylabel('Coherence');
        xlim([0 50])
        plot_title = sprintf('Info = %g (InfoLow=%g InfoUp=%g) bits/s',...
            info.info{index}{handles.stdIndex}, info.infodown{index}{handles.stdIndex}, info.infoup{index}{handles.stdIndex});
        title(plot_title);

    case 4
        % Get proper axes as current
        set(handles.smooth, 'Visible', 'off')
        child = get(handles.smooth, 'Children');
        set(child, 'Visible', 'Off')
        % Get proper axes as current
        set(handles.ccratiotolval, 'Visible', 'off')
        child = get(handles.ccratiotolval, 'Children');
        set(child, 'Visible', 'Off')
        % Get proper axes as current
        set(handles.tolfixed, 'Visible', 'off')
        child = get(handles.tolfixed, 'Children');
        set(child, 'Visible', 'Off')

        % get info data
        info = handles.info;

        % Display overall CC
        %set(handles.overallCC, 'String', info.raw_r{index}{handles.stdIndex});

        numtol = length(handles.tol_val);
        rawr = zeros(numtol, 1);
        infor = zeros(numtol, 1);

        %axes(handles.ccratiotolval)
        for ii = 1:numtol
            rawr(ii) = info.raw_r{ii}{handles.stdIndex};
            infor(ii) = info.infopre{ii}{handles.stdIndex};
        end
        % plot precited info
        axes(handles.predict_axes);
        semilogx(handles.tol_val, infor, 'o--');
        xlabel('Tol Val')
        ylabel('Predicted info value');

        % plot predict cc
        axes(handles.orignal_axes);
        semilogx(handles.tol_val, rawr, 'o--');
        xlabel('Tol Val')
        ylabel('uncorrected r');

    case 4
        % Display CC_ration and Smooth Factor and Tol Values together
        % Get proper axes as current
        set(handles.predict_axes, 'Visible', 'off')
        child = get(handles.predict_axes, 'Children');
        set(child, 'Visible', 'Off')
        set(handles.orignal_axes, 'Visible', 'off')
        child = get(handles.orignal_axes, 'Children');
        set(child, 'Visible', 'Off')

        % get cc data
        cc = handles.cc;

        numtol = length(handles.tol_val);
        ccratio = zeros(numtol, 1);
        smoothfilter = zeros(numtol, 1);

        %axes(handles.ccratiotolval)
        for ii = 1:numtol
            ccratio(ii) = cc.cc_ratio_max{ii}{handles.stdIndex};
            smoothfilter(ii) = cc.tmax_ratio{ii}{handles.stdIndex};
        end

        % Display overall CC
        set(handles.overallCC, 'String', cc.raw_r{index}{handles.stdIndex});

        % Now displaying using ginput
        % Plot CC_ratio and tol value

        axes(handles.ccratiotolval)

        plot(handles.tol_val, ccratio);
        xlabel('Tol Val')
        ylabel('Max CC Ratio')
        %title('Max CC Ratio vs. Tol Val')

        axes(handles.smooth)
        plot(handles.tol_val, smoothfilter)
        %axis off
        xlabel('Tol Val')
        ylabel('Smooth filter width (ms)')
        %title('Smooth filter width vs. Tol Val')
        %legend(['Max CCRatio = ', sprintf('%4.2f', ccratio(index))]);

        chooseTol = index;
        axes(handles.tolfixed)
        plot(width, cc.cc_spike_pre_corrval{chooseTol}{handles.stdIndex} ./...
            cc.cc_two_halves_corrval{chooseTol}{handles.stdIndex})
        xlabel('Smoothing Filter Width (ms)');
        ylabel('CC Ratio');
        legend(['Tol Val = ', sprintf('%f',handles.tol_val(index))])

end



% --------------------------------------------------------------------
function varargout = tolval_show_ButtondownFcn(h, eventdata, handles, varargin)
% Stub for ButtondownFcn of the uicontrol handles.tolval_show.
disp('tolval_show ButtondownFcn not implemented yet.')


% --- Executes on button press in prevstdval.
function prevstdval_Callback(hObject, eventdata, handles)
% hObject    handle to prevstdval (see GCBO)
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

displayinfocc(handles);
guidata(hObject, handles);

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

displayinfocc(handles);
guidata(hObject, handles);

% --------------------------------------------------------------------
function varargout = help_button_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
ttlStr = get(handles.figure1, 'Name');
hlpStr = [...
    '                                                                            '
    '    STRFPAK: Computing goodness of fit and display window                   '
    '                                                                            '
    ' The window shows the validation results. The current version of STRFPAK    '
    ' implements two validation methods: correlation coeff and information value.'
    ' From this window, the user can display different results and visually      '
    ' compare them.                                                              '
    '                                                                            '
    '                                                                            '
    ' On Right Panel:                                                            '
    '     Display Option:                                                        '
    '        Choose display option:                                              '
    '            1. Corr Coeff(r): If the validation has been done, the user     '
    '                choose this option and the right panel will show the corr   '
    '                coeff results. If the validation has not been done, it gives'
    '                warning sign: Click the following "Compute goodness of fit" '
    '                first.                                                      '
    '            2. Info Values: If the validation has been done, the user       '
    '                choose this option and the right panel will show the info   '
    '                values results. If the validation has not been done, it     '
    '                gives warning sign: Click the following                     '
    '                "Compute goodness of fit" first.                            '
    '            3. Info/r vs Tol value: This option displays how info values    '
    '                vary with Tol values and how corr coeff (r) vary with Tol   '
    '                values. It is useful since it shows how these two validation'
    '                differ as Tol values changes.                               '
    '         Tol Val box: shows the tolerance value used for the STRF which do  '
    '                validation.                                                 '
    '         Tol val slider: use mouse can choose next or prev tol val          '
    '                  and the corresponding STRF and the right panel will show  '
    '                  the corresponding validation results.                     '
    '         Spareness box: shows the sparseness values used for the            '
    '                  STRF.                                                     '
    '         Spareness slider: use mouse can choose next or prev                '
    '               sparseness value and the corresponding STRF.                 '
    '                                                                            '
    '      Compute goodness of fit: This button should be clicked before choosing'
    '         different display option. But if the validation has been done, the '
    '         user dont need to click this button before choosing display option.'
    ' On Left Panel:                                                             '
    '      DISPLAY OPTION: Display CC                                            '
    '         1. Plot of Correlation coefficents with its lower and upper bound  '
    '            and predicted corr coeff with its lower and upper bound.        '
    '         2. Plot of Corr Coeff ratio as smoothing filter width (ms).        '
    '      DISPLAY OPTION: Display Info                                          '
    '         1. Plot of predicted coherence as frequency.                       '
    '         2. Plot of coherence as frequency.                                 '
    '                                                                            '
    '      DISPLAY OPTION: Info/corr. coeff vs. Tol val                          '
    '         1. Plot of Info values vs. Tol val.                                '
    '         2. Plot of unsmoothed corr coeff vs. Tol val.                      '
    '                                                                            '
    ' Updated by Junli, May 2006.                                                '
    '                                                                            '];

myFig = handles.figure1;
helpwin(hlpStr, ttlStr);

% --------------------------------------------------------------------
% END of DISPLAYINFOCC_GUI Window
% --------------------------------------------------------------------

% --- Executes on button press in compsave.
function compsave_Callback(hObject, eventdata, handles)

global outputPath;
% First check whether you have proper trial numbers (must > 1) for validation.
spike_psth = Check_And_Load(fullfile(outputPath,'spike_psth_count.mat'));
ntrials_proper=spike_psth(end);

if ntrials_proper <= 1
    global novalidation
    novalidation = 1;
    errordlg(['Since the most common trial number of prediction data is only one,',...
        'the current validation methods are not allowed. Sorry!']);
    return;
end

spike_psth=spike_psth(1:end-1);
ntrials_index=find(spike_psth>=ntrials_proper);
ntrials = spike_psth;

%fprintf(fid, 'ntrials_index = [%f];\n', ntrials_index);

% Variable binWindow used for validation
% The binWindow really depends on  width of autocorrelation of psth
global sDim Std_val Tol_val ampsamprate
if ~isempty(sDim)
    if strcmp(sDim, '1-D')
        binWindow = 128;
    else
        binWindow = 20;
    end
end

disp('Now computing goodness of fit.');
% JXZ, 8/23/2005
%  Ask the user to input smoothing range for cc
prompt={'Smallest smoothing window width (in ms):', 'Smoothing window width step (in ms):',...
    'Largest smoothing window width (in ms):','Fixed smoothing window width (in ms):'};
global smoothVect psth_smoothconst

if isempty(smoothVect)
    smoothVect = [5 8 29];
end
if isempty(psth_smoothconst)
    psth_smoothconst = 13;
end

def = {num2str(smoothVect(1)),num2str(smoothVect(2)),num2str(smoothVect(3)),num2str(psth_smoothconst)};
dlgTitle='Hanning window range for correlation coefficient';
lineNo=1;

% picture feature
AddOpts.Resize='on';
AddOpts.WindowStyle='normal';
AddOpts.Interpreter='tex';
smooth_V =inputdlg(prompt,dlgTitle,lineNo,def,AddOpts);
% show the waitbar to see how calculate process goes.
pause(.3);

smoothVect = round([str2double(smooth_V{1}) str2double(smooth_V{2}) str2double(smooth_V{3})]);
psth_smoothconst = round(str2double(smooth_V{4}));

%%%$$$ begin goodness_of_fit

% Loop through all the results for all the tol values
% JXZ 7/13/2005
%   Add extra loop for Std_val list

% Load part I of PSTH of predict response files
spredresult = fullfile(outputPath, 'predResult_avgSpike1.mat');
avgSpike1 = Check_And_Load(spredresult);

% Load part II of PSTH of predict response files
spredresult = fullfile(outputPath, 'predResult_avgSpike2.mat');
avgSpike2 = Check_And_Load(spredresult);
clear spike_psth;
global predDS
nrec = length(predDS);

% Load all predicted results and form them to a big matrix
spike_est1 = [];
spike_est2 = [];
for k=1:nrec

    %irec = ntrials_index(k);
    % Load avg_spike1.mat
    spike_est1temp = avgSpike1{k};

    % Load avg_spike2
    spike_est2temp = avgSpike2{k};

    spike_est1 = [spike_est1; spike_est1temp];
    spike_est2 = [spike_est2; spike_est2temp];
end

%figure; plot(spike_est1);

%figure; plot(spike_est2);

clear infopre infouppre infodownpre info...
    infoup infodown...
    cc_spike_pre_best cc_spike_pre_constant...
    cc_two_halves_tmax...
    cc_two_halves_constant tmax_pre...
    cc_ratio_max cc_ratio_consttmax_ratio...
    fpxypre cxy cxypre cxyup cxydown...
    cxyuppre cxydownpre...
    cc_spike_pre_corrval cc_spike_pre_corrval_low cc_spike_pre_corrval_high...
    cc_two_halves_corrval cc_two_halves_corrval_low cc_two_halves_corrval_high...
    raw_r;

global running_in_script_mode
%pack;
tempWait = waitbar(0,...
    'Computing goodness of fit, please wait... ');
%     patch([0 istd istd 0], [0 0 1 1], 'r');
global  allow_negative_rates
if isempty(allow_negative_rates)
    allow_negative_rates = 0;
end

for ntol = 1:length(Tol_val)
    spredresult = fullfile(outputPath,...
        sprintf('predResult_EstSpike_Tol%d.mat', ntol));
    estSpike = Check_And_Load(spredresult);
    for istd = 1:length(Std_val)
        waitbar((ntol-1)/length(Tol_val) + istd/length(Std_val)/length(Tol_val), tempWait);
        this_est = {};
        for kk = 1:length(estSpike)
            this_est{kk} = estSpike{kk}{istd};
        end
        % calculate Info and CC
        cal_V_checksum = checksum(load_function_text('cal_Validate'),this_est,spike_est1, spike_est2,ampsamprate, binWindow,smoothVect,psth_smoothconst,allow_negative_rates);
        if 1%~strcmp(running_in_script_mode,'yes')
            [infopre{ntol}{istd}, infouppre{ntol}{istd}, infodownpre{ntol}{istd}, info{ntol}{istd},...
                infoup{ntol}{istd}, infodown{ntol}{istd},...
                cc_spike_pre_best{ntol}{istd}, cc_spike_pre_constant{ntol}{istd},...
                cc_two_halves_tmax{ntol}{istd},...
                cc_two_halves_constant{ntol}{istd}, tmax_pre{ntol}{istd},...
                cc_ratio_max{ntol}{istd}, cc_ratio_const{ntol}{istd},tmax_ratio{ntol}{istd},...
                fpxypre{ntol}{istd}, cxy{ntol}{istd}, cxypre{ntol}{istd}, cxyup{ntol}{istd}, cxydown{ntol}{istd},...
                cxyuppre{ntol}{istd}, cxydownpre{ntol}{istd},...
                cc_spike_pre_corrval{ntol}{istd}, cc_spike_pre_corrval_low{ntol}{istd}, cc_spike_pre_corrval_high{ntol}{istd},...
                cc_two_halves_corrval{ntol}{istd},cc_two_halves_corrval_low{ntol}{istd}, cc_two_halves_corrval_high{ntol}{istd},...
                raw_r{ntol}{istd}]...
                = do_locally_cached_calc_checksum_known(get_local_cache_dir,'cal_Validate',cal_V_checksum,this_est,spike_est1, spike_est2,ntol, istd, ampsamprate, binWindow); %ampsamprate);
        else
            [infopre{ntol}{istd}, infouppre{ntol}{istd}, infodownpre{ntol}{istd}, info{ntol}{istd},...
                infoup{ntol}{istd}, infodown{ntol}{istd},...
                cc_spike_pre_best{ntol}{istd}, cc_spike_pre_constant{ntol}{istd},...
                cc_two_halves_tmax{ntol}{istd},...
                cc_two_halves_constant{ntol}{istd}, tmax_pre{ntol}{istd},...
                cc_ratio_max{ntol}{istd}, cc_ratio_const{ntol}{istd},tmax_ratio{ntol}{istd},...
                fpxypre{ntol}{istd}, cxy{ntol}{istd}, cxypre{ntol}{istd}, cxyup{ntol}{istd}, cxydown{ntol}{istd},...
                cxyuppre{ntol}{istd}, cxydownpre{ntol}{istd},...
                cc_spike_pre_corrval{ntol}{istd}, cc_spike_pre_corrval_low{ntol}{istd}, cc_spike_pre_corrval_high{ntol}{istd},...
                cc_two_halves_corrval{ntol}{istd},cc_two_halves_corrval_low{ntol}{istd}, cc_two_halves_corrval_high{ntol}{istd},...
                raw_r{ntol}{istd}]...
                = cal_Validate(this_est,spike_est1, spike_est2,ntol, istd, ampsamprate, binWindow); %ampsamprate);
        end

    end % END of ntol
end  % END of istd
%set(handles.figure1,'Pointer', 'arrow');
close(tempWait);
% set(handles.completed, 'visible', 'on');
% set(handles.inprogress, 'visible', 'off');

% save to the file for each data pair
% For displaying, save them to display files
save(fullfile(outputPath,'info_r_result.mat'), 'Tol_val', 'infopre', 'infouppre','smoothVect',...
    'infodownpre', 'info', 'infoup', 'infodown',...
    'cc_spike_pre_best', 'cc_spike_pre_constant','cc_ratio_const',...
    'cc_two_halves_tmax', 'cc_two_halves_constant', 'tmax_pre',...
    'cc_ratio_max', 'tmax_ratio');

save(fullfile(outputPath,'display_CC_result.mat'),'cc_spike_pre_corrval','cc_spike_pre_corrval_low','cc_spike_pre_corrval_high','raw_r',...
    'cc_two_halves_corrval','cc_two_halves_corrval_low','cc_two_halves_corrval_high','cc_ratio_max', 'tmax_ratio',...
    'cc_ratio_const');

save(fullfile(outputPath,'display_INFO_result.mat'),'fpxypre','cxy','cxyup','cxydown',...
    'cxypre', 'cxyuppre', 'cxydownpre', 'infopre', 'infouppre',...
    'infodownpre', 'info', 'infoup', 'infodown', 'raw_r');
global best_std_index finalDS predDS temp_DS DS

[strf,best_tol_index,best_std_index,savefile] = copy_best_strf;

tempWait = waitbar(0,...
    ['Now re-computing goodness of fit for the' char(10) 'best parameters using untouched data.']);
if ~isempty(finalDS)
    temp_DS = predDS;  %Switch predDS and finalDS
    predDS = finalDS;
end
nstim = length(estSpike);
global now_do_untouched Std_val NBAND
best_std = Std_val(best_std_index);
now_do_untouched = 'yes';

% Calculate stim_avg if not calculated before
stim_avg = cal_AVG(predDS, NBAND,1);

global TimeLagUnit TimeLag NBAND matchflg
if strcmp(TimeLagUnit, 'msec')
    global ampsamprate
    if isempty(ampsamprate)
        ampsamprate = 1000;
    end
    twindow = [-round(TimeLag*ampsamprate/1000) round(TimeLag*ampsamprate/1000)];
else
    twindow = [-TimeLag TimeLag];
end


if ~strcmp(running_in_script_mode,'yes')
    % Call cal_PredStrf to do prediction
    errFlg = cal_PredStrf(1,predDS, stim_avg,...
        fullfile(outputPath,['strfResult_Tol',num2str(best_tol_index),'.mat']),...
        best_tol_index, twindow, NBAND);
else
    errFlg = cal_PredStrf(1,predDS, stim_avg,...
        fullfile(outputPath,['strfResult_Tol',num2str(best_tol_index),'.mat']),...
        best_tol_index, twindow, NBAND,running_in_script_mode);
end

waitbar(1/3,tempWait);

clear infopre infouppre infodownpre info...
    infoup infodown...
    cc_spike_pre_best cc_spike_pre_constant...
    cc_two_halves_tmax...
    cc_two_halves_constant tmax_pre...
    cc_ratio_max cc_ratio_consttmax_ratio...
    fpxypre cxy cxypre cxyup cxydown...
    cxyuppre cxydownpre...
    cc_spike_pre_corrval cc_spike_pre_corrval_low cc_spike_pre_corrval_high...
    cc_two_halves_corrval cc_two_halves_corrval_low cc_two_halves_corrval_high...
    raw_r;

if matchflg == 0
    spredresult = fullfile(outputPath, 'finalResult_avgSpike1.mat');
    avgSpike1 = Check_And_Load(spredresult);

    % Load part II of PSTH of predict response files
    spredresult = fullfile(outputPath, 'finalResult_avgSpike2.mat');
    avgSpike2 = Check_And_Load(spredresult);
end
clear spike_psth;
global predDS
nrec = length(predDS);

% Load all predicted results and form them to a big matrix
spike_est1 = [];
spike_est2 = [];
for k=1:nrec

    %irec = ntrials_index(k);
    % Load avg_spike1.mat
    spike_est1temp = avgSpike1{k};

    % Load avg_spike2
    spike_est2temp = avgSpike2{k};

    spike_est1 = [spike_est1; spike_est1temp];
    spike_est2 = [spike_est2; spike_est2temp];
end

for ntol = best_tol_index
    spredresult = fullfile(outputPath,...
        sprintf('predResult_EstSpike_untouched.mat'));
    estSpike = Check_And_Load(spredresult);
    for istd = best_std_index
        this_est = {};
        for kk = 1:length(predDS)
            this_est{kk} = estSpike{kk}{best_std_index}; % This is the double-jackknifed prediction
        end
        % calculate Info and CC
        if 1 %~strcmp(running_in_script_mode,'yes')
            cal_V_checksum = checksum(load_function_text('cal_Validate'),now_do_untouched,this_est,spike_est1, spike_est2,ampsamprate, binWindow,smoothVect,psth_smoothconst);
            [infopre, infouppre, infodownpre, info,...
                infoup, infodown,...
                cc_spike_pre_best, cc_spike_pre_constant,...
                cc_two_halves_tmax,...
                cc_two_halves_constant, tmax_pre,...
                cc_ratio_max, cc_ratio_const,tmax_ratio,...
                fpxypre, cxy, cxypre, cxyup, cxydown,...
                cxyuppre, cxydownpre,...
                cc_spike_pre_corrval, cc_spike_pre_corrval_low, cc_spike_pre_corrval_high,...
                cc_two_halves_corrval,cc_two_halves_corrval_low, cc_two_halves_corrval_high,...
                raw_r]...
                = do_locally_cached_calc_checksum_known(get_local_cache_dir,'cal_Validate',cal_V_checksum,this_est,spike_est1, spike_est2,ntol, istd, ampsamprate, binWindow); %ampsamprate);
        else
            [infopre, infouppre, infodownpre, info,...
                infoup, infodown,...
                cc_spike_pre_best, cc_spike_pre_constant,...
                cc_two_halves_tmax,...
                cc_two_halves_constant, tmax_pre,...
                cc_ratio_max, cc_ratio_const,tmax_ratio,...
                fpxypre, cxy, cxypre, cxyup, cxydown,...
                cxyuppre, cxydownpre,...
                cc_spike_pre_corrval, cc_spike_pre_corrval_low, cc_spike_pre_corrval_high,...
                cc_two_halves_corrval,cc_two_halves_corrval_low, cc_two_halves_corrval_high,...
                raw_r]...
                = cal_Validate(this_est,spike_est1, spike_est2,ntol, istd, ampsamprate, binWindow); %ampsamprate);
        end

    end % END of ntol
end  % END of istd

save(savefile,  'infopre', 'infouppre','smoothVect', 'psth_smoothconst', ...
    'infodownpre', 'info', 'infoup', 'infodown',...
    'cc_spike_pre_best', 'cc_spike_pre_constant','cc_ratio_const',...
    'cc_two_halves_tmax', 'cc_two_halves_constant', 'tmax_pre',...
    'cc_ratio_max', 'tmax_ratio','-APPEND');
close(tempWait);
now_do_untouched = 'no';

if ~isempty(finalDS)
    predDS = temp_DS; %Switch predDS back
end

%%%$$$ end goodness_of_fit
tt = msgbox(['Done fit evaluation.  Best STRF and validation info' char(10) 'using untouched data have been saved to the file:' char(10) savefile], 'modal');
uiwait(tt);
global outputPath
save(fullfile(outputPath,'STRFPAK_script_parameters.mat'), 'smoothVect' , 'psth_smoothconst' ,'binWindow', '-APPEND');
add_to_STRFPAK_script('displayinfocc_GUI.m','goodness_of_fit');

if exist('handles','var')
    set(handles.figure1, 'Pointer', 'Arrow');
end
disp('Done measurement of goodness of fit.')

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
displayinfocc(handles);

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
displayinfocc(handles);

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
displayinfocc(handles);

% --- Executes during object creation, after setting all properties.
function stdvalSlider_CreateFcn(hObject, eventdata, handles)

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
    handles.stdIndex = newIndex;
    %set(hObject, 'String', num2str(handles.Index));
    % Update slider value based on dataSet value
    set(handles.stdvalSlider, 'String', newIndex);
end
guidata(handles.figure1, handles);
displayinfocc(handles);

% --- Executes during object creation, after setting all properties.
function stdvalshow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


