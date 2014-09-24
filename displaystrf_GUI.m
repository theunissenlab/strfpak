function varargout = displaystrf_GUI(varargin)
% DISPLAYSTRF_GUI Application M-file for displaystrf_GUI.fig
%    FIG = DISPLAYSTRF_GUI launch displaystrf_GUI GUI.
%    DISPLAYSTRF_GUI('callback_name', ...) invoke the named callback.
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
% Last Modified by GUIDE v2.5 06-Apr-2006 13:33:50
%
% Modifiedy by JXZ, 7/13/2005.
%      1. Add global variable 'Std_val' for smoothing STRFs
%      2. Add function 'plot_smoothedstrf'.
%      3. Add one more display option "display all strfs".
%

if nargin == 0  % LAUNCH GUI

    fig = openfig(mfilename);

    % for resize property
    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
        hUIControls],'fontname', 'Times New Roman','units','normalized','fontunits','normalized');

    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);

    % Initialize index of picture display
    handles.Index = 1;
    handles.imageIndex = 0;
    handles.nt = 0;
    handles.numList = 35;

    set(handles.small_strf, 'Visible', 'off')
    set(handles.small_strf_x, 'Visible', 'off')
    set(handles.small_strf_y, 'Visible', 'Off')
    set(handles.small_mod_strf, 'Visible', 'off')
    set(handles.STRFaxes, 'visible', 'off');
    set(handles.STAaxes, 'visible', 'off');
    set(handles.strfonly, 'visible', 'off')


    % total numbers of image on the screen
    for ii = 1:handles.numList
        axesStr = strcat('subaxes', num2str(ii-1));
        subhandle = getfield(handles, axesStr);
        set(subhandle, 'Visible', 'off');
    end

    % Need save all handles to the data
    guidata(fig, handles);

    global Tol_val sDim ampsamprate initialFreq endFreq TimeLag DS NBAND
    global outputPath Std_val
    if isempty(outputPath) | isempty(Tol_val) | isempty(sDim)| isempty(Std_val)
        anws = questdlg('Have you done estimation of STRF?',...
            'Calculation Question','Yes', 'No', 'Yes');
        switch anws
            case 'Yes'
                prompt={['Enter the path of output results:']};
                defaultdir = fullfile(pwd, 'Output');
                def={defaultdir};
                dlgTitle='Please Enter Path for Output';
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

                % Load old copy of global variables
                tempGV = fullfile(outputPath, 'GlobalVariables.mat');
                if not(exist(tempGV, 'file'))
                    errordlg('Wrong path for output files')
                    return
                end
                load(tempGV,'Tol_val', 'Std_val','sDim','ampsamprate','DS','TimeLag',...
                    'NBAND','initialFreq','endFreq')

            case 'No'
                msgbox('You need do calculation first', 'Warning', 'modal');
                return
        end
    end

    handles.tolindex = 1;
    handles.stdindex = 1;
    handles.imageIndex = 0;
    handles.numList = 35;  % subaxes for 2-D strf

    global Std_val Tol_val
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
    guidata(fig, handles);

    % Need save all handles to the data
    displayall(handles);

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
function varargout = displayoption_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
displayall(handles);

% --------------------------------------------------------------------
function displayall(handles)
% --------------------------------------------------------------------
global outputPath DS Tol_val Std_val
strfresultfile = fullfile(outputPath,...
    ['strfResult_Tol',num2str(handles.tolindex),'.mat']);
if not(exist(strfresultfile))
    errordlg('No strf results yet, please do calculation first.','Results missing','modal');
    return;
end
strf_result = load(strfresultfile);

% Based on displayoption to display
v = get(handles.displayoption, 'value');
switch v
    case 1 % STRF only

        plot_strfonly( handles);

    case 2  % Display STRF and Smoothed STRFs for each Tol

        % Check if only choosing one data set for calculation, smoothing is
        % not available since it need base on JN results.
        if length(DS) ==1
            set(handles.displayoption, 'value', 1);
            errordlg('Only one data set is chosen, then can not smooth the STRF.')
            return;
        end
        forwardJN = strf_result.STRFJN_Cell;
        forwardJNstd = strf_result.STRFJNstd_Cell;
        plot_smoothedstrf(forwardJN, forwardJNstd, handles);

    case 3  % Display STRF and Smoothed STRFs for all

        if length(DS) ==1 
            set(handles.displayoption, 'value', 1);
            errordlg('Only one data set is chosen, then can not smooth the STRF.')
            return;
        end
        
        global Std_val Tol_val preprocessOption initialFreq endFreq ampsamprate sDim
        nstd = length(Std_val);
        ntol = length(Tol_val);
        
        [nx, nt, nJN]= size(strf_result.STRFJN_Cell);
        clear strf_result;
        pack
        t = -(nt-1)/2:(nt-1)/2;
        if isempty(initialFreq) | isempty(endFreq)
            f = 1:nx;
            flabel = logspace(log10(1), log10(nx), nx);
            xtext = 'Frequency Band';
        else
            f = (linspace(initialFreq, endFreq, nx))/1000;
            flabel = logspace(log10(initialFreq), log10(endFreq), nx)/1000;
            xtext = 'Frequency (kHz)';
        end
        
        % Try to save cache space in memory
        %forwardJN_s = zeros(nx, nt, nJN, nstd, ntol);
        forwardJN_s = zeros(nx, nt, nJN);
        if strcmp(sDim, '1-D') | strcmp(sDim, '0-D')
            figure('name', 'Display STRFs as Tol val and Std val(sparsness)');
        end
        for itol=1:ntol
            
            % Load all the results for each tol value
            clear strfresult forwardJN forwardJNstd ttH maxforward minforward absforward forward
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
                    forwardJN_s(:,:,iJN) = fast_filter_filter(mean(forwardJN(:,:,[1:(iJN-1) (iJN+1):nJN]),3),...
                        forwardJNstd(:,:,iJN), stdfilt);  
                end
                forward = squeeze(mean(forwardJN_s, 3));
                maxforward = max(max(forward));
                minforward = min(min(forward));
                absforward = max(abs(minforward),abs(maxforward));
                
                % Display based on spatial dimension size
                if strcmp(sDim, '1-D')
                    
                    % Now plotting all filtered strfs for all tol values
                    subplot(nstd,ntol,itol + (istd-1)*ntol)
                    if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')         
                        ttH =pcolor(ceil(t*1000/ampsamprate), flabel, forward(:,:)); shading interp;
                        caxis([-absforward absforward]);   
                    else
                       
                        ttH = imagesc(ceil(t*1000/ampsamprate),f,forward(:,:),[-absforward absforward]);
                        axis xy;  axis tight;
                    end
                    
                    xlim([-20 80]);
                    set(get(ttH,'Parent'),'YTickLabel',[]);
                    set(get(ttH,'Parent'),'XTickLabel',[]);       
                    
                    curTitle = sprintf('T=%4.3f,S=%3.2f', Tol_val(itol),Std_val(istd));
                    title(curTitle);
                    drawnow;
                elseif strcmp(sDim, '2-D')
                    forwardcat = cat(3,forwardcat, forward);
                elseif strcmp(sDim, '0-D')
                     % Now plotting all filtered strfs for all tol values
                     subplot(nstd,ntol,itol + (istd-1)*ntol)

                     ttH = plot(ceil(t*1000/ampsamprate),forward(:,:));
                     axis tight;
                     set(get(ttH,'Parent'),'YTickLabel',[]);
                     set(get(ttH,'Parent'),'XTickLabel',[]);

                     curTitle = sprintf('T=%4.3f,S=%3.2f', Tol_val(itol),Std_val(istd));
                     title(curTitle);
                     drawnow;
                end
                
            end  % END of istd
            
            % Now display all 2-D smoothed STRFs
            if strcmp(sDim, '2-D')
                % display
                splitX = floor(sqrt(nx));
                global ampsamprate
                if isempty(ampsamprate)
                    binsize = 1;
                else
                    binsize = ceil((1/ampsamprate)*1000);
                end
                 figure('name', ['All smoothed strfs for tol val = ', num2str(Tol_val(itol))]);
                global preprocessOption
                if strcmp(preprocessOption, 'Fourier Power transform')
                    showkern(forwardcat,'pfft',[splitX splitX],titlestr, 0, binsize);
                else
                    showkern(forwardcat,'space',[splitX splitX],titlestr, 0, binsize);
                end
                
            end
        end      % END of itol
        
    case 4  % STRF with projection on both axes
        global sDim
        if strcmp(sDim, '1-D')
            strfonly = strf_result.STRF_Cell;
            strfonly_std = strf_result.STRFJNstd_Cell;
            plot_strfProjection(strfonly,strfonly_std, handles);
        else
            set(handles.displayoption, 'value', 1);
            errordlg('This option is only for 1-D')
            return
        end

    case 5   % STRF and STA

        plot_strfSTA(strf_result.STRF_Cell, handles);
    case 6  % Display STRFs for all tol values
        
        global Std_val Tol_val preprocessOption initialFreq endFreq ampsamprate sDim
        nstd = length(Std_val);
        ntol = length(Tol_val);
        
        [nx, nt, nJN]= size(strf_result.STRFJN_Cell);
        clear strf_result;
        pack
        t = -(nt-1)/2:(nt-1)/2;
        if isempty(initialFreq) | isempty(endFreq)
            f = 1:nx;
            flabel = logspace(log10(1), log10(nx), nx);
            xtext = 'Frequency Band';
        else
            f = (linspace(initialFreq, endFreq, nx))/1000;
            flabel = logspace(log10(initialFreq), log10(endFreq), nx)/1000;
            xtext = 'Frequency (kHz)';
        end
        
        % Try to save cache space in memory
        %forwardJN_s = zeros(nx, nt, nJN, nstd, ntol);
        forwardJN_s = zeros(nx, nt, nJN);
        if strcmp(sDim, '1-D') | strcmp(sDim, '0-D')
            figure('name', 'Display all tolvaled STRFs  ');
        end
        for itol=1:ntol

            % Load all the results for each tol value
            clear strfresult forwardJN forwardJNstd ttH maxforward minforward absforward forward
            strfresult = load(fullfile(outputPath,['strfResult_Tol',num2str(itol),'.mat']));
            forward = strfresult.STRF_Cell;
            maxforward = max(max(forward));
            minforward = min(min(forward));
            absforward = max(abs(minforward),abs(maxforward));
            forwardcat = [];


            if strcmp(sDim, '1-D')

                % Now plotting all filtered strfs for all tol values
                subplot(ntol/2,2,itol)
                if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
                    ttH =pcolor(ceil(t*1000/ampsamprate), flabel, forward(:,:)); shading interp;
                    caxis([-absforward absforward]);
                    xlim([-20 100]);
                    set(get(ttH,'Parent'),'YTickLabel',[]);
                    set(get(ttH,'Parent'),'XTickLabel',[]);
                    %
                else
                    ttH = imagesc(ceil(t*1000/ampsamprate),f,forward(:,:),[-absforward absforward]);
                    axis xy;  axis tight;
                    xlim([-20 100]);
                    set(get(ttH,'Parent'),'YTickLabel',[]);
                    set(get(ttH,'Parent'),'XTickLabel',[]);
                    %
                end

                curTitle = sprintf('TolVal=%f', Tol_val(itol));
                title(curTitle);
                drawnow;
            elseif strcmp(sDim, '2-D')
                forwardcat = cat(3,forwardcat, forward);
            elseif strcmp(sDim, '0-D')
                % Now plotting all filtered strfs for all tol values
                subplot(ntol/2,2,itol)
                plot(ceil(t*1000/ampsamprate),forward);
                axis tight;
                curTitle = sprintf('TolVal=%f', Tol_val(itol));
                title(curTitle);
                drawnow;
            end


            % Now display all 2-D smoothed STRFs
            if strcmp(sDim, '2-D')
                % display
                splitX = floor(sqrt(nx));
                global ampsamprate
                if isempty(ampsamprate)
                    binsize = 1;
                else
                    binsize = ceil((1/ampsamprate)*1000);
                end
                figure('name', ['All smoothed strfs for tol val = ', num2str(Tol_val(itol))]);
                global preprocessOption
                if strcmp(preprocessOption, 'Fourier Power transform')
                    showkern(forwardcat,'pfft',[splitX splitX],titlestr, 0, binsize);
                else
                    showkern(forwardcat,'space',[splitX splitX],titlestr, 0, binsize);
                end

            end
        end      % END of itol

end

set(handles.tolvaltext, 'visible', 'on');
set(handles.stdvaltext, 'visible', 'on');
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
set(handles.outdirBrowser, 'visible', 'on');

% --------------------------------------------------------------------
function varargout = displaystrf_close_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% clean other separate figures
close all;
%delete(handles.figure1);


% --------------------------------------------------------------------
function plot_strfonly(handles)
% --------------------------------------------------------------------

% Get calc parameter value from global variables
global sDim initialFreq endFreq ampsamprate preprocessOption
global Std_val

% Get proper axes as current
set(handles.small_strf, 'Visible', 'off')
child = get(handles.small_strf, 'Children');
set(child, 'Visible', 'Off')
set(handles.small_strf_x, 'Visible', 'off')
child = get(handles.small_strf_x, 'Children');
set(child, 'Visible', 'Off')
set(handles.small_strf_y, 'Visible', 'Off')
child = get(handles.small_strf_y, 'Children');
set(child, 'Visible', 'Off')
set(handles.small_mod_strf, 'Visible', 'off')
child = get(handles.small_mod_strf, 'Children');
set(child, 'Visible', 'Off')
set(handles.STRFaxes, 'visible', 'off');
child = get(handles.STRFaxes, 'Children');
set(child, 'Visible', 'Off')
itsTitle = get(handles.STRFaxes, 'Title');
set(itsTitle, 'Visible', 'Off')
set(handles.STAaxes, 'visible', 'off');
child = get(handles.STAaxes, 'Children');
set(child, 'Visible', 'Off')
itsTitle = get(handles.STAaxes, 'Title');
set(itsTitle, 'Visible', 'Off')
set(handles.strfonly, 'visible', 'off');
child = get(handles.strfonly, 'Children');
set(child, 'Visible', 'Off')
itsTitle = get(handles.strfonly, 'Title');
set(itsTitle, 'Visible', 'Off')

% total numbers of image on the screen
for ii = 1:handles.numList
    axesStr = strcat('subaxes', num2str(ii-1));
    subhandle = getfield(handles, axesStr);
    set(subhandle, 'Visible', 'off');
    child = get(subhandle, 'Children');
    set(child, 'Visible', 'off')

end

global outputPath
strf_result = load(fullfile(outputPath,...
    ['strfResult_Tol',num2str(handles.tolindex),'.mat']));
forwardJN = strf_result.STRFJN_Cell;
forwardJNstd = strf_result.STRFJNstd_Cell;
[nx, nt, nJN]= size(forwardJN);
t = -(nt-1)/2:(nt-1)/2;

% Filter the filters
forwardJN_s = zeros(nx, nt, nJN);
stdfilt = Std_val(handles.stdindex);
if nJN == 1
    forwardJN_s = strf_result.STRF_Cell;
else
    for iJN =1:nJN
        %  find the filtered JN STRF
        forwardJN_s(:,:,iJN) = fast_filter_filter(mean(forwardJN(:,:,[1:(iJN-1) (iJN+1):nJN]),3),...
            forwardJNstd(:,:,iJN), stdfilt);
    end
end
forward = squeeze(mean(forwardJN_s, 3));
maxforward = max(max(forward));
minforward = min(min(forward));
absforward = max(abs(minforward),abs(maxforward));

handles.nt = nt;
guidata(handles.figure1,handles);

% Display estimated STRF based on sptial dimension
if strcmp(sDim, '1-D')

    global initialFreq endFreq
    if isempty(initialFreq) | isempty(endFreq)
        f = 1:nx;
        flabel = logspace(log10(1), log10(nx), nx);
        xtext = 'Frequency Band';
    else
        f = (linspace(initialFreq, endFreq, nx))/1000;
        flabel = logspace(log10(initialFreq), log10(endFreq), nx)/1000;
        xtext = 'Frequency (kHz)';
    end

    % 1-D display
    set(handles.strfonly, 'visible', 'on')
    axes(handles.strfonly);

    if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')

        pcolor(ceil(t*1000/ampsamprate), flabel, forward); shading interp;
        caxis([-absforward absforward]);
        axis([ceil(t(1)*1000/ampsamprate) ceil(t(end)*1000/ampsamprate) flabel(1) flabel(end)])
    else
        imagesc(ceil(t*1000/ampsamprate),f,forward,[-absforward absforward]);
        axis xy; axis tight;
        axis([ceil(t(1)*1000/ampsamprate) ceil(t(end)*1000/ampsamprate) f(1) f(end)])
    end
    %xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
    xlabel('Time (ms)')
    ylabel(xtext);

    title('STRF');
elseif strcmp(sDim, '2-D')

    % Set other axes to invisible
    set(handles.strfonly, 'visible', 'off');
    global preprocessOption
    if strcmp(preprocessOption, 'Fourier Power transform')
        phasecount=1;
        tbincount = nt;
        kcount = 1;
        spacebincount = nx;
        chancount=spacebincount/phasecount;
        Xmax=sqrt(chancount*2);
        [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);

        tsf=zeros(Xmax*Xmax,tbincount,kcount);
        tsf(cfilt,:,:)=squeeze(sum(reshape(forward,chancount,phasecount,...
            tbincount,kcount),2));

        tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
        STA=reshape(tsf,Xmax,Xmax,tbincount);
    else
        % 2-D display
        splitX = floor(sqrt(nx));
        STA = reshape(forward, splitX, splitX, nt);

    end
    
    amax = max(abs(STA(:)));
    amin = -amax;
    numList = nt;

    global stimsamprate
    if isempty(stimsamprate)
        binsize = 1;
        titlestr = 'Latency(frame)';
    else
        binsize = ceil((1/stimsamprate)*1000);
        titlestr = 'Latency(msec)';
    end

    % display total 35 images on one screen
    for fr = handles.imageIndex * numList+1:numList+handles.imageIndex * numList-(nt-1)/2

        if fr > nt | fr -handles.imageIndex*numList > handles.numList
            break
        end

        % display them in appropriate axes
        axesStr = strcat('subaxes', num2str(fr-handles.imageIndex * numList-1));
        axes(getfield(handles, axesStr));

        % Display them in specific format

        if amin~=amax,
            ttH = imagesc(STA(:,:,fr+(nt-1)/2),[amin,amax]);

        else
            %ttH = imagesc(zeros(size(STA,1),size(STA,2)));
            ttH = imagesc(STA(:,:,fr+(nt-1)/2));
        end

        % Display them without label
        if size(STA,1)>1 & size(STA,2)>1,
            axis image
            set(get(ttH,'Parent'),'YTickLabel',[]);
            %set(get(ttH,'Parent'),'YTick',[]);
            set(get(ttH,'Parent'),'XTickLabel',[]);
            %set(get(ttH,'Parent'),'XTick',[]);
            %axis off
        end

        %set(h2,'YTickLabel',[]);
        %set(h2,'XTickLabel',[]);

        if (fr-handles.imageIndex * numList -1 ) ~= 0
            %curTitle = num2str((fr-1-(nt-1)/2)*binsize);
            curTitle = num2str((fr-1)*binsize);
        else
            curTitle = titlestr;
        end
        title(curTitle);
        colormap(redblue);
    end

elseif strcmp(sDim, '0-D')
    set(handles.strfonly, 'visible', 'on')
    axes(handles.strfonly);
    plot(ceil(t*1000/ampsamprate),forward);
    axis tight;
    xlabel('Time (ms)')
    title('STRF');
end

% Display appropriate tol-value on the window
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function plot_smoothedstrf(forwardJN, forwardJNstd, handles)
% --------------------------------------------------------------------

% Get calc parameter value from global variables
global sDim initialFreq endFreq ampsamprate preprocessOption


% Prepare for displaying axes
guidata(handles.figure1,handles);
% Get proper axes as current
set(handles.small_strf, 'Visible', 'off')
child = get(handles.small_strf, 'Children');
set(child, 'Visible', 'Off')
set(handles.small_strf_x, 'Visible', 'off')
child = get(handles.small_strf_x, 'Children');
set(child, 'Visible', 'Off')
set(handles.small_strf_y, 'Visible', 'Off')
child = get(handles.small_strf_y, 'Children');
set(child, 'Visible', 'Off')
set(handles.small_mod_strf, 'Visible', 'off')
child = get(handles.small_mod_strf, 'Children');
set(child, 'Visible', 'Off')
set(handles.STRFaxes, 'visible', 'off');
child = get(handles.STRFaxes, 'Children');
set(child, 'Visible', 'Off')
itsTitle = get(handles.STRFaxes, 'Title');
set(itsTitle, 'Visible', 'Off')
set(handles.STAaxes, 'visible', 'off');
child = get(handles.STAaxes, 'Children');
set(child, 'Visible', 'Off')
itsTitle = get(handles.STAaxes, 'Title');
set(itsTitle, 'Visible', 'Off')

% Set other axes to invisible
set(handles.strfonly, 'visible', 'off');
child = get(handles.strfonly, 'Children');
set(child, 'Visible', 'Off')
% Display appropriate tol-value on the window
global Tol_val
tol = Tol_val(handles.tolindex);

% total numbers of image on the screen
for ii = 1:handles.numList
    axesStr = strcat('subaxes', num2str(ii-1));
    subhandle = getfield(handles, axesStr);
    set(subhandle, 'Visible', 'off');
    child = get(subhandle, 'Children');
    set(child, 'Visible', 'off')

end
% 7/13/2005
% Smoothed Jackknifed version STRFs using Std_val
% Refer to filters4.m
global Std_val
nstd = length(Std_val);
[nx, nt, nJN]= size(forwardJN);
% Filter the filters
forwardJN_s = zeros(nx, nt, nJN, nstd);

for istd=1:nstd

    stdfilt = Std_val(istd);
    titlestr{istd} = sprintf('std = %3.2f\n',...
        Std_val(istd));

    for iJN =1:nJN
        %  find the filtered JN STRF
        forwardJN_s(:,:,iJN, istd) = fast_filter_filter(mean(forwardJN(:,:,[1:(iJN-1) (iJN+1):nJN]),3),...
            forwardJNstd(:,:,iJN), stdfilt);
    end
  
end

forward = squeeze(mean(forwardJN_s, 3));
maxforward = max(max(max(forward)));
minforward = min(min(min(forward)));
absforward = max(abs(minforward),abs(maxforward));

t = -(nt-1)/2:(nt-1)/2;
if isempty(initialFreq) | isempty(endFreq)
    f = 1:nx;
else
    f = linspace(initialFreq, endFreq, nx);
end


% Display estimated STRF based on sptial dimension
if strcmp(sDim, '1-D')

    % 7/13/2005
    % Need add smoothed STRFs display

    % Smoothed STRFs display
    nt = length(Std_val);  % how many smoothed STRFs for each Tol val
    numList = handles.numList;  % total numbers of STRFs can display in the window
    if numList > nt
        numList = nt;
    end

    % Smoothed STRFs display
    nt = length(Std_val);  % how many smoothed STRFs for each Tol val
    numList = handles.numList;  % total numbers of STRFs can display in the window
    if numList > nt
        numList = nt;
    end

    % display total 35 images on one screen
    for fr = handles.imageIndex * numList+1:numList+handles.imageIndex * numList

        if fr > nt
            break
        end

        % display them in appropriate axes
        axesStr = strcat('subaxes', num2str(fr-handles.imageIndex * numList-1));
        axes(getfield(handles, axesStr));

        % Display them in specific format
        if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
            flabel = logspace(log10(initialFreq), log10(endFreq), nx);
            ttH =pcolor(ceil(t*1000/ampsamprate), flabel/1000, forward(:,:,fr)); shading interp;
            caxis([-absforward absforward]);

        else
            ttH = imagesc(ceil(t*1000/ampsamprate),f/1000,forward(:,:,fr),[-absforward absforward]);
            axis xy;

        end
        %xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
        xlim([-20 80]);
        set(get(ttH,'Parent'),'YTickLabel',[]);
        %set(get(ttH,'Parent'),'YTick',[]);
        set(get(ttH,'Parent'),'XTickLabel',[]);

        title(titlestr(fr));

    end

elseif strcmp(sDim, '2-D') % 2-D display
    global stimsamprate
    if isempty(stimsamprate)
        binsize = 1;
        %titlestr = 'Latency(frame)=';
    else
        binsize = ceil((1/stimsamprate)*1000);
        %titlestr = 'Latency(msec)=';
    end
    global preprocessOption
    if strcmp(preprocessOption, 'Fourier Power transform')
        phasecount=1;
        tbincount = nt;
        kcount = 1;
        spacebincount = nx;
        chancount=spacebincount/phasecount;
        Xmax=sqrt(chancount*2);
        [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);

        tsf=zeros(Xmax*Xmax,tbincount,kcount);
        tsf(cfilt,:,:)=squeeze(sum(reshape(forward(:,:,handles.stdindex),chancount,phasecount,...
            tbincount,kcount),2));

        tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
        tsf=reshape(tsf,Xmax,Xmax,tbincount);
        ii = 1;
        while ii+4 < nstd
            figure('name', ['All smoothed strfs for tol val = ', num2str(tol)]);
            showkern(tsf(:,:,ii:ii+4),'space',[splitX splitX],titlestr(ii:ii+4), 0, binsize);
            ii = ii+5;
        end
        figure('name', ['All smoothed strfs for tol val = ', num2str(tol)]);
        showkern(tsf(:,:,ii:nstd),'space',[splitX splitX],titlestr(ii:nstd), 0, binsize);
        
    else
        % 2-D display
        splitX = floor(sqrt(nx));
        %tsf = reshape(forward, splitX, splitX, nt);
        ii = 1;
        while ii+4 < nstd
            figure('name', ['All smoothed strfs for tol val = ', num2str(tol)]);
            showkern(forward(:,:,ii:ii+4),'space',[splitX splitX],titlestr(ii:ii+4), 0, binsize);
            ii = ii+5;
        end
        figure('name', ['All smoothed strfs for tol val = ', num2str(tol)]);
        showkern(forward(:,:,ii:nstd),'space',[splitX splitX],titlestr(ii:nstd), 0, binsize);
    end

elseif strcmp(sDim, '0-D')
      % Smoothed STRFs display
    nt = length(Std_val);  % how many smoothed STRFs for each Tol val
    numList = handles.numList;  % total numbers of STRFs can display in the window
    if numList > nt
        numList = nt;
    end

    % Smoothed STRFs display
    nt = length(Std_val);  % how many smoothed STRFs for each Tol val
    numList = handles.numList;  % total numbers of STRFs can display in the window
    if numList > nt
        numList = nt;
    end

    % display total 35 images on one screen
    for fr = handles.imageIndex * numList+1:numList+handles.imageIndex * numList

        if fr > nt
            break
        end

        % display them in appropriate axes
        axesStr = strcat('subaxes', num2str(fr-handles.imageIndex * numList-1));
        axes(getfield(handles, axesStr));

        % Display them in specific format   
        ttH=plot(ceil(t*1000/ampsamprate),forward(:,fr));
        title(titlestr(fr));
        set(get(ttH,'Parent'),'YTickLabel',[]);
        %set(get(ttH,'Parent'),'YTick',[]);
        set(get(ttH,'Parent'),'XTickLabel',[]);

    end
end

guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function plot_strfProjection(forward,forwardJN_std, handles)
% --------------------------------------------------------------------

global sDim ampsamprate

maxforward = max(max(forward));
minforward = min(min(forward));
absforward = max(abs(minforward),abs(maxforward));

[nx, nt] = size(forward);
handles.nt = nt;
guidata(handles.figure1, handles);

t = -(nt-1)/2:(nt-1)/2;
f = 1:nx;

global initialFreq endFreq preprocessOption
if ~isempty(initialFreq) & ~isempty(endFreq)
    fstep = (endFreq - initialFreq)/nx;
    faxis = initialFreq+fstep/2:fstep:endFreq;
else
    faxis = f*1000;
end

% Set other axes to invisible
for ii = 1:handles.numList
    axesStr = strcat('subaxes', num2str(ii-1));
    subhandle = getfield(handles, axesStr);
    set(subhandle, 'Visible', 'off');
    child = get(subhandle, 'Children');
    set(child, 'Visible', 'Off');
    itsTitle = get(subhandle, 'Title');
    set(itsTitle, 'Visible', 'Off')

end
set(handles.strfonly, 'Visible', 'off');
child = get(handles.strfonly, 'Children');
set(child, 'Visible', 'Off')
itsTitle = get(handles.strfonly, 'Title');
set(itsTitle, 'Visible', 'Off')
set(handles.STRFaxes, 'visible', 'off');
child = get(handles.STRFaxes, 'Children');
set(child, 'Visible', 'Off')
itsTitle = get(handles.STRFaxes, 'Title');
set(itsTitle, 'Visible', 'Off')
set(handles.STAaxes, 'visible', 'off');
child = get(handles.STAaxes, 'Children');
set(child, 'Visible', 'Off')
itsTitle = get(handles.STAaxes, 'Title');
set(itsTitle, 'Visible', 'Off')

% Set other axes as current handles
set(handles.small_strf, 'Visible', 'On')
set(handles.small_strf_x, 'Visible', 'On')
set(handles.small_strf_y, 'Visible', 'On')
set(handles.small_mod_strf, 'Visible', 'On')

% Display STRF
% subplot(2,2,1);

if  strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
    axes(handles.small_strf)
    flabel = logspace(log10(initialFreq), log10(endFreq), nx);
    pcolor(ceil(t*1000/ampsamprate), flabel/1000, forward); shading interp;
    caxis([-absforward absforward]);
else
    axes(handles.small_strf)
    imagesc(ceil(t*1000/ampsamprate),faxis/1000,forward,[-absforward absforward]);
    axis xy;
end
if strcmp(sDim, '1-D')
    %axis([-20 100 initialFreq/1000 endFreq/1000]);
    xlim([-20 100])
    %xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
    xlabel('Time (ms)')
    ylabel('Frequency (kHz)');

end
title('STRF');

if ( abs(maxforward) > abs(minforward) )
    [fpeak tpeak] = find(forward==maxforward);
else
    [fpeak tpeak] = find(forward==minforward);
end

% Display projection on time domain
% subplot(2,2,2);
axes(handles.small_strf_x)

if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
    %'SongWave->Spectrogram(STFFT)')

    % Confidence Interval centered at mean
    
    plot(forward(:,tpeak)',flabel'/1000, forward(:,tpeak)'+forwardJN_std(:,tpeak)'*2.0,...
        flabel'/1000, forward(:,tpeak)'+forwardJN_std(:,tpeak)'*-2.0,flabel'/1000);
else
    plot(forward(:,tpeak)',faxis'/1000, forward(:,tpeak)'+forwardJN_std(:,tpeak)'*2.0,...
        faxis'/1000, forward(:,tpeak)'+forwardJN_std(:,tpeak)'*-2.0,faxis'/1000);
end

if ~isempty(initialFreq) & ~isempty(endFreq)
    ylim([initialFreq/1000 endFreq/1000])
    xlabel('Gain');
    ylabel('Frequency (kHz)');
end
titlestring=sprintf('   At Time =%6.2f (ms)', t(tpeak));
title(titlestring);

if 0
    % fit with a guassian to get bandwidth
    beta(1)= forward(fpeak,tpeak);
    beta(2)= faxis(fpeak);
    beta(3)= 500.0;
    options = foptions;
    options(9)=1;
    [new_beta, options, res, jac] = curvefit('gaussfit',beta,...
        f,forward(:,tpeak)',options,'gaussgrad');
    forward_fpeak = gaussfit(new_beta,f);
    hold on
    plot(faxis/1000,forward_fpeak,'k');
    titlestring=sprintf('   CF=%6.2f  BW=%6.2f', faxis(fpeak), new_beta(3));

    %hold off
end

%subplot(2,2,3);
axes(handles.small_strf_y)
forwardstd2 = squeeze(mean(forwardJN_std, 3));
plot(ceil(t*1000/ampsamprate),forward(fpeak,:),...
    ceil(t*1000/ampsamprate),forward(fpeak,:)+2*forwardstd2(fpeak,:),...
    ceil(t*1000/ampsamprate),forward(fpeak,:)-2*forwardstd2(fpeak,:));
axis([-20 100 minforward maxforward]);

if strcmp(preprocessOption, 'SongWave->Spectrogram(STFFT)')
    titlestring=sprintf(' At Freq (Hz)=%6.2f', faxis(fpeak));
elseif strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
    titlestring=sprintf(' At Freq (Hz)=%6.2f', flabel(fpeak));
else
    titlestring=sprintf(' At Freq (Hz)=%6.2f', faxis(fpeak));
end

title(titlestring);
%xlabel(sprintf('Time ( in %2.2f ms)', 1000/ampsamprate))
xlabel('Time (ms)')
ylabel('Gain');

%subplot(2,2,4)
axes(handles.small_mod_strf)
[pow fpow]=psd(forward(fpeak,:), length(forward(fpeak,:)), 1000);
plot(fpow, pow);
powmax=max(pow);
bmf_index=find(pow==powmax);
bmf=fpow(bmf_index);

titlestring=sprintf('BMF=%6.2f', bmf);
title(titlestring);
xlabel('Modulation Frequency Hz');
ylabel('Power');


% --------------------------------------------------------------------
function plot_strfSTA(forward, handles)
% --------------------------------------------------------------------

% Axes preparation for displaying
% Get proper axes as current
set(handles.small_strf, 'Visible', 'off')
child = get(handles.small_strf, 'Children');
set(child, 'Visible', 'Off')
set(handles.small_strf_x, 'Visible', 'off')
child = get(handles.small_strf_x, 'Children');
set(child, 'Visible', 'Off')
set(handles.small_strf_y, 'Visible', 'Off')
child = get(handles.small_strf_y, 'Children');
set(child, 'Visible', 'Off')
set(handles.small_mod_strf, 'Visible', 'off')
child = get(handles.small_mod_strf, 'Children');
set(child, 'Visible', 'Off')
set(handles.strfonly, 'Visible', 'off');
child = get(handles.strfonly, 'Children');
set(child, 'Visible', 'Off')
% variable numList defined as total num of images on one screen
for ii = 1:handles.numList
    axesStr = strcat('subaxes', num2str(ii-1));
    subhandle = getfield(handles, axesStr);
    set(subhandle, 'Visible', 'off');
    child = get(subhandle, 'Children');
    set(child, 'Visible', 'off');
    titlethis = get(subhandle, 'Title');
    set(titlethis, 'Visible', 'off');

end

set(handles.STRFaxes, 'visible', 'on');
set(handles.STAaxes, 'visible', 'on');

% Get STA values
global outputPath
crossfile = fullfile(outputPath,'StimResp_crosscorr.mat');
stim_spike = Check_And_Load(crossfile);

[nb, nt] = size(stim_spike);
stim_spike = fliplr(stim_spike);
w = hanning(nt);
for ib=1:nb
    stim_spike(ib,:) = stim_spike(ib,:).*w';
end

handles.nt = nt;
guidata(handles.figure1, handles);

halft = floor(nt/2);
t = -halft:halft;
global sDim ampsamprate

if strcmp(sDim, '1-D')
    global initialFreq endFreq preprocessOption
    if isempty(initialFreq) | isempty(endFreq)
        f = 1:nb;
        label_yaxis = 'Frequency Band';
    else
        f = initialFreq:endFreq;
        label_yaxis = 'Frequency (Hz)';
    end

    %STRF
    maxforward = max(max(forward));
    minforward = min(min(forward));
    absforward = max(abs(minforward),abs(maxforward));

    if strcmp(preprocessOption, 'SongWave->Scalogram(wavelet)')
        axes(handles.STRFaxes);
        flabel = logspace(log10(initialFreq), log10(endFreq), nb);
        pcolor(ceil(t*1000/ampsamprate), flabel, forward); shading interp;
        caxis([-absforward absforward]);
        ylabel(label_yaxis)
        title('STRF');
        axes(handles.STAaxes);
        pcolor(ceil(t*1000/ampsamprate), flabel, stim_spike); shading interp;
        ylabel(label_yaxis)

        xlabel('Time (ms)')
        title('STA')
    else
        axes(handles.STRFaxes);
        imagesc(ceil(t*1000/ampsamprate),f,forward, [-absforward absforward]);
        axis xy;
        ylabel(label_yaxis)
        title('STRF');
        %STA
        axes(handles.STAaxes);
        imagesc(ceil(t*1000/ampsamprate),f,stim_spike);
        axis xy;
        ylabel(label_yaxis)

        xlabel('Time (ms)')
        title('STA')
    end



elseif strcmp(sDim, '2-D')
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
        tsf(cfilt,:,:)=squeeze(sum(reshape(forward,chancount,phasecount,...
            tbincount,kcount),2));

        tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
        STRF=reshape(tsf,Xmax,Xmax,tbincount);
        clear tsf
        tsf=zeros(Xmax*Xmax,tbincount,kcount);
        tsf(cfilt,:,:)=squeeze(sum(reshape(stim_spike,chancount,phasecount,...
            tbincount,kcount),2));

        tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
        STA=reshape(tsf,Xmax,Xmax,tbincount);
    else
        % 2-D display
        splitX = floor(sqrt(nb));
        STA = reshape(stim_spike, splitX, splitX, nt);
        STRF = reshape(forward, splitX, splitX, nt);
        
    end
%     splitX = floor(sqrt(nb));
%     STA = reshape(stim_spike, splitX, splitX, nt);
%     STRF = reshape(forward, splitX, splitX, nt);
    amaxSTA = max(abs(STA(:)));
    aminSTA = -amaxSTA;
    amaxSTRF = max(abs(STRF(:)));
    aminSTRF = -amaxSTRF;

    % For general case, we set numList as handles.numList
    % But for better display, we only display 10 frames
    % Need work on this in the future. -- JXZ, 5/30/03
    if 0
        numList = handles.numList;
    else
        numList = 10+(nt-1)/2;
    end

    if numList > nt
        numList = nt;
    end

    wholeSTA = [];
    wholeSTRF = [];

    % Only display video frame from time = 0
    for i = 1+(nt-1)/2:numList
        wholeSTA = [wholeSTA STA(:,:,i)];
        wholeSTRF = [wholeSTRF STRF(:,:,i)];
    end
    axes(handles.STRFaxes);
    ttH = imagesc(wholeSTRF);

    %ylabel('Latency (ms)');
    title('STRF');
    colormap(redblue);
    axis image
    set(get(ttH,'Parent'),'YTickLabel',[]);
    set(get(ttH,'Parent'),'YTick',[]);
    set(get(ttH,'Parent'),'XTickLabel',[]);
    set(get(ttH,'Parent'),'XTick',[]);

    axes(handles.STAaxes);
    ttH = imagesc(wholeSTA);

    %ylabel('Latency (ms)');
    title('STA');
    %colormap(redblue);
    axis image
    set(get(ttH,'Parent'),'YTickLabel',[]);
    set(get(ttH,'Parent'),'YTick',[]);
    set(get(ttH,'Parent'),'XTickLabel',[]);
    set(get(ttH,'Parent'),'XTick',[]);

elseif strcmp(sDim, '0-D')
    set(handles.STRFaxes, 'visible', 'on');
    axes(handles.STRFaxes);
    plot(ceil(t*1000/ampsamprate),forward)
    axis tight;
    title('STRF')
    set(handles.STAaxes, 'visible', 'on');
    axes(handles.STAaxes);
    plot(ceil(t*1000/ampsamprate), stim_spike)
    axis tight;
    xlabel('Time (ms)')
    title('STA')
end

% --------------------------------------------------------------------
function varargout = nextlag_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
if isfield(handles,'strf_result') & ~isempty(handles.strf_result)
    global outputPath

    strf_result = load(fullfile(outputPath,...
        ['strfResult_Tol',num2str(handles.Index),'.mat']));

    handles.imageIndex = handles.imageIndex +1;
    if handles.imageIndex*handles.numList> handles.nt
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

    % Only display STRF
    v = get(handles.displayoption, 'value');
    switch v
        case 1 % STRF only option
            plot_strfonly(handles);
    end
end

% --------------------------------------------------------------------
function varargout = prevlag_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
if isfield(handles,'strf_result') & ~isempty(handles.strf_result)
    global outputPath

    strf_result = load(fullfile(outputPath,...
        ['strfResult_Tol',num2str(handles.Index),'.mat']));
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

    % Only display STRF
    v = get(handles.displayoption, 'value');
    switch v
        case 1 % STRF only option
            plot_strfonly(handles);
    end
end

% --------------------------------------------------------------------
function varargout = displaystrf_help_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
ttlStr = get(handles.figure1, 'Name');
hlpStr = [...
    '                                                               '
    '    STRFPAK: DISPLAY STRF Window                               '
    '                                                               '
    ' This window visually displays the estimation results. There   '
    ' are five display options. The information panel shows the     '
    ' associated parameters of the STRF, e.g. Tol Val and sparseness'
    ' parameters and etc.                                           '
    '                                                               '
    ' Five display options include:                                 '
    '    1. STRF only: For 1-D data, STRF is image plot as time vs  '
    '        frequency; for 2-D data, STRF is frames verse latency. '
    '        The total frames on the window is up to 35 frames.     '
    '    2. Filtered STRFs for each Tol value: Here use the second  '
    '      regularization parameter: sparseness to filter the STRF. '
    '    3. All filtered STRFs (display in separate window): For 1-D'
    '       case, the filtered STRFs are displayed in a separate    '
    '       window. For 2-D case, there are number of figures       '
    '       showing the filtered STRFs.                             '
    '    4. STRF with projection on two axes (for 1-D only): On left'
    '       panel, there are four plots: unfiltered STRF,  the time '
    '       slice with its confidence interval, the frequency slice '
    '       with its confidence interval and the best modulation    '
    '       frequency.                                              '
    '    5. STRF and STA: Display strf and STA together to get some '
    '        idea how different they look like visually.            '
    '                                                               '
    ' Information panel includes:                                   '
    '    Output Dir: the output directory is shown up.              '
    '    Tol Val box: shows the tolerance values used for the STRF. '
    '    Tol val slider: use mouse can choose next or prev tol val  '
    '                  and the corresponding STRF.                  '
    '    Spareness box: shows the sparseness values used for the    '
    '                  STRF.                                        '
    '    Spareness slider: use mouse can choose next or prev        '
    '               sparseness value and the corresponding STRF.    '
    '    Help: help window shows up.                                '
    '    Close: close this window if it is clicked.                 '
    '                                                               '
    'Updated by Junli, May 2006.                                    '
    '                                                               '];
myFig = handles.figure1;
helpwin(hlpStr, ttlStr);

% --------------------------------------------------------------------
% END of displaystrf_GUI
% --------------------------------------------------------------------


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
guidata(handles.figure1, handles);
displayall(handles);

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
guidata(handles.figure1, handles);
displayall(handles);

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
displayall(handles)

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
    handles.stdindex = newIndex;
    %set(hObject, 'String', num2str(handles.Index));
    % Update slider value based on dataSet value
    set(handles.stdvalSlider, 'String', newIndex);
end
guidata(handles.figure1, handles);
displayall(handles);

% --- Executes during object creation, after setting all properties.
function stdvalshow_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in outdirBrowser.
function outdirBrowser_Callback(hObject, eventdata, handles)
% hObject    handle to outdirBrowser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function outdirshow_Callback(hObject, eventdata, handles)
% hObject    handle to outdirshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outdirshow as text
%        str2double(get(hObject,'String')) returns contents of outdirshow as a double


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


