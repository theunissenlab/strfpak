global rawDS
if isempty(rawDS)
    set(handles.displayrawstim_new,'Enable','off');
    set(handles.displayrawstimpsth_new,'Enable','off');
    set(handles.preprocess_new,'Enable','off');
    set(handles.Special_options,'Enable','off');
    set(handles.Special_options_warning,'Visible','off');
else
    set(handles.Special_options,'Enable','on');
    set(handles.displayrawstim_new,'Enable','on');
    set(handles.displayrawstimpsth_new,'Enable','on');
    global rawData allow_negative_rates use_alien_space
    if rawData
        set(handles.preprocess_new,'Enable','on');
    else
        set(handles.preprocess_new,'Enable','off');
        set(handles.displayrawstim_new,'Enable','off');
        set(handles.displayrawstimpsth_new,'Enable','off');
    end
    if isempty(allow_negative_rates)
        allow_negative_rates = 0;
    end
    if isempty(use_alien_space)
        use_alien_space = 0;
    end
    if allow_negative_rates | use_alien_space
        set(handles.Special_options_warning,'Visible','on');
    else
        set(handles.Special_options_warning,'Visible','off');
    end
end

global DS
if isempty(DS)
    set(handles.ppdata_display,'Enable','off');
    set(handles.strf_calculate,'Enable','off');
else
    set(handles.ppdata_display,'Enable','on');
    set(handles.strf_calculate,'Enable','on');
end
global outputPath
if ~isempty(outputPath)
    if exist(fullfile(outputPath,'predResult_avgSpike1.mat'),'file') | exist(fullfile(outputPath,'predResult_avgSpike2.mat'),'file')
        handles.predictionOK = 1;
    end
    if exist(fullfile(outputPath,'best_strf.mat'),'file')
        handles.validationOK = 1;
    end
end

if ~isfield(handles,'predictionOK')
    handles.predictionOK = 0;
end
if ~handles.predictionOK
    set(handles.strf_displaystimstat,'Enable','off');
    set(handles.strf_displaystrf,'Enable','off');
    set(handles.strf_displaypredstrf,'Enable','off');
    set(handles.goodnessfitting,'Enable','off');
else
    set(handles.strf_displaystimstat,'Enable','on');
    set(handles.strf_displaystrf,'Enable','on');
    set(handles.strf_displaypredstrf,'Enable','on');
    set(handles.goodnessfitting,'Enable','on');
end

if ~isfield(handles,'validationOK')
    handles.validationOK = 0;
end

if ~handles.validationOK
    set(handles.strf_displaybestfit,'Enable','off');
else
    set(handles.strf_displaybestfit,'Enable','on');
end