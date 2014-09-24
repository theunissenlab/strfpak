function varargout = select_validation_files(varargin)
% SELECT_VALIDATION_FILES M-file for select_validation_files.fig
%      SELECT_VALIDATION_FILES, by itself, creates a new SELECT_VALIDATION_FILES or raises the existing
%      singleton*.
%
%      H = SELECT_VALIDATION_FILES returns the handle to a new SELECT_VALIDATION_FILES or the handle to
%      the existing singleton*.
%
%      SELECT_VALIDATION_FILES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_VALIDATION_FILES.M with the given input arguments.
%
%      SELECT_VALIDATION_FILES('Property','Value',...) creates a new SELECT_VALIDATION_FILES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_validation_files_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_validation_files_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_validation_files

% Last Modified by GUIDE v2.5 14-Sep-2006 13:11:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @select_validation_files_OpeningFcn, ...
    'gui_OutputFcn',  @select_validation_files_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
if nargin == 0
    fig = openfig(mfilename,'reuse');


    hAxes = findall(fig,'type','axes');
    hText  = findall(hAxes,'type','text');
    hUIControls = findall(fig,'type','uicontrol');
    set([hAxes; hText;...
        hUIControls],'fontname', 'Times New Roman','units','normalized','fontunits','normalized');


    global DS predDS rawDS for_validation for_final finalDS
    %  This next part is becasue DS, etc, are struct cells, not struct
    [DS_names, predDS_names,finalDS_names,for_validation,for_final,DS,predDS,finalDS] = process_DSs(DS,predDS,finalDS,rawDS);
%     if length(DS) > 0
%         [DS_names,order] = sort(DS_names);
%         DS = {DS{order}};
%     end
%     if length(predDS) > 0
%         [predDS_names,order] = sort(predDS_names);
%         predDS = {predDS{order}};
%     end
%     if length(finalDS) > 0
%         [finalDS_names,order] = sort(finalDS_names);
%         finalDS = {finalDS{order}};
%     end
sort_DSs;

    handles = guihandles(fig);
    finalDS_names = beautify_filenames(finalDS);
    predDS_names = beautify_filenames(predDS);
    DS_names = beautify_filenames(DS);
    set(handles.Calculation_data,'String',DS_names);
    set(handles.Validation_data,'String',predDS_names);
    set(handles.Final_Validation_data,'String',finalDS_names);
    set(handles.Validation_data,'Value',[]);
    set(handles.Calculation_data,'Value',[]);
    set(handles.Final_Validation_data,'Value',[]);


end


% --- Executes just before select_validation_files is made visible.
function select_validation_files_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_validation_files (see VARARGIN)

% Choose default command line output for select_validation_files
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes select_validation_files wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = select_validation_files_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Calculation_data.
function Calculation_data_Callback(hObject, eventdata, handles)
% hObject    handle to Calculation_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Calculation_data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Calculation_data


% --- Executes during object creation, after setting all properties.
function Calculation_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Calculation_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Validation_data.
function Validation_data_Callback(hObject, eventdata, handles)
% hObject    handle to Validation_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Validation_data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Validation_data


% --- Executes during object creation, after setting all properties.
function Validation_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Validation_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Add_selected_files.
function Add_selected_files_Callback(hObject, eventdata, handles)
% hObject    handle to Add_selected_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DS predDS rawDS for_validation

index_selected = get(handles.Calculation_data,'Value');
to_remain = ones(1,length(DS));
to_remain(index_selected) = 0;
index_remain = find(to_remain);
if length(index_selected) > 0
    if length(predDS) > 0
        predDS = {predDS{:} DS{index_selected}};
    else
        predDS = {DS{index_selected}};
    end
    DS = {DS{index_remain}};
    set(handles.Calculation_data,'Value',[])

    select_validation_files;
end

% --- Executes on button press in Remove_selected_files.
function Remove_selected_files_Callback(hObject, eventdata, handles)
% hObject    handle to Remove_selected_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DS predDS rawDS for_validation

index_selected = get(handles.Validation_data,'Value');
to_remain = ones(1,length(predDS));
to_remain(index_selected) = 0;
index_remain = find(to_remain);
if length(index_selected) > 0

    if length(DS) > 0
        DS = {DS{:} predDS{index_selected}};
    else
        DS = {predDS{index_selected}};
    end
    predDS = {predDS{index_remain}};
    set(handles.Validation_data,'Value',[])


    select_validation_files;
end

% --- Executes on button press in get_help.
function get_help_Callback(hObject, eventdata, handles)
% hObject    handle to get_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
message = ['You are now being asked to select:' char(10) ...
    '1: which files you want to use for STRF calculation' char(10) ...
    '2: which files you will reserve for determining the best STRF, and' char(10) ...
    '3: which ones are left for validation only.' char(10) ...
    'For best results, leave a minimum of three files ' char(10) ...
    'in each category.'];
msgbox(message,'Validation Help');
% --- Executes on button press in Done_data_selection.
function Done_data_selection_Callback(hObject, eventdata, handles)
% hObject    handle to Done_data_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
add_to_STRFPAK_script('select_validation_files.m','select_validation');
global outputPath for_validation for_final
save(fullfile(outputPath,'STRFPAK_script_dataset.mat'),'for_validation','for_final','-APPEND');


delete(handles.figure1);

function sort_DSs(nothing)


%%%$$$ begin select_validation
%  This code should be  used exclusively in script mode; it's just to help
%  select data out of the "rawDS" which is found automatically.
global DS predDS rawDS for_validation for_final finalDS originalDS

DS = originalDS;
predDS = {};
finalDS = {};

%if isempty(predDS) | isempty(finalDS)

    %  If the for_validation vector generated by the GUI run is too long,
    %  it gets cut, and if it's too short, it gets repeated.
    if length(for_validation) > length(DS)
        for_validation = for_validation(1:length(DS));
        for_final = for_final(1:length(DS));
    end
    if length(for_validation) < length(DS)
        for_validation = for_validation(1+mod(0:(length(DS)-1),length(for_validation)));
        for_final = for_final(1+mod(0:(length(DS)-1),length(for_final)));
    end

    index_selected = find(for_validation);
    to_remain = ones(1,length(DS));
    to_remain(index_selected) = 0;
    index_remain = find(to_remain);
    if length(index_selected) > 0
        if length(predDS) > 0
            predDS = {predDS{:} DS{index_selected}};
        else
            predDS = {DS{index_selected}};
        end        
    end
        index_selected_final = find(for_final);
    to_remain_final = ones(1,length(DS));
    to_remain_final(index_selected_final) = 0;
    index_remain_final = find(to_remain_final);
    if length(index_selected_final) > 0
        if length(finalDS) > 0
            finalDS = {finalDS{:} DS{index_selected_final}};
        else
            finalDS = {DS{index_selected_final}};
        end     
    end
    
    
    
    DS = {DS{intersect(index_remain,index_remain_final)}};

%end
global predinitialFreq predendFreq predampsamprate sDim;
global  initialFreq endFreq ampsamprate;

predinitialFreq = initialFreq;
predendFreq = endFreq;
predampsamprate = ampsamprate;


%%%$$$ end select_validation


% --- Executes on selection change in Final_Validation_data.
function Final_Validation_data_Callback(hObject, eventdata, handles)
% hObject    handle to Final_Validation_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Final_Validation_data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Final_Validation_data


% --- Executes during object creation, after setting all properties.
function Final_Validation_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Final_Validation_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in add_selected_final.
function add_selected_final_Callback(hObject, eventdata, handles)
% hObject    handle to add_selected_final (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global finalDS predDS rawDS

index_selected = get(handles.Validation_data,'Value');
to_remain = ones(1,length(predDS));
to_remain(index_selected) = 0;
index_remain = find(to_remain);
if length(index_selected) > 0
    if length(finalDS) > 0
        finalDS = {finalDS{:} predDS{index_selected}};
    else
        finalDS = {predDS{index_selected}};
    end
    predDS = {predDS{index_remain}};
    set(handles.Validation_data,'Value',[])

    select_validation_files;
end


% --- Executes on button press in remove_selected_final.
function remove_selected_final_Callback(hObject, eventdata, handles)
% hObject    handle to remove_selected_final (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global finalDS predDS rawDS

index_selected = get(handles.Final_Validation_data,'Value');
to_remain = ones(1,length(finalDS));
to_remain(index_selected) = 0;
index_remain = find(to_remain);
if length(index_selected) > 0

    if length(predDS) > 0
        predDS = {predDS{:} finalDS{index_selected}};
    else
        predDS = {finalDS{index_selected}};
    end
    finalDS = {finalDS{index_remain}};
    set(handles.Final_Validation_data,'Value',[])


    select_validation_files;
end
