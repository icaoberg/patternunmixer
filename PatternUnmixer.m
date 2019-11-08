function varargout = PatternUnmixer(varargin)
% PATTERNUNMIXER M-file for PatternUnmixer.fig
%      PATTERNUNMIXER, by itself, creates a new PATTERNUNMIXER or raises the existing
%      singleton*.
%
%      H = PATTERNUNMIXER returns the handle to a new PATTERNUNMIXER or the handle to
%      the existing singleton*.
%
%      PATTERNUNMIXER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PATTERNUNMIXER.M with the given input arguments.
%
%      PATTERNUNMIXER('Property','Value',...) creates a new PATTERNUNMIXER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PatternUnmixer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PatternUnmixer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PatternUnmixer

% Last Modified by GUIDE v2.5 11-Jul-2010 10:59:33

% Authors: Tao Peng,  Abdul-Saboor Sheikh and Robert F. Murphy
% Past Contributors:Ting Zhao, Ivan E. Cao-Berg
% Created: April 6, 2009
%
% Copyright (C) 2009 Center for Bioimage Informatics/Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PatternUnmixer_OpeningFcn, ...
                   'gui_OutputFcn',  @PatternUnmixer_OutputFcn, ...
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

addpath(genpath('./functions'))
addpath(genpath('./functions/omero'))
% End initialization code - DO NOT EDIT


% --- Executes just before PatternUnmixer is made visible.
function PatternUnmixer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PatternUnmixer (see VARARGIN)
handles.PUversion = '2.0';
handles.traindata = [];
handles.traindatachs = [];
handles.lastdir = pwd;
handles.defaultdir = pwd;
handles.models = [];
handles = setResultTableRowNames(handles);
set(handles.modellist,'String','');
handles.modelFileList = [];
datapath = strcat(handles.defaultdir,filesep,'data',filesep,'train');
if exist(datapath,'dir')
    handles.defaultdir = datapath;
end
    storedModels = dir(strcat(handles.defaultdir,filesep,'model*.mat'));
    modelFiles = [];
    for modelIdx = 1:length(storedModels)
        modelFiles{modelIdx} = strcat(handles.defaultdir,filesep,storedModels(modelIdx).name);
    end
    handles = loadModels(modelFiles,handles);


if(exist('ptrn.mat','file'))
    s = load('ptrn');
    if(~isempty(s.preDefPatternList))
        set(handles.predefPatternList,'String',s.preDefPatternList);
        index = get(handles.predefPatternList,'Value');
        if (index(end)>0)
            ptrns = get(handles.predefPatternList,'String');
            set(handles.newPatternName,'String',ptrns(index(end)));
            enablePatternListControls(handles, 'on');
        end
    end
end
% Choose default command line output for PatternUnmixer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PatternUnmixer wait for user response (see UIRESUME)
% uiwait(handles.addedPatternList);

function handles = loadModels(modelFiles,handles)
modelslist = get(handles.modellist,'String');
for modelIdx = 1:length(modelFiles)
    curFileName = char(modelFiles{modelIdx});
    try
        load(curFileName);
        if(isstruct(model))
            if isempty(handles.models)
                handles.models = model;
            else
                handles.models(end+1) = model;
            end
            model = [];
            handles.modelFileList{end+1} = curFileName;
            modelslist{end+1} = handles.models(end).name;
        end
    catch exception
        disp(strcat(curFileName,' is not an unmixer model'));
    end
end
set(handles.modellist,'String',modelslist);
set(handles.modellist,'Value',length(modelslist));
%set(handles.models,'Struct',handles.models);
%guidata(hObject, handles);
setTableColumns(handles);

% --- Outputs from this function are returned to the command line.
function varargout = PatternUnmixer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in protimgs.
function protimgs_Callback(hObject, eventdata, handles)
% hObject    handle to protimgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns protimgs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protimgs


% --- Executes during object creation, after setting all properties.
function protimgs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protimgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in add1.
function add1_Callback(hObject, eventdata, handles)
% hObject    handle to add1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (isInvalidCh(handles.protCh))
    return
end

filelist = get(handles.protimgs,'String');
[filelistadd,fpath] = getimgurls('Select protein images',handles.lastdir);
if fpath~=0
    handles.lastdir = fpath;
end
filelist = [filelist;filelistadd];
filelist = unique(filelist);
set(handles.protimgs,'Value',length(filelist));
set(handles.protimgs,'Max',length(filelist));
set(handles.protimgs,'String',filelist);
guidata(hObject, handles);

function invalidCh = isInvalidCh(hObject)
invalidCh = 0;
ch = uint8(str2double(get(hObject,'String')));
if (ch < 1)
    uicontrol(hObject);
    msgbox('Channel to be imported must be a non-zero integer','Invalid Input','error');
    invalidCh = 1;
end

% --- Executes on button press in omeroAddProt.
function omeroAddProt_Callback(hObject, eventdata, handles)
% hObject    handle to omeroAddProt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isInvalidCh(handles.protCh))
    return
end

filelist = get(handles.protimgs,'String');
[filelistadd,fpath] = omeroImageSelectionDialog(uint8(get(handles.protCh,'Value')));
if fpath~=0
    handles.lastdir = fpath;
end
filelist = [filelist;filelistadd];
filelist = unique(filelist)';
set(handles.protimgs,'Value',length(filelist));
set(handles.protimgs,'Max',length(filelist));
set(handles.protimgs,'String',filelist);
guidata(hObject, handles);



% --- Executes on button press in delete1.
function delete1_Callback(hObject, eventdata, handles)
% hObject    handle to delete1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filelist = get(handles.protimgs,'String');
delidx = get(handles.protimgs,'Value');
if delidx
    filelist(delidx) = [];
    set(handles.protimgs,'Value',length(filelist));
    set(handles.protimgs,'Max',length(filelist));
    set(handles.protimgs,'String',filelist);
end
guidata(hObject, handles);

% --- Executes on button press in clear1.
function clear1_Callback(hObject, eventdata, handles)
% hObject    handle to clear1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.protimgs,'Value',0);
set(handles.protimgs,'Max',0);
set(handles.protimgs,'String','');
guidata(hObject, handles);


% --- Executes on selection change in dnaimgs.
function dnaimgs_Callback(hObject, eventdata, handles)
% hObject    handle to dnaimgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns dnaimgs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dnaimgs


% --- Executes during object creation, after setting all properties.
function dnaimgs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dnaimgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in add2.
function add2_Callback(hObject, eventdata, handles)
% hObject    handle to add2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (isInvalidCh(handles.dnaCh))
    return
end

filelist = get(handles.dnaimgs,'String');
[filelistadd,fpath] = getimgurls('Select DNA images',handles.lastdir);
if fpath~=0
    handles.lastdir = fpath;
end
filelist = [filelist;filelistadd];
filelist = unique(filelist);
set(handles.dnaimgs,'Value',length(filelist));
set(handles.dnaimgs,'Max',length(filelist));
set(handles.dnaimgs,'String',filelist);
guidata(hObject, handles);


% --- Executes on button press in delete2.
function delete2_Callback(hObject, eventdata, handles)
% hObject    handle to delete2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filelist = get(handles.dnaimgs,'String');
delidx = get(handles.dnaimgs,'Value');
if delidx
    filelist(delidx) = [];
    set(handles.dnaimgs,'Value',length(filelist));
    set(handles.dnaimgs,'Max',length(filelist));
    set(handles.dnaimgs,'String',filelist);
end
guidata(hObject, handles);


% --- Executes on button press in clear2.
function clear2_Callback(hObject, eventdata, handles)
% hObject    handle to clear2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.dnaimgs,'Value',0);
set(handles.dnaimgs,'Max',0);
set(handles.dnaimgs,'String','');
guidata(hObject, handles);

% --- Executes on button press in addpattern.
function addpattern_Callback(hObject, eventdata, handles)
% hObject    handle to addpattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = addPatternToModel(hObject,handles);
guidata(hObject, handles);

function handles = addPatternToModel(hObject,handles)
protlist = get(handles.protimgs,'String');
dnalist = get(handles.dnaimgs,'String');
if isempty(protlist)
    msgbox('No images selected for pattern to be added.','Fail to add pattern','error');
    return
end
if ~isempty(dnalist) && length(protlist)~=length(dnalist)
    msgbox('Protein images and DNA images mismatch.','Fail to add pattern','error');
    return
end

nameadded = strtrim(get(handles.newPatternName,'String'));
if isempty(nameadded)
    uicontrol(handles.newPatternName);
    msgbox('Pattern name cannot be empty.','Fail to add pattern','error');
    return
end
for i = 1:length(handles.traindata)
    if strcmp(handles.traindata(i).name,nameadded)
        msgbox('Duplicated pattern name.','Fail to add pattern','error');
        return
    end
end
if ~iscell(nameadded)
    nameadded = {nameadded};
end
patternlist = get(handles.patternlist,'String');
patternlist = [patternlist;nameadded];
patternlist = unique(patternlist);
set(handles.patternlist,'Value',length(patternlist));
set(handles.patternlist,'Max',length(patternlist));
set(handles.patternlist,'String',patternlist);
set(handles.trainmodel,'Enable','on');
endex = length(handles.traindata) + 1;
handles.traindata(endex).name = nameadded{1};
handles.traindata(endex).prot = protlist;
handles.traindata(endex).dna = dnalist;
handles.traindata(endex).protch = uint8(get(handles.protCh,'Value'));
handles.traindata(endex).dnach = uint8(get(handles.dnaCh,'Value'));
set(handles.protimgs,'Value',0);
set(handles.protimgs,'Max',0);
set(handles.protimgs,'String',{});
set(handles.dnaimgs,'Value',0);
set(handles.dnaimgs,'Max',0);
set(handles.dnaimgs,'String',{});
set(handles.newPatternName,'String','');
enablePatternListControls(handles, 'off');
set(handles.predefPatternList,'Value',1);
%guidata(hObject, handles);

% --- Executes on selection change in predefPatternList.
function predefPatternList_Callback(hObject, eventdata, handles)
% hObject    handle to predefPatternList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns predefPatternList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from predefPatternList
index = get(hObject,'Value');
if (length(index) && index(end)>0)
    ptrns = get(hObject,'String');
    set(handles.newPatternName,'String',ptrns(index(end)));
    enablePatternListControls(handles, 'on');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function predefPatternList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to predefPatternList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function newPatternName_Callback(hObject, eventdata, handles)
% hObject    handle to newPatternName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of newPatternName as text
%        str2double(get(hObject,'String')) returns contents of newPatternName as a double

% --- Executes during object creation, after setting all properties.
status = 'off';
if(~isempty((strtrim(get(hObject,'String')))))
    status = 'on'; 
end
enablePatternListControls(handles, status);
guidata(hObject, handles);

function newPatternName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to newPatternName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in patternlist.
function patternlist_Callback(hObject, eventdata, handles)
% hObject    handle to patternlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns patternlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from patternlist


% --- Executes during object creation, after setting all properties.
function patternlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to patternlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function modelname_Callback(hObject, eventdata, handles)
% hObject    handle to modelname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of modelname as text
%        str2double(get(hObject,'String')) returns contents of modelname as a double

% --- Executes during object creation, after setting all properties.
function modelname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modelname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in delpattern.
function delpattern_Callback(hObject, eventdata, handles)
% hObject    handle to delpattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
patternlist = get(handles.patternlist,'String');
delidx = get(handles.patternlist,'Value');

if delidx
    if(length(patternlist)==length(delidx))
        set(handles.trainmodel,'Enable','off');
    end
    patterndel = patternlist{delidx};
    patternlist(delidx) = [];
    set(handles.patternlist,'Value',length(patternlist));
    set(handles.patternlist,'Max',length(patternlist));
    set(handles.patternlist,'String',patternlist);
    for i = 1:length(handles.traindata)
        if strcmp(handles.traindata(i).name,patterndel)
            handles.traindata(i) = [];
        end
    end
end
guidata(hObject, handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over protimgs.
function protimgs_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to protimgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'String');
imgurl = contents{get(hObject,'Value')};
img = imread(imgurl);
figure, imshow(img,[]);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over dnaimgs.
function dnaimgs_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to dnaimgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'String');
imgurl = contents{get(hObject,'Value')};
img = imread(imgurl);
figure, imshow(img,[]);


% --- Executes on button press in clearall.
function clearall_Callback(hObject, eventdata, handles)
% hObject    handle to clearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.protimgs,'Value',0);
set(handles.protimgs,'Max',0);
set(handles.protimgs,'String','');
set(handles.dnaimgs,'Value',0);
set(handles.dnaimgs,'Max',0);
set(handles.dnaimgs,'String','');
set(handles.patternlist,'Value',0);
set(handles.patternlist,'Max',0);
set(handles.patternlist,'String','');
set(handles.trainmodel,'Enable','off');
set(handles.newPatternName,'String','');
set(handles.predefPatternList,'Value',1);
set(handles.modelname,'String','');
set(handles.maxNumObjTypes,'String','');
set(handles.outlierDetThreshold,'String','');
handles.traindata = [];
guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over newPatternName.
function newPatternName_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to newPatternName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in trainmodel.
function trainmodel_Callback(hObject, eventdata, handles)
% hObject    handle to trainmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
modelname = get(handles.modelname,'String');
if (strcmp(strtrim(modelname),''))
    modelname = datestr(now);
end

maxObjects = uint8(str2double(get(handles.maxNumObjTypes,'String')));
if (maxObjects < 2)
    uicontrol(handles.maxNumObjTypes);
    msgbox('Maximum number of object types must be an integer >= 2','Invalid Input','error');
    
    return
end

set(hObject,'Enable','off');

[modelFile, retrain] = gnf_trainmixmodel(handles.traindata,modelname,handles.defaultdir,maxObjects);

if(retrain)
    set(hObject,'Enable','on');
    msgbox('Please increase the maximum number of object types','Low Maximum Number of Object Types','error');
    return
else
    handles = loadModels({modelFile},handles);
    set(handles.patternlist,'Value',0);
    set(handles.patternlist,'Max',0);
    set(handles.patternlist,'String','');
    if(get(handles.showBarGraphs,'Value'))
        
        curFileName = char(modelFile);
        try
            objlearned = load(curFileName);
            gnf_plotobjfractions(objlearned.model);
        catch exception
            %disp();
        end
    end
end
guidata(hObject, handles);


% --- Executes on button press in unmixing.
function unmixing_Callback(hObject, eventdata, handles)
% hObject    handle to unmixing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
protlist = get(handles.protimgs,'String');
dnalist = get(handles.dnaimgs,'String');
if isempty(protlist)
    msgbox('No images selected for pattern to be added.','Data failure','error');
    return
end
if ~isempty(dnalist) && length(protlist)~=length(dnalist)
    msgbox('Protein images and DNA images mismatch.','Data failure','error');
    return
end

objDetThreshold = str2double(get(handles.outlierDetThreshold,'String'));
if isnan(objDetThreshold)
    uicontrol(handles.objDetThreshold);
    msgbox('Threshold for object detection must be a number','Invalid Input','error');
    return
end

mUnmixRMSerr = str2double(get(handles.maxUnmixRMSerror,'String'));
if isnan(mUnmixRMSerr)
    uicontrol(handles.maxUnmixRMSerror);
    msgbox('Maximum unmixing RMS error must be a number','Invalid Input','error');
    return
end

modellist = get(handles.modellist,'String');
index = get(handles.modellist,'Value');
modelname = modellist{index};
for i = 1:length(handles.models)
    if strcmp(modelname,handles.models(i).name)
        index = i;
    end
end

type = get(handles.unmixType,'Value');
if(type == length(get(handles.unmixType,'String')))
    type = 1:length(get(handles.unmixType,'String'))-1;
end
resdata = [];
handles.mixturedatafiles = [];
protch = uint8(str2double(get(handles.protCh,'String')));
dnach = uint8(str2double(get(handles.dnaCh,'String')));
if(get(handles.individualUnmix,'Value'))
    for ind = 1:length(protlist)
        data.prot = {protlist{ind}};
        data.protch = protch;
        data.dna = [];
        datafiles = protlist{ind};
        if(~isempty(dnalist))
            data.dna = {dnalist{ind}};
            data.dnach = dnach;
            datafiles = strcat(protlist{ind},',',dnalist{ind});
        end
        [alpha,concentration] = gnf_unmixing(handles.models(index),data,objDetThreshold,mUnmixRMSerr);
        
        for rowIdx = 1:length(type)
            resdata = [resdata; alpha{rowIdx}];
            handles.mixturedatafiles = [handles.mixturedatafiles; datafiles];
        end
    end
    handles.resultRowNames = repmat(get(handles.result,'RowName'),length(protlist),1);
else
    handles.testdata.prot = protlist;
    handles.testdata.dna = dnalist;
    handles.testdata.protch = protch;
    handles.testdata.dnach = dnach;
    
    [alpha,concentration] = gnf_unmixing(handles.models(index),handles.testdata,objDetThreshold,mUnmixRMSerr);
    for rowIdx = 1:length(type)
        resdata = [resdata; alpha{rowIdx}];
    end
    datafiles = protlist;
    if(~isempty(dnalist))
        datafiles = strcat(protlist,',',dnalist);
    end
    handles.mixturedatafiles = char(datafiles);
    handles.resultRowNames = get(handles.result,'RowName');
end
handles.resdata = resdata;
set(handles.result,'RowName',handles.resultRowNames);
set(handles.result,'Data',resdata);
set(handles.protimgs,'Value',0);
set(handles.protimgs,'Max',0);
set(handles.protimgs,'String',{});
set(handles.dnaimgs,'Value',0);
set(handles.dnaimgs,'Max',0);
set(handles.dnaimgs,'String',{});
set(handles.save,'Enable','on');
guidata(hObject, handles);


% --- Executes on button press in clearresult.
function clearresult_Callback(hObject, eventdata, handles)
% hObject    handle to clearresult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.result,'ColumnName',{'Pattern 1';'Pattern 2';'...'});
set(handles.result,'Data',['','';'','';'','']);
set(handles.save,'Enable','off');
guidata(hObject, handles);


% --- Executes on selection change in modellist.
function modellist_Callback(hObject, eventdata, handles)
% hObject    handle to modellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns modellist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modellist
setTableColumns(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function modellist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in delmodel.
function delmodel_Callback(hObject, eventdata, handles)
% hObject    handle to delmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = deleteModels(get(handles.modellist,'Value'),handles,1);

guidata(hObject, handles);

function handles = deleteModels(delidx,handles,perm)
modellist = get(handles.modellist,'String');
if delidx
    modeldel = modellist{delidx};
    modeldelidx = [];
    for i = 1:length(handles.models)
        if find(strcmp(handles.models(i).name,{modellist{delidx}}))
            modeldelidx(end+1) = i;
            filename = char(handles.modelFileList{i});
            if (perm && exist(filename,'file'))
                try
                    delete(filename);
                catch exception
                    disp('Error: Could not delete the model file:');
                    disp(filename);
                end
            end
        end
    end
    modellist(delidx) = [];
    set(handles.modellist,'Value',length(modellist));
    set(handles.modellist,'Max',length(modellist));
    set(handles.modellist,'String',modellist);
    handles.models(modeldelidx) = [];
    handles.modelFileList{modeldelidx} = [];
end

% --- Executes on button press in omeroAddDna.
function omeroAddDna_Callback(hObject, eventdata, handles)
% hObject    handle to omeroAddDna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (isInvalidCh(handles.dnaCh))
    return
end

filelist = get(handles.dnaimgs,'String');
[filelistadd,fpath] = omeroImageSelectionDialog('dna',uint8(get(handles.dnaCh,'Value')));
if fpath~=0
    handles.lastdir = fpath;
end
filelist = [filelist;filelistadd];
filelist = unique(filelist);
set(handles.dnaimgs,'Value',length(filelist));
set(handles.dnaimgs,'Max',length(filelist));
set(handles.dnaimgs,'String',filelist);
guidata(hObject, handles);





% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over omeroAddProt.
function omeroAddProt_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to omeroAddProt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function maxNumObjTypes_Callback(hObject, eventdata, handles)
% hObject    handle to maxNumObjTypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxNumObjTypes as text
%        str2double(get(hObject,'String')) returns contents of maxNumObjTypes as a double


% --- Executes during object creation, after setting all properties.
function maxNumObjTypes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxNumObjTypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outlierDetThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to outlierDetThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outlierDetThreshold as text
%        str2double(get(hObject,'String')) returns contents of outlierDetThreshold as a double


% --- Executes during object creation, after setting all properties.
function outlierDetThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outlierDetThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in deleteAll.
function deleteAll_Callback(hObject, eventdata, handles)
% hObject    handle to deleteAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = deleteModels(1:length(handles.models),handles,1);

guidata(hObject, handles);




function maxUnmixRMSerror_Callback(hObject, eventdata, handles)
% hObject    handle to maxUnmixRMSerror (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxUnmixRMSerror as text
%        str2double(get(hObject,'String')) returns contents of maxUnmixRMSerror as a double


% --- Executes during object creation, after setting all properties.
function maxUnmixRMSerror_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxUnmixRMSerror (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addToList.
function addToList_Callback(hObject, eventdata, handles)
% hObject    handle to addToList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = addPatternToList(handles);
guidata(hObject, handles);

function handles = addPatternToList(handles)
nameadded = strtrim(get(handles.newPatternName,'String'));
if isempty(nameadded)
    uicontrol(handles.newPatternName);
    msgbox('Pattern name cannot be empty.','Fail to add pattern','error');
    return
end
patternlist = get(handles.predefPatternList,'String');
patternlist{end+1} = char(nameadded);
%patternlist = [patternlist;nameadded];
patternlist = unique(patternlist);
set(handles.predefPatternList,'Value',length(patternlist));
set(handles.predefPatternList,'Max',length(patternlist));
set(handles.predefPatternList,'String',patternlist);
set(handles.saveList,'Enable','on');

% --- Executes on button press in addToBoth.
function addToBoth_Callback(hObject, eventdata, handles)
% hObject    handle to addToBoth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = addPatternToList(handles);
handles = addPatternToModel(hObject,handles);
guidata(hObject, handles);


% --- Executes on button press in deleteFromList.
function deleteFromList_Callback(hObject, eventdata, handles)
% hObject    handle to deleteFromList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
patternlist = get(handles.predefPatternList,'String');
delidx = get(handles.predefPatternList,'Value');

if delidx
    patternlist(delidx) = [];
    set(handles.predefPatternList,'Value',length(patternlist));
    set(handles.predefPatternList,'Max',length(patternlist));
    set(handles.predefPatternList,'String',patternlist);
end
set(handles.saveList,'Enable','on');
guidata(hObject, handles);


% --- Executes on button press in saveList.
function saveList_Callback(hObject, eventdata, handles)
% hObject    handle to saveList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
preDefPatternList = get(handles.predefPatternList,'String');
save('ptrn','preDefPatternList');
set(handles.saveList,'Enable','off');
guidata(hObject, handles);



% --- Executes on key press with focus on newPatternName and none of its controls.
function newPatternName_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to newPatternName (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
enablePatternListControls(handles, 'on');
guidata(hObject, handles);

function setTableColumns(handles)
index = get(handles.modellist,'Value');
if(index > 0)
    set(handles.result,'ColumnName',[handles.models(index).patternlist; 'Unknown']);
    set(handles.result,'Data',[]);
end



function enablePatternListControls(handles, status)
    set(handles.addToList,'Enable',status);
    set(handles.addpattern,'Enable',status);
    set(handles.addToBoth,'Enable',status);


% --- Executes on button press in about.
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
About('version',handles.PUversion);


% --- Executes on button press in setDir.
function setDir_Callback(hObject, eventdata, handles)
% hObject    handle to setDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dirname = uigetdir(handles.defaultdir,'Default directory for loading and storing trained models');

if (ischar(dirname))
    handles = deleteModels(1:length(handles.models),handles,0);
    handles.defaultdir = dirname;
    storedModels = dir(strcat(handles.defaultdir,filesep,'*.mat'));
    modelFiles = [];
    for modelIdx = 1:length(storedModels)
        modelFiles{modelIdx} = strcat(handles.defaultdir,filesep,storedModels(modelIdx).name);
    end
    handles = loadModels(modelFiles,handles);
end


guidata(hObject, handles);




% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile({'*.*'},'Save results as');
if (ischar(filename))
    filename = regexprep(filename,'\.csv','','ignorecase');
    filename = regexprep(filename,'\.\*','','ignorecase');
    
    [fid,msg] = fopen(strcat(pathname,filesep,filename,'.tsv'),'w');
    colNames = get(handles.result,'ColumnName');
    for colIdx = 1:length(colNames)
        fprintf(fid, '%s \t',colNames{colIdx});
    end
    fprintf(fid, 'Unmix \t Image\n');
    
    [numR,d] = size(handles.resdata);
    for dataIdx = 1:numR
        fprintf(fid, '%s \t %s \t %s \n', num2str(handles.resdata(dataIdx,:))...
            ,handles.resultRowNames{dataIdx},handles.mixturedatafiles(dataIdx,:));
    end
    [numF,d] = size(handles.mixturedatafiles);
    if(numF > numR)
        for dataIdx = 2:numF
            fprintf(fid, '  \t   \t %s \n',handles.mixturedatafiles(dataIdx,:));
        end
    end
    fclose(fid);
    
    csvfile = strcat(pathname,filesep,filename,'.csv');
    csvwrite(csvfile,handles.resdata);
end

% --- Executes on selection change in unmixType.
function unmixType_Callback(hObject, eventdata, handles)
% hObject    handle to unmixType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = setResultTableRowNames(handles);
guidata(hObject, handles);

function handles = setResultTableRowNames(handles)
selection = get(handles.unmixType,'Value');
options = get(handles.unmixType,'String') ;

rows = [];
if(selection == length(options))
    for opInd = 1:length(options)-1
        rows{end+1} = options{opInd};
    end
else
    rows = {options{selection}};
end
set(handles.result,'RowName',rows);
set(handles.result,'Data',[]);
set(handles.save,'Enable','off');


% Hints: contents = get(hObject,'String') returns unmixType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from unmixType


% --- Executes during object creation, after setting all properties.
function unmixType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unmixType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in individualUnmix.
function individualUnmix_Callback(hObject, eventdata, handles)
% hObject    handle to individualUnmix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of individualUnmix
set(handles.result,'Data',[]);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function result_CreateFcn(hObject, eventdata, handles)
% hObject    handle to result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes during object creation, after setting all properties.
function addedPatternList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to addedPatternList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on selection change in dnaCh.
function dnaCh_Callback(hObject, eventdata, handles)
% hObject    handle to dnaCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns dnaCh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dnaCh


% --- Executes during object creation, after setting all properties.
function dnaCh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dnaCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in protCh.
function protCh_Callback(hObject, eventdata, handles)
% hObject    handle to protCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns protCh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protCh


% --- Executes during object creation, after setting all properties.
function protCh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in predefPatternList.
function listbox20_Callback(hObject, eventdata, handles)
% hObject    handle to predefPatternList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns predefPatternList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from predefPatternList


% --- Executes during object creation, after setting all properties.
function listbox20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to predefPatternList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deleteFromList.
function pushbutton93_Callback(hObject, eventdata, handles)
% hObject    handle to deleteFromList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in saveList.
function pushbutton94_Callback(hObject, eventdata, handles)
% hObject    handle to saveList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to newPatternName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of newPatternName as text
%        str2double(get(hObject,'String')) returns contents of newPatternName as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to newPatternName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addpattern.
function pushbutton90_Callback(hObject, eventdata, handles)
% hObject    handle to addpattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in addToList.
function pushbutton91_Callback(hObject, eventdata, handles)
% hObject    handle to addToList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in addToBoth.
function pushbutton92_Callback(hObject, eventdata, handles)
% hObject    handle to addToBoth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox23.
function listbox23_Callback(hObject, eventdata, handles)
% hObject    handle to listbox23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox23 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox23


% --- Executes during object creation, after setting all properties.
function listbox23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton103.
function pushbutton103_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton103 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton104.
function pushbutton104_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton104 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton105.
function pushbutton105_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton105 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in protimgs.
function listbox24_Callback(hObject, eventdata, handles)
% hObject    handle to protimgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns protimgs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protimgs


% --- Executes during object creation, after setting all properties.
function listbox24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protimgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in omeroAddProt.
function pushbutton106_Callback(hObject, eventdata, handles)
% hObject    handle to omeroAddProt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in add1.
function pushbutton107_Callback(hObject, eventdata, handles)
% hObject    handle to add1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in delete1.
function pushbutton111_Callback(hObject, eventdata, handles)
% hObject    handle to delete1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clear1.
function pushbutton112_Callback(hObject, eventdata, handles)
% hObject    handle to clear1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Help('version',handles.PUversion);


% --- Executes on button press in showBarGraphs.
function showBarGraphs_Callback(hObject, eventdata, handles)
% hObject    handle to showBarGraphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showBarGraphs


