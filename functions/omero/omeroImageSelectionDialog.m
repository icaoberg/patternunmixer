function varargout = omeroImageSelectionDialog(varargin)
% OMEROIMAGESELECTIONDIALOG M-file for omeroImageSelectionDialog.fig
%      OMEROIMAGESELECTIONDIALOG, by itself, creates a new OMEROIMAGESELECTIONDIALOG or raises the existing
%      singleton*.
%
%      H = OMEROIMAGESELECTIONDIALOG returns the handle to a new OMEROIMAGESELECTIONDIALOG or the handle to
%      the existing singleton*.
%
%      OMEROIMAGESELECTIONDIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OMEROIMAGESELECTIONDIALOG.M with the given input arguments.
%
%      OMEROIMAGESELECTIONDIALOG('Property','Value',...) creates a new OMEROIMAGESELECTIONDIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before omeroImageSelectionDialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to omeroImageSelectionDialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help omeroImageSelectionDialog

% Last Modified by GUIDE v2.5 05-Jul-2010 13:54:49

% Author: Abdul-Saboor Sheikh 
% Created: May, 2010
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
                   'gui_OpeningFcn', @omeroImageSelectionDialog_OpeningFcn, ...
                   'gui_OutputFcn',  @omeroImageSelectionDialog_OutputFcn, ...
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

if exist('omero.client','class') == 0
	try	
		loadOmero;
	catch ME
		errordlg('Could not find or load OMERO Matlab library','Error','modal');
		return;
	end
end



% --- Executes just before omeroImageSelectionDialog is made visible.
function omeroImageSelectionDialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to omeroImageSelectionDialog (see VARARGIN)

% Choose default command line output for omeroImageSelectionDialog
handles.output = hObject;

handles.PWG = '';
handles.cachedImages = {};
handles.channel = 'prot';
handles.channelNum = 0;
if(nargin >= 3)
    for index = 1:2:(nargin-3),
        if strcmpi(varargin{index},'dna')
             handles.channel = 'dna';
        elseif isinteger(varargin{index})
            handles.channelNum = uint8(varargin{index}) - 1;
        end
         if nargin-3==index, break, end
         
    end
end


if(exist('omero-conn-conf.mat','file'))
    s = load('omero-conn-conf');
    set(handles.hostname,'String',s.h);
    set(handles.username,'String',s.u);
    set(handles.port,'String',s.p);
end

% Update handles structure
guidata(hObject, handles);

% Make the GUI modal
%set(handles.Load,'WindowStyle','modal')

% UIWAIT makes omeroImageSelectionDialog wait for user response (see UIRESUME)
uiwait(handles.Load);


% --- Outputs from this function are returned to the command line.
function varargout = omeroImageSelectionDialog_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.cachedImages;
varargout{2} = 0;

% The figure can be deleted now
delete(handles.Load);

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Load_CloseRequestFcn(handles.Load, eventdata, handles);


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imgs = get(handles.images, 'String');

if (~isempty(imgs))
    selProj = get(handles.projects, 'Value');
    projNames = get(handles.projects, 'String');
    selDataset = get(handles.datasets, 'Value');
    dsNames = get(handles.datasets, 'String');
    cacheDir = regexprep(strcat(projNames(selProj), filesep, dsNames(selDataset), filesep,handles.channel), '\s', '_');
    if(~exist('@cacheDir','dir'))
        mkdir(char(cacheDir));
    end
    
    %retrieveCh = 1;
    %if((strfind(lower(handles.channel), 'dna')))
    %    retrieveCh = 0;
    %end

    curDatasetsLinks = handles.projectsList.get(get(handles.projects, 'Value')-1).copyDatasetLinks();
    curImagesLinks = curDatasetsLinks.get(get(handles.datasets, 'Value')-1).getChild().linkedImageList;
    handles.cachedImages = {};
    selImagesIdxes = get(handles.images, 'Value');
    [dummy, numIdx] = size(selImagesIdxes);
    for selIdx = 1:numIdx,
        imageName = curImagesLinks.get(selImagesIdxes(selIdx)-1).getName().getValue();
        imageId = curImagesLinks.get(selImagesIdxes(selIdx)-1).getId().getValue();
        pixels = handles.gateway.getPixels(imageId);
        for z=1:pixels.getSizeZ().getValue(),
            for t=1:pixels.getSizeT().getValue(),
                rawPlane = handles.gateway.getPlane(imageId, z-1, handles.channelNum, t-1);
                plane = toMatrix(rawPlane, pixels);
                handles.cachedImages = [handles.cachedImages; strcat(pwd,filesep,cacheDir,filesep,char(imageName))];
                imwrite(plane,char(strcat(cacheDir,filesep,char(imageName))));
            end
        end
    end
    guidata(hObject, handles);
    
    if(get(handles.stayConnected,'Value'))
        omeroKeepAlive(handles.client);
    else
        unloadOmero;
    end
    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.Load);
    %closereq;
end



% --- Executes on button press in connect.
function connect_Callback(hObject, eventdata, handles)
% hObject    handle to connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hName = strtrim(get(handles.hostname,'String'));
p = strtrim(get(handles.port,'String'));
uName = strtrim(get(handles.username,'String'));
pWord = strtrim(handles.PWG);
if (isempty(hName) || isempty(p) +...
	isempty(uName) || isempty(pWord))
	errordlg('All the fields are required for establishing a connection','Error','modal');
	return;
end


[pNum, status] = str2num(p);
if (~status)
	errordlg('Port value has to be numerical','Error','modal');
	return;
end

try
	handles.client = omero.client(hName,pNum );
	session = handles.client.createSession(uName, pWord);
	handles.gateway = session.createGateway();
    
catch exception
	errordlg('Error occured while trying to connect to the specified OMERO server. Please verify your connection information','Error','modal');
	return;
end

% Get the projects 
handles.projectsList = handles.gateway.getProjects([], true);
%handles.cachedImages = {};
projectNames = {};
if (handles.projectsList.size()>0)
    projectNames = getProjectNames(handles.projectsList);
else
    errordlg('No images available. Your user has no projects on the server','Error','modal');
	%return;
end;
set(handles.projects,'String',projectNames); 
guidata(hObject, handles);
loadDatasets(hObject,handles);


function loadDatasets(hObject,handles)
datasetsNames = {};
if(handles.projectsList.size()>0)
    handles.curDatasetsLinks = handles.projectsList.get(get(handles.projects, 'Value')-1).copyDatasetLinks();
    if (handles.curDatasetsLinks.size()> 0)
        datasetsNames = getDatasetsNames(handles.curDatasetsLinks);
    else
        errordlg('No images available. The selected project contains no dataset','Error','modal');
    end;
end;
set(handles.datasets,'String',datasetsNames);
loadImages(hObject,handles);
guidata(hObject, handles);

function loadImages(hObject,handles)
imagesNames = {};
if(handles.curDatasetsLinks.size()>0)
    handles.curImagesLinks = handles.curDatasetsLinks.get(get(handles.datasets, 'Value')-1).getChild().linkedImageList;
    if (handles.curImagesLinks.size()> 0)
        imagesNames = getImagesNames(handles.curImagesLinks);
    else
        errordlg('No images were found in the selected dataset','Error','modal');
    end;
end;
set(handles.images,'String',imagesNames);
guidata(hObject, handles);

% --- Executes on button press in rememberSettings.
function rememberSettings_Callback(hObject, eventdata, handles)
% hObject    handle to rememberSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = get(handles.hostname,'String');
u = get(handles.username,'String');
p = get(handles.port,'String');
save('omero-conn-conf','h','u','p');




function hostname_Callback(hObject, eventdata, handles)
% hObject    handle to hostname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hostname as text
%        str2double(get(hObject,'String')) returns contents of hostname as a double


% --- Executes during object creation, after setting all properties.
function hostname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hostname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function username_Callback(hObject, eventdata, handles)
% hObject    handle to hostname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hostname as text
%        str2double(get(hObject,'String')) returns contents of hostname as a double


% --- Executes during object creation, after setting all properties.
function username_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hostname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function password_Callback(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of password as text
%        str2double(get(hObject,'String')) returns contents of password as a double


% --- Executes during object creation, after setting all properties.
function password_CreateFcn(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function port_Callback(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of port as text
%        str2double(get(hObject,'String')) returns contents of port as a double


% --- Executes during object creation, after setting all properties.
function port_CreateFcn(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stayConnected.
function stayConnected_Callback(hObject, eventdata, handles)
% hObject    handle to stayConnected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stayConnected


% --- Executes on selection change in projects.
function projects_Callback(hObject, eventdata, handles)
% hObject    handle to projects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loadDatasets(hObject,handles);
guidata(hObject, handles);

% Hints: contents = get(hObject,'String') returns projects contents as cell array
%        contents{get(hObject,'Value')} returns selected item from projects


% --- Executes during object creation, after setting all properties.
function projects_CreateFcn(hObject, eventdata, handles)
% hObject    handle to projects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in datasets.
function datasets_Callback(hObject, eventdata, handles)
% hObject    handle to datasets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.projectsList.size()>0)
    handles.curDatasetsLinks = handles.projectsList.get(get(handles.projects, 'Value')-1).copyDatasetLinks();
    loadImages(hObject,handles);
    guidata(hObject, handles);
end;


% Hints: contents = get(hObject,'String') returns datasets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from datasets


% --- Executes during object creation, after setting all properties.
function datasets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datasets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in connect.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in rememberSettings.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to rememberSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to hostname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hostname as text
%        str2double(get(hObject,'String')) returns contents of hostname as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hostname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of username as text
%        str2double(get(hObject,'String')) returns contents of username as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of password as text
%        str2double(get(hObject,'String')) returns contents of password as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of port as text
%        str2double(get(hObject,'String')) returns contents of port as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stayConnected.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to stayConnected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stayConnected


% --- Executes on selection change in images.
function images_Callback(hObject, eventdata, handles)
% hObject    handle to images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns images contents as cell array
%        contents{get(hObject,'Value')} returns selected item from images


% --- Executes during object creation, after setting all properties.
function images_CreateFcn(hObject, eventdata, handles)
% hObject    handle to images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on password and none of its controls.
function password_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% Deals with user input.
CC = get(handles.Load,'currentcharacter'); % The character user entered.
num = int8(CC);

if num == 13  % This is a carriage return.
    return
end

E = get(hObject,'string');  % the string of the edit box.
% Any key handling other than the return key should be handled
% in the following if else block.
if num == 8  % Backspace pressed, update password and screen.
    set(hObject,'string',E(1:end-1));
    handles.PWG = handles.PWG(1:end-1);
elseif num == 127  % The Delete Key: do nothing.
% On some systems this will delete the symbols.  How would you
% prevent this?
elseif ~isempty(num)
    set(hObject,'string',[E,'*'])  ;  % Print out an asterisk in gui.
    handles.PWG = [handles.PWG CC];
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Load_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes when user attempts to close Load.
function Load_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(handles.Load, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.Load);
else
    % The GUI is no longer waiting, just close it
    delete(handles.Load);
end

% --- Executes on selection change in cachedImagesList.
function cachedImagesList_Callback(hObject, eventdata, handles)
% hObject    handle to cachedImagesList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns cachedImagesList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cachedImagesList


% --- Executes during object creation, after setting all properties.
function cachedImagesList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cachedImagesList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on key press with focus on username and none of its controls.
function username_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to username (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
CC = get(handles.Load,'currentcharacter'); % The character user entered.
num = int8(CC);

if num == 13  % This is a carriage return.
    
end


% --- Executes on key press with focus on hostname and none of its controls.
function hostname_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to hostname (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
CC = get(handles.Load,'currentcharacter'); % The character user entered.
num = int8(CC);

if num == 13  % This is a carriage return.
    
end


% --- Executes on key press with focus on port and none of its controls.
function port_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to port (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
CC = get(handles.Load,'currentcharacter'); % The character user entered.
num = int8(CC);

if num == 13  % This is a carriage return.
    
end

% --- Executes on key press with focus on images and none of its controls.
function images_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to images (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
CC = get(handles.Load,'currentcharacter'); % The character user entered.
num = int8(CC);

if num == 13  % This is a carriage return.
    ok_Callback(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function uipanel7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on key press with focus on Load and none of its controls.
function Load_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

CC = get(handles.Load,'currentcharacter'); % The character user entered.
num = int8(CC);

if num == 27  % This is an esc.
     Load_CloseRequestFcn(hObject, eventdata, handles);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over rememberSettings.
function rememberSettings_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to rememberSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on connect and none of its controls.
function connect_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to connect (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


