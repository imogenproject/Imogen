function varargout = Slice3D_GUI(varargin)
% SLICE3D_GUI M-file for Slice3D_GUI.fig
%      SLICE3D_GUI, by itself, creates a new SLICE3D_GUI or raises the existing
%      singleton*.
%
%      H = SLICE3D_GUI returns the handle to a new SLICE3D_GUI or the handle to
%      the existing singleton*.
%
%      SLICE3D_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SLICE3D_GUI.M with the given input arguments.
%
%      SLICE3D_GUI('Property','Value',...) creates a new SLICE3D_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Slice3D_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Slice3D_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Slice3D_GUI

% Last Modified by GUIDE v2.5 12-Jul-2008 14:41:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Slice3D_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Slice3D_GUI_OutputFcn, ...
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


% --- Executes just before Slice3D_GUI is made visible.
function Slice3D_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Slice3D_GUI (see VARARGIN)

% Choose default command line output for Slice3D_GUI
handles.output = hObject;
handles.sourceFigure = varargin{1};
handles.plotMode = varargin{2};

set(hObject,'CloseRequestFcn',{@partneredCloseCB,handles.sourceFigure});

pName = get(handles.sourceFigure,'Tag');
data.pTag = pName;
set(hObject,'Tag',[pName '_GUI'],'Name',[pName '_GUI']);

switch (handles.plotMode)
    case 'slicePlotter3x3'
        maxVals = varargin{3};
        curVals = varargin{4};
        minVal = 0;
        data.array = varargin{5};
    case 'plot3DArray'
        maxVals = size(varargin{4});
        curVals = varargin{3};
        minVal = 1;
        data.array = varargin{4};
        data.contrast = varargin{5};
        data.scaleMax = varargin{6};
end

set(handles.Slice3DFig, 'UserData', data);
set(handles.slX,'Min',minVal,'Max',maxVals(1),'Value',curVals(1),'SliderStep',[1/(maxVals(1) - minVal) 0.1]);
set(handles.slY,'Min',minVal,'Max',maxVals(2),'Value',curVals(2),'SliderStep',[1/(maxVals(2) - minVal) 0.1]);
set(handles.slZ,'Min',minVal,'Max',maxVals(3),'Value',curVals(3),'SliderStep',[1/(maxVals(3) - minVal) 0.1]);

set(handles.txtX,'String',num2str(curVals(1)));
set(handles.txtY,'String',num2str(curVals(2)));
set(handles.txtZ,'String',num2str(curVals(3)));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Slice3D_GUI wait for user response (see UIRESUME)
% uiwait(handles.Slice3DFig);


% --- Outputs from this function are returned to the command line.
function varargout = Slice3D_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slX_Callback(hObject, eventdata, handles)
% hObject    handle to slX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = round(get(hObject,'Value'));
set(hObject,'Value', val);
set(handles.txtX,'String',num2str(val));

% --- Executes during object creation, after setting all properties.
function slX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slY_Callback(hObject, eventdata, handles)
% hObject    handle to slY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = round(get(hObject,'Value'));
set(hObject,'Value', val);
set(handles.txtY,'String',num2str(val));

% --- Executes during object creation, after setting all properties.
function slY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slZ_Callback(hObject, eventdata, handles)
% hObject    handle to slZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = round(get(hObject,'Value'));
set(hObject,'Value', val);
set(handles.txtZ,'String',num2str(val));

% --- Executes during object creation, after setting all properties.
function slZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnUpdate.
function btnUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to btnUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
slice = zeros(1,3);
slice(1) = round(get(handles.slX,'Value'));
slice(2) = round(get(handles.slY,'Value'));
slice(3) = round(get(handles.slZ,'Value'));

% Protect from complete zero
if (sum(slice) < 1); slice(1) = 1; set(handles.slX,'Value',1); set(handles.txtX,'String',1); end

data = get(handles.Slice3DFig,'UserData');
set(0,'CurrentFigure',handles.sourceFigure);
clf(handles.sourceFigure);

switch (handles.plotMode)
    case 'slicePlotter3x3'; slicePlotter3x3(data.pTag,[],slice,data.array);
    case 'plot3DArray'; plot3DArray(data.array, data.contrast, data.scaleMax, data.pTag, slice);
end


