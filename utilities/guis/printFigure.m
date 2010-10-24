function varargout = printFigure(varargin)
% PRINTFIGURE M-file for printFigure.fig
%      PRINTFIGURE, by itself, creates a new PRINTFIGURE or raises the existing
%      singleton*.
%
%      H = PRINTFIGURE returns the handle to a new PRINTFIGURE or the handle to
%      the existing singleton*.
%
%      PRINTFIGURE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRINTFIGURE.M with the given input arguments.
%
%      PRINTFIGURE('Property','Value',...) creates a new PRINTFIGURE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before printFigure_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to printFigure_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help printFigure

% Last Modified by GUIDE v2.5 18-Jul-2008 08:30:44

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @printFigure_OpeningFcn, ...
                       'gui_OutputFcn',  @printFigure_OutputFcn, ...
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
end

% --- Executes just before printFigure is made visible.
function printFigure_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to printFigure (see VARARGIN)


    % Choose default command line output for printFigure
    handles.output = hObject;

    % UIWAIT makes printFigure wait for user response (see UIRESUME)
    % uiwait(handles.printFigure);

    loc_refreshFigureList(handles);
    
    % Update handles structure
    guidata(hObject, handles);
    
end

% --- Outputs from this function are returned to the command line.
function varargout = printFigure_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end



% --- Executes on button press in bSave.
function bSave_Callback(hObject, eventdata, handles)
    % hObject    handle to bSave (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    pdf = get(handles.chPDF,'Value');
    hFig = loc_handleOfSelectedFigure(handles);
    
    fileName = get(handles.txtFileName,'String');

    sizeMode = 'standard';
    if get(handles.rbStandard,'Value'), sizeMode = 'standard'; end
    if get(handles.rbWeb,'Value'), sizeMode = 'web'; end
    if get(handles.rbPrint,'Value'), sizeMode = 'print'; end
    if get(handles.rbHD,'Value'), sizeMode = 'HD'; end

    printResultsFigure(fileName,sizeMode,hFig,pdf);
end

% --- Executes on button press in bCancel.
function bCancel_Callback(hObject, eventdata, handles)
% hObject    handle to bCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    delete(handles.printFigure);
end

% --- Executes on button press in chPDF.
function chPDF_Callback(hObject, eventdata, handles)
% hObject    handle to chPDF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chPDF
end

% --- Executes on selection change in puFigures.
function puFigures_Callback(hObject, eventdata, handles)
% hObject    handle to puFigures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns puFigures contents as cell array
%        contents{get(hObject,'Value')} returns selected item from puFigures
end

% --- Executes during object creation, after setting all properties.
function puFigures_CreateFcn(hObject, eventdata, handles)
% hObject    handle to puFigures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function txtFileName_Callback(hObject, eventdata, handles)
% hObject    handle to txtFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFileName as text
%        str2double(get(hObject,'String')) returns contents of txtFileName as a double
end

% --- Executes during object creation, after setting all properties.
function txtFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


% --- Executes on button press in bView.
function bView_Callback(hObject, eventdata, handles)
% hObject    handle to bView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    hFig = loc_handleOfSelectedFigure(handles);
    if ~isempty(hFig);
        figure(hFig);
    end
end

function hFig = loc_handleOfSelectedFigure(handles)
% This function determines the selected index of the puFigures object and returns the handle of the
% slected value.
    hPU = handles.puFigures;
    hList = get(hPU,'UserData');
    strVal = get(hPU,'Value');
    hFig = hList(strVal);
end

function loc_refreshFigureList(handles)

    set(handles.bView,'Enable','on');
    set(handles.bSave,'Enable','on');

    hList = findall(0,'Type','figure');
    puTags = {};                    puHandles = {};
    hTags = [];                     hHandles = [];
    iTags = 1;                      iHandles = 1;
    
    for i=1:length(hList)
       tempStr = get(hList(i),'Tag');
       thisTest = regexpi(tempStr,'printFigure');
       guiTest = regexpi(tempStr,'GUI');
       if ( ~isempty(tempStr) && isempty(thisTest) && isempty(guiTest) )
           puTags{iTags,1} = tempStr; 
           hTags(iTags) = hList(i);
           iTags = iTags + 1; 
       elseif isempty(tempStr)
          puHandles{iHandles,1} = ['Handled Figure: ' num2Str(hList(i))];
          hHandles(iHandles) = hList(i);
          iHandles = iHandles + 1;
       end
    end

    if ~isempty(puTags),        figList = puTags;    handleList = hTags;        
    elseif ~isempty(puHandles), figList = puHandles; handleList = hHandles;
    else
        figList = 'No valid figures'; handleList = []; 
        set(handles.bSave,'Enable','off');
        set(handles.bView,'Enable','off');
    end

    set(handles.puFigures,'String',figList);
    set(handles.puFigures,'UserData',handleList);
end


% --- Executes on button press in bRefresh.
function bRefresh_Callback(hObject, eventdata, handles)
% hObject    handle to bRefresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    loc_refreshFigureList(handles)
end
