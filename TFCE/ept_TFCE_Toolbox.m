function varargout = ept_TFCE_Toolbox(varargin)
% EPT_TFCE_TOOLBOX MATLAB code for ept_TFCE_Toolbox.fig
%      EPT_TFCE_TOOLBOX, by itself, creates a new EPT_TFCE_TOOLBOX or raises the existing
%      singleton*.
%
%      H = EPT_TFCE_TOOLBOX returns the handle to a new EPT_TFCE_TOOLBOX or the handle to
%      the existing singleton*.
%
%      EPT_TFCE_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EPT_TFCE_TOOLBOX.M with the given input arguments.
%
%      EPT_TFCE_TOOLBOX('Property','Value',...) creates a new EPT_TFCE_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ept_TFCE_Toolbox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ept_TFCE_Toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ept_TFCE_Toolbox

% Last Modified by GUIDE v2.5 03-Sep-2013 12:23:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ept_TFCE_Toolbox_OpeningFcn, ...
                   'gui_OutputFcn',  @ept_TFCE_Toolbox_OutputFcn, ...
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


% --- Executes just before ept_TFCE_Toolbox is made visible.
function ept_TFCE_Toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ept_TFCE_Toolbox (see VARARGIN)

% Choose default command line output for ept_TFCE_Toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ept_TFCE_Toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ept_TFCE_Toolbox_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;






% --- Executes on button press in PB_CreateSummary.
function PB_CreateSummary_Callback(hObject, eventdata, handles)
% hObject    handle to PB_CreateSummary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ept_CreateSummary

% Update handles structure
guidata(hObject, handles);





% --- Executes on button press in PB_EditSummary.
function PB_EditSummary_Callback(hObject, eventdata, handles)
% hObject    handle to PB_EditSummary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ept_EditSummary

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PB_AverageSummaries.
function PB_AverageSummaries_Callback(hObject, eventdata, handles)
% hObject    handle to PB_AverageSummaries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ept_AverageSummaries

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in PB_SubtractSummaries.
function PB_SubtractSummaries_Callback(hObject, eventdata, handles)
% hObject    handle to PB_SubtractSummaries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ept_SubtractSummaries

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PB_ChangeBaseline.
function PB_ChangeBaseline_Callback(hObject, eventdata, handles)
% hObject    handle to PB_ChangeBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ept_ChangeBaseline

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in PB_TFCE.
function PB_TFCE_Callback(hObject, eventdata, handles)
% hObject    handle to PB_TFCE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ept_TFCE

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PB_TFCE_ANOVA.
function PB_TFCE_ANOVA_Callback(hObject, eventdata, handles)
% hObject    handle to PB_TFCE_ANOVA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ept_TFCE_ANOVA

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in PB_ResultViewer.
function PB_ResultViewer_Callback(hObject, eventdata, handles)
% hObject    handle to PB_ResultViewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ept_ResultViewer

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in PB_ExtractROI.
function PB_ExtractROI_Callback(hObject, eventdata, handles)
% hObject    handle to PB_ExtractROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ept_ExtractROI

% Update handles structure
guidata(hObject, handles);
