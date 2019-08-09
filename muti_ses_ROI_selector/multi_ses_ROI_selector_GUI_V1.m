function varargout = ROI_selector_GUI_V1(varargin)
% ROI_SELECTOR_GUI_V1 MATLAB code for ROI_selector_GUI_V1.fig
%      ROI_SELECTOR_GUI_V1, by itself, creates a new ROI_SELECTOR_GUI_V1 or raises the existing
%      singleton*.
%
%      H = ROI_SELECTOR_GUI_V1 returns the handle to a new ROI_SELECTOR_GUI_V1 or the handle to
%      the existing singleton*.
%
%      ROI_SELECTOR_GUI_V1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROI_SELECTOR_GUI_V1.M with the given input arguments.
%
%      ROI_SELECTOR_GUI_V1('Property','Value',...) creates a new ROI_SELECTOR_GUI_V1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ROI_selector_GUI_V1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ROI_selector_GUI_V1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROI_selector_GUI_V1

% Last Modified by GUIDE v2.5 03-Apr-2018 12:36:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ROI_selector_GUI_V1_OpeningFcn, ...
                   'gui_OutputFcn',  @ROI_selector_GUI_V1_OutputFcn, ...
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

% --- Executes just before ROI_selector_GUI_V1 is made visible.
function ROI_selector_GUI_V1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ROI_selector_GUI_V1 (see VARARGIN)

% Choose default command line output for ROI_selector_GUI_V1
handles.output = hObject;
handles.vars = varargin{1};

%Intiial values
handles.initComp = 1;

%Create a logical array for keep track of removed ROIs
global removedROI;
global somaROI;
global dendriteROI;

switch nargin
    case  4
        removedROI = logical(zeros(size(handles.vars.A_keep,2),1));
        somaROI = logical(zeros(size(handles.vars.A_keep,2),1));
        dendriteROI = logical(zeros(size(handles.vars.A_keep,2),1));
    case 5
        removedROI = varargin{end};
        somaROI = logical(zeros(size(handles.vars.A_keep,2),1));
        dendriteROI = logical(zeros(size(handles.vars.A_keep,2),1));
    case 6
        removedROI = varargin{end-1};
        somaROI = varargin{end};
        dendriteROI = logical(zeros(size(handles.vars.A_keep,2),1));
    otherwise
        removedROI = varargin{end-2};
        somaROI = varargin{end-1};
        dendriteROI = varargin{end};
end

%initialize slider
handles.compNb = size(handles.vars.A_keep,2);

%set min and max values of slider
handles.slider_ROI.Min = 1;
handles.slider_ROI.Max = handles.compNb;
handles.slider_ROI.Value = 1;
handles.slider_ROI.SliderStep = [1/(handles.compNb-1) , 10/(handles.compNb-1)];

%initialize component text box
handles.compText.String = [num2str(1), '/',num2str(handles.compNb)];

% Update handles structure
guidata(hObject, handles);

%Initialize plots
updatePlots(handles.initComp, handles);

% UIWAIT makes ROI_selector_GUI_V1 wait for user response (see UIRESUME)
 uiwait(handles.figure1);

function updatePlots(initComp,handles)
%plot BW outline of the component
axes(handles.stackAvg);
imagesc(handles.vars.template);
hold on
grayMap = brighten(gray,0.6);
colormap(handles.stackAvg,grayMap)

%plot centers on top from D1
scatter(handles.vars.centers(initComp,2),handles.vars.centers(initComp,1),'y*');
handles.xWidthSca(1) = handles.vars.centers(initComp,2)-60;
handles.xWidthSca(2) = handles.vars.centers(initComp,2)+60;
handles.yWidthSca(1) = handles.vars.centers(initComp,1)-60;
handles.yWidthSca(2) = handles.vars.centers(initComp,1)+60;
% 
xlim([handles.xWidthSca]);
ylim([handles.yWidthSca]);

%plot componenet outline
plot(handles.vars.Coor_kp{initComp}(1,:),handles.vars.Coor_kp{initComp}(2,:),'r', 'LineWidth',1);

hold off

%plot max correlation image
axes(handles.cn_ROI);
imagesc(handles.vars.Cn);
hold on
colormap(handles.cn_ROI, parula);

%plot centers on top from D1
scatter(handles.vars.centers(initComp,2),handles.vars.centers(initComp,1),'m*');
handles.xWidthSca(1) = handles.vars.centers(initComp,2)-60;
handles.xWidthSca(2) = handles.vars.centers(initComp,2)+60;
handles.yWidthSca(1) = handles.vars.centers(initComp,1)-60;
handles.yWidthSca(2) = handles.vars.centers(initComp,1)+60;
% 
xlim([handles.xWidthSca]);
ylim([handles.yWidthSca]);

%plot componenet outline
plot(handles.vars.Coor_kp{initComp}(1,:),handles.vars.Coor_kp{initComp}(2,:),'r', 'LineWidth',1);

hold off 

%plot exponentially smoothed calcium trace of the component
hold on
axes(handles.trace);
trace_length = size(handles.vars.expDffMedZeroed,2);
plot(handles.vars.expDffMedZeroed(initComp,:));
ylim([-0.5 2]);
xlim([0 trace_length]);
hold off


% --- Outputs from this function are returned to the command line.
function varargout = ROI_selector_GUI_V1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
global removedROI;
global somaROI;
global dendriteROI;
varargout{1} = removedROI;
varargout{2} = somaROI;
varargout{3} = dendriteROI;




% --- Executes on button press in next_ROI.
function next_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to next_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.initComp = handles.initComp + 1;

%Update componnent string
handles.compText.String = [num2str(handles.initComp), '/',num2str(handles.compNb)];

%Update slider position
handles.slider_ROI.Value = handles.slider_ROI.Value + 1;

% Update handles structure
guidata(hObject, handles);

%Update plots to reflect selected component
updatePlots(handles.initComp,handles);

% --- Executes on button press in remove_ROI.
function remove_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to remove_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global removedROI
%remove selected component
removedROI(handles.initComp) = 1;
disp(sprintf('Removing ROI: %d',handles.initComp));

%update text box with action
handles.actionText.String = sprintf('Removing ROI: %d',handles.initComp);

%disp(removedROI);
% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function slider_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to slider_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
positionSld = round(hObject.Value);
handles.initComp = single(positionSld);

%Update componnent string
handles.compText.String = [num2str(handles.initComp), '/',num2str(handles.compNb)];

% Update handles structure
guidata(hObject, handles);

%Update plot to reflect position of slider
updatePlots(handles.initComp, handles);



% --- Executes during object creation, after setting all properties.
function slider_ROI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in keep_ROI.
function keep_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to keep_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global removedROI
%remove selected component
removedROI(handles.initComp) = 0;
disp(sprintf('Keeping ROI: %d',handles.initComp));

%update text box with action
handles.actionText.String = sprintf('Keeping ROI: %d',handles.initComp);

%disp(removedROI);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in soma_ROI_button.
function soma_ROI_button_Callback(hObject, eventdata, handles)
% hObject    handle to soma_ROI_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global somaROI
%remove selected component
somaROI(handles.initComp) = 1;
disp(sprintf('Adding soma ROI: %d',handles.initComp));

%update text box with action
handles.actionText.String = sprintf('Adding soma ROI: %d',handles.initComp);

%disp(removedROI);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in dendrite_ROI_button.
function dendrite_ROI_button_Callback(hObject, eventdata, handles)
% hObject    handle to dendrite_ROI_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dendriteROI
%remove selected component
dendriteROI(handles.initComp) = 1;
disp(sprintf('Adding dendrite ROI: %d',handles.initComp));

%update text box with action
handles.actionText.String = sprintf('Adding dendrite ROI: %d',handles.initComp);

%disp(removedROI);
% Update handles structure
guidata(hObject, handles);
