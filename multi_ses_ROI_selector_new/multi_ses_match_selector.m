function varargout = multi_ses_match_selector(varargin)
% MULTI_SES_MATCH_SELECTOR MATLAB code for multi_ses_match_selector.fig
%      MULTI_SES_MATCH_SELECTOR, by itself, creates a new MULTI_SES_MATCH_SELECTOR or raises the existing
%      singleton*.
%
%      H = MULTI_SES_MATCH_SELECTOR returns the handle to a new MULTI_SES_MATCH_SELECTOR or the handle to
%      the existing singleton*.
%
%      MULTI_SES_MATCH_SELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTI_SES_MATCH_SELECTOR.M with the given input arguments.
%
%      MULTI_SES_MATCH_SELECTOR('Property','Value',...) creates a new MULTI_SES_MATCH_SELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before multi_ses_match_selector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to multi_ses_match_selector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help multi_ses_match_selector

% Last Modified by GUIDE v2.5 10-Aug-2019 10:48:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @multi_ses_match_selector_OpeningFcn, ...
                   'gui_OutputFcn',  @multi_ses_match_selector_OutputFcn, ...
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


% --- Executes just before multi_ses_match_selector is made visible.
function multi_ses_match_selector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to multi_ses_match_selector (see VARARGIN)

% Choose default command line output for multi_ses_match_selector
handles.output = hObject;

%load variablesi into handles
handles.vars = varargin{1};
handles.assign = varargin{2};
handles.input_log = varargin{3};

%Intiial values
handles.initComp = 1;

%generate 2D logical equal to size of assignment matrix
global input_logical
input_logical = handles.input_log;

%initialize slider
handles.compNb = size(handles.input_log,1);

%set min and max values of slider
handles.slider_ROI.Min = 1;
handles.slider_ROI.Max = handles.compNb;
handles.slider_ROI.Value = 1;
handles.slider_ROI.SliderStep = [1/(handles.compNb-1) , 10/(handles.compNb-1)];

%initialize component text box
handles.text2.String = [num2str(1), '/',num2str(handles.compNb)];

% Update handles structure
guidata(hObject, handles);

%Initialize plots
updatePlots(handles.initComp, handles);

% UIWAIT makes ROI_selector_GUI_V1 wait for user response (see UIRESUME)
uiwait(handles.figure1);
%
function updatePlots(initComp,handles)
%clear all axes
cla(handles.axes1)
cla(handles.axes2)
cla(handles.axes3)
cla(handles.axes4)
cla(handles.axes5)
cla(handles.axes6)
cla(handles.axes7)

%plot BW outline of the component
for ii=1:size(handles.assign,2)
    %check if ROI is nan
    if ~isnan(handles.assign(initComp,ii))
        
        switch ii
            case 1
                axes(handles.axes1);
                
            case 2
                axes(handles.axes2);
                %axis OFF 
            case 3
                axes(handles.axes3);
                %axis OFF 
            case 4
                axes(handles.axes4);
                %axis OFF 
            case 5
                axes(handles.axes5);
                %axis OFF 
            case 6
                axes(handles.axes6);
                %axis OFF 
            case 7
                axes(handles.axes7);
                %axis OFF 
        end
        
        %plot the template
        imagesc(handles.vars(ii).template);
        axis off
        
        hold on
        grayMap = brighten(gray,0.6);
        switch ii
            case 1
                colormap(handles.axes1,grayMap)
            case 2
                colormap(handles.axes2,grayMap)
            case 3
                colormap(handles.axes3,grayMap)
            case 4
                colormap(handles.axes4,grayMap)
            case 5
                colormap(handles.axes5,grayMap)
            case 6
                colormap(handles.axes6,grayMap)
            case 7
                colormap(handles.axes7,grayMap)
        end
        
        %plot centers on top from D1
        %scatter(handles.vars(ii).centers_filt(handles.assign(initComp,ii),2),handles.vars(ii).centers_filt(handles.assign(initComp,ii),1),'y*');
        
        %generate zoom widths around ROI
        handles.xWidthSca(1) = handles.vars(ii).centers_filt(handles.assign(initComp,ii),2)-30;
        handles.xWidthSca(2) = handles.vars(ii).centers_filt(handles.assign(initComp,ii),2)+30;
        handles.yWidthSca(1) = handles.vars(ii).centers_filt(handles.assign(initComp,ii),1)-30;
        handles.yWidthSca(2) = handles.vars(ii).centers_filt(handles.assign(initComp,ii),1)+30;
        %
        xlim([handles.xWidthSca]);
        ylim([handles.yWidthSca]);
        
        %plot componenet outline
        plot(handles.vars(ii).Coor_kp_filt{handles.assign(initComp,ii)}(1,:),handles.vars(ii).Coor_kp_filt{handles.assign(initComp,ii)}(2,:),'r', 'LineWidth',1);
        
    end
end
hold off



% --- Outputs from this function are returned to the command line.
function varargout = multi_ses_match_selector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in remove_ROI_1.
function remove_ROI_1_Callback(hObject, eventdata, handles)
% hObject    handle to remove_ROI_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in next_ROI.
function next_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to next_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.initComp = handles.initComp + 1;

%Update componnent string
handles.text2.String = [num2str(handles.initComp), '/',num2str(handles.compNb)];

%Update slider position
handles.slider_ROI.Value = handles.slider_ROI.Value + 1;

% Update handles structure
guidata(hObject, handles);

%Update plots to reflect selected component
updatePlots(handles.initComp,handles);


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
handles.text2.String = [num2str(handles.initComp), '/',num2str(handles.compNb)];

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


% --- Executes on button press in remove_all_matches.
function remove_all_matches_Callback(hObject, eventdata, handles)
% hObject    handle to remove_all_matches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
