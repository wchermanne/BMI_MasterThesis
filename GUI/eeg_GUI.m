function varargout = eeg_GUI(varargin)
% EEG_GUI MATLAB code for eeg_GUI.fig
%      EEG_GUI, by itself, creates a new EEG_GUI or raises the existing
%      singleton*.
%
%      H = EEG_GUI returns the handle to a new EEG_GUI or the handle to
%      the existing singleton*.
%
%      EEG_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EEG_GUI.M with the given input arguments.
%
%      EEG_GUI('Property','Value',...) creates a new EEG_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before eeg_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eeg_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eeg_GUI

% Last Modified by GUIDE v2.5 10-Feb-2018 14:44:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eeg_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @eeg_GUI_OutputFcn, ...
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


% --- Executes just before eeg_GUI is made visible.
function eeg_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eeg_GUI (see VARARGIN)

filename = 'pic.png';
y = imread(filename);
handles.axes1 = imshow(y);

% Choose default command line output for eeg_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes eeg_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = eeg_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Gen_button.
function Gen_button_Callback(hObject, eventdata, handles)
% hObject    handle to Gen_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fs = str2double(get(handles.fs_edit,'String'));
% if (Fs<= 120)
%     Fs = 120;
% end
% if (Fs>= 500)
%     Fs = 500;
% end 
% handles.fs_edit = Fs;
duration = str2double(get(handles.time_edit,'String'));
SNR = str2double(get(handles.SNR_edit,'String'));
time_start = str2double(get(handles.move_time_edit,'String'));
eye_blink_check = get(handles.eye_blink_check,'Value');
emg_check = get(handles.emg_effects_check,'Value');
power_line = get(handles.power_check,'Value');
movement = get(handles.popupmenu,'Value');
[C3,C4,Cz] = eeg_generator(Fs,duration,SNR,movement, time_start,eye_blink_check,power_line);
save('channels.mat','C3','C4','Cz')

parameters=[];
parameters.duration=duration;
parameters.fsample=Fs;
parameters.timestart=time_start;
paremeters.movement=movement;
paremeters.eye_blink_check=eye_blink_check;
paremeters.emg_check=emg_check;
paremeters.power_line=power_line;
save('parameters.mat', 'parameters')



% --- Executes on button press in eye_blink_check.
function eye_blink_check_Callback(hObject, eventdata, handles)
% hObject    handle to eye_blink_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eye_blink_check


% --- Executes on button press in emg_effects_check.
function emg_effects_check_Callback(hObject, eventdata, handles)
% hObject    handle to emg_effects_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of emg_effects_check


% --- Executes on button press in power_check.
function power_check_Callback(hObject, eventdata, handles)
% hObject    handle to power_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of power_check



function fs_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fs_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_edit as text
       %str2double(get(hObject,'String'))% returns contents of fs_edit as a double


% --- Executes during object creation, after setting all properties.
function fs_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_edit as text
       %str2double(get(hObject,'String'));% returns contents of time_edit as a double


% --- Executes during object creation, after setting all properties.
function time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SNR_edit_Callback(hObject, eventdata, handles)
% hObject    handle to SNR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SNR_edit as text
%        str2double(get(hObject,'String')) returns contents of SNR_edit as a double


% --- Executes during object creation, after setting all properties.
function SNR_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SNR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function move_time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to move_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of move_time_edit as text
%        str2double(get(hObject,'String')) returns contents of move_time_edit as a double


% --- Executes during object creation, after setting all properties.
function move_time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to move_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu.
function popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu


% --- Executes during object creation, after setting all properties.
function popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EOG_check.
function EOG_check_Callback(hObject, eventdata, handles)
% hObject    handle to EOG_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EOG_check


% --- Executes on button press in C3_check.
function C3_check_Callback(hObject, eventdata, handles)
% hObject    handle to C3_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of C3_check


% --- Executes on button press in Cz_check.
function Cz_check_Callback(hObject, eventdata, handles)
% hObject    handle to Cz_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Cz_check


% --- Executes on button press in C4_check.
function C4_check_Callback(hObject, eventdata, handles)
% hObject    handle to C4_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of C4_check
