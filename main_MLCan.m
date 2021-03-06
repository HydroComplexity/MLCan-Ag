function varargout = main_MLCan(varargin)
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                       MAIN INTERFACE PROGRAM                          %%
%%           Canopy-Root-Soil-Atmosphere Exchange Model                  %%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
% This interface is used to connect and facilitate the MLCan matlab code  %
% written by Darren Drewry. See Drewry et al, 2010.                    	  %
% GUI is developed based on MLCan_v01 code written by Phong Le.           %
% See Le 2012.                                                            %
%                                                                         %
% This interface include the following components:                        %
%   + Model Setup   : Used for setting up model information               %
%   + Option        : Used for choosing sub-models                        %
%   + Forcings / Initial conditions                                       %
%                   : Used for creating/choosing forcings and initital    %
%                     conditions for sub-models                           %
%   + Parameters    : Used for entering parameters values                 %
%   + Results       : Showing plot and results from files                 %
%-------------------------------------------------------------------------%
% MAIN_MLCAN M-file for main_MLCan.fig                                    %
%-------------------------------------------------------------------------%
%   Created by      : Dongkook Woo                                        %
%   Date            : May 20, 2012                                        %
%   Last Modified   : June 23, 2016                                       %

%   Added N-allocation by   : Rohit Nandan
%   Date                    : Oct 02, 2021                                        %
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%% Begin initialization code - DO NOT EDIT                               %%
%-------------------------------------------------------------------------%
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @main_MLCan_OpeningFcn, ...
    'gui_OutputFcn',  @main_MLCan_OutputFcn, ...
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
%-------------------------------------------------------------------------%
%% End initialization code - DO NOT EDIT                                 %%
%-------------------------------------------------------------------------%

% --- Executes just before main_MLCan is made visible.
function main_MLCan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_MLCan (see VARARGIN)

% Choose default command line output for main_MLCan

% Add Users directory to call sub-interface
addpath('./users/');
workpath = pwd;
setappdata(0,'workpath',workpath);

addpath('./users/icons/');


handles.output = hObject;
set(handles.mnu_model,'Enable','off');
set(handles.mnu_result,'Enable','off');
set(handles.mnu_save,'Enable','off');
set(handles.mnu_saveas,'Enable','off');
set(handles.toolbar_save,'Enable','off');

set(handles.pan_start_screen,'Visible','on');
set(handles.button_setup_model,'Visible','off');
set(handles.pan_layer,'Visible','off');
set(handles.txt_model_mlcan,'Visible','on');
set(handles.txt_model_full_name,'Visible','on');
set(handles.Main,'Color',[0.502 0.502 0.502]);

% Update handles structure
guidata(hObject, handles);
%-------------------------------------------------------------------------%
% UIWAIT makes main_MLCan wait for user response (see UIRESUME)           %
% uiwait(handles.Main);                                                   %
%-------------------------------------------------------------------------%


% --- Outputs from this function are returned to the command line.
function varargout = main_MLCan_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function mnu_file_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_model_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_result_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_help_sub_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_help_sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model_help()

% --------------------------------------------------------------------
function mnu_undo_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_redo_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_redo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_copy_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_cut_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_paste_Callback(hObject, eventdata, handles)

%-------------------------------------------------------------------------%
% New in menu editor connected to "New function"                          %
%-------------------------------------------------------------------------%
function mnu_new_Callback(hObject, eventdata, handles)
button_new_main_Callback(hObject, eventdata, handles)


%-------------------------------------------------------------------------%
% Open in menu editor connected to "Open function"                        %
%-------------------------------------------------------------------------%
function mnu_open_Callback(hObject, eventdata, handles)
button_open_main_Callback(hObject, eventdata, handles)


%-------------------------------------------------------------------------%
% "Save function" - Save all variables to the current working files       %
%-------------------------------------------------------------------------%
function mnu_save_Callback(hObject, eventdata, handles)
load './Temps/temporary.mat' 'working_name';
copyfile('./Temps/temp_model.mat', working_name,'f');


%-------------------------------------------------------------------------%
% "Save as function" - Save all variables to the other files              %
%-------------------------------------------------------------------------%
function mnu_saveas_Callback(hObject, eventdata, handles)
[filename,pathname] = uiputfile('*.mat','Save your files');
if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
end

fullpath_filename = [pathname, filename];
workpath = getappdata(0,'workpath');
rem_path = strrep(fullpath_filename,workpath,'');
if strcmp(rem_path,fullpath_filename) == 1
    Msgbox('Do NOT save files outside the working folder, please save the file again','MLCan Error');
    return
else
    working_name = ['.',rem_path];
    hgsave(working_name);
    copyfile('./Temps/temp_model.mat', fullpath_filename,'f');
end


%-------------------------------------------------------------------------%
% "Exit function" - Close the program and exit                            %
%-------------------------------------------------------------------------%
function mnu_exit_Callback(hObject, eventdata, handles)
close


% -------------------------------------------------------------------------
function mnu_print_Callback(hObject, eventdata, handles)


%-------------------------------------------------------------------------%
% "Run function" - Connect and run Darren Drewry model                    %
%-------------------------------------------------------------------------%
function mnu_run_Callback(hObject, eventdata, handles)
load './Temps/temp_model.mat' 'working_forcings' 'root_init' 'root_init_litter';

try
    load(working_forcings)
    
catch exception
    msgbox('Please load forcings - go back to FORCINGS & INITIAL CONDITIONS','MLCan Error','Error');
    return
end

load('.\Temps\temporary.mat', 'config_path');

ResPath=fullfile(config_path, 'Results');

if ~isdir(ResPath)
    mkdir (ResPath)
end


DRIVER_CROP_3_0;


% -------------------------------------------------------------------------
function mnu_forcings_condition_Callback(hObject, eventdata, handles)
button_group_condition_Callback(hObject, eventdata, handles)

%-------------------------------------------------------------------------%
% Call parameter interface and connect the model                          %
%-------------------------------------------------------------------------%
function mnu_parameters_Callback(hObject, eventdata, handles)
model_parameters();

% -------------------------------------------------------------------------
function Context_test_Callback(hObject, eventdata, handles)



% -------------------------------------------------------------------------
function toolbar_print_ClickedCallback(hObject, eventdata, handles)



% -------------------------------------------------------------------------
function toolbar_open_ClickedCallback(hObject, eventdata, handles)
button_open_main_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function axe_start_screen_CreateFcn(hObject, eventdata, handles)
background = imread('./users/background.jpg');
image(background);
axis off

%imshow('./users/background.jpg');                                                   % Load image background in axes


%-------------------------------------------------------------------------%
% NEW BUTTON ON START SCREEN TO CREATE NEW FILES...                       %
%-------------------------------------------------------------------------%
function button_new_main_Callback(hObject, eventdata, handles)
[filename, pathname] = uiputfile('*.mat', 'Create new files');
if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
else
    fullpath_filename = [pathname, filename];
    workpath = getappdata(0,'workpath');
    rem_path = strrep(fullpath_filename,workpath,'');
    if strcmp(rem_path,fullpath_filename) == 1
        Msgbox('Do NOT creat project files outside the working folder','MLCan Error');
        return
    else
        working_name = ['.',rem_path];
        config_path = pathname;
        
        set(handles.pan_start_screen,'Visible','off');
        set(handles.Main,'Color',[0.941 0.941 0.941]);
        set(handles.pan_layer,'Visible','on');
        main_menu_enable(hObject, eventdata, handles);
        save './Temps/temporary.mat' 'working_name' 'config_path';
        
        % Empty information used for model setup form.
        
        LAImin_face = '';
        
        lat_face    = 0;
        long_face   = 0;
        elev_face   = 0;
        
        DOY_start   = 1;
        DOY_end   = 365;
        
        % Sub-Model
        CGM=0;
        Turbulence=0;
        HR=0;
        Soil_heat=0;
        Soil_nutrient=0;
        
        N_denit=0;
        N_Adepo=0;
        N_Adepo_amm=0.31;
        N_Adepo_nit=0.25;
        NCanAlMod=0;
        ParOn='Off';
        
        % Canopy Process
        Sim_species     = 1;
        Sim_species_con = 1;
        num_species=Sim_species_con;
        vanGen=1;
        RHC=1;
        CO2_Ambient=400;
        CO2_Elev=0;
        CO2_Elev_con=550;
        Temp_Elev=0;
        Temp_Elev_con=1;
        
        % Yearwise Data
        Start_Y=2008;
        End_Y=2008;
        
        YrFilesList=zeros(-1);
        
        % Soil CN Process
        Soil_C_pool1=1;
        N_Uptake_RB     = 1;
        opt_root_para = {...
            ' diam'     , 0.25,  0.3, 0.35, 0.45,   ; ...
            ' dist'     ,  0.5, 0.24, 0.13, 0.13,   ; ...
            ' SRL'      ,   22,   28,   38,   20,     ...
            };
        N_Fix           = 0;
        N_Remo          = 0;
        
        
        %% Parameters
        root_init           = zeros(3,3);
        root_name           = {'Depth', 'Moisture', 'Temperature'};
        nutrient_int        = zeros(3,3);
        nutrient_name1      = {'Depth', 'Organic Carbon',                'Microbial Carbon', 'Ammonium-N'...
            'Nitrate-N', 'Organic C:N' };
        nutrient_name2      = {'Depth', 'Litter Carbon', 'Humus Carbon', 'Microbial Carbon', 'Ammonium-N'...
            'Nitrate-N', 'Litter C:N'  };
        root_init_litter    = zeros(3,3);
        nutrient_int_litter = zeros(3,3);
        
        dat_decom        = zeros(3,4);
        dat_decom_litter = zeros(1,4);
        
        %Empty information used for parameters.
        fullpath_forcings   = '';
        working_forcings    = '';
        
        %%
        % save empty information for new project-------------------------------
        save (fullpath_filename, ...
            ...%LocationInfo
            'lat_face', 'long_face', 'elev_face', 'DOY_start', 'DOY_end',...
            ...%SubModel
            'CGM', 'Turbulence', 'HR', 'Soil_heat', 'Soil_nutrient', ...
            'N_denit', 'N_Adepo', 'N_Adepo_amm', 'N_Adepo_nit', 'NCanAlMod','ParOn',...
            ...%CanopyProcess
            'Sim_species', 'Sim_species_con', 'num_species', 'vanGen', 'RHC',...
            'CO2_Ambient', 'CO2_Elev', 'CO2_Elev_con', 'Temp_Elev', ...
            'Temp_Elev_con', ...
            ...%YearWise
            'Start_Y', 'End_Y', 'YrFilesList', ...
            ...%SoilCNModel
            'Soil_C_pool1', 'N_Uptake_RB', 'opt_root_para', 'N_Fix',...
            'N_Remo',...
            ...%Others
            'root_init', 'root_name', 'nutrient_int',...
            'nutrient_name1', 'nutrient_name2', 'root_init_litter', ...
            'nutrient_int_litter', 'dat_decom', 'dat_decom_litter',...
            'fullpath_forcings', 'working_forcings')
        
        
        save ('./Temps/temp_model.mat',...
            ...%LocationInfo
            'lat_face', 'long_face', 'elev_face', 'DOY_start', 'DOY_end',...
            ...%SubModel
            'CGM', 'Turbulence', 'HR', 'Soil_heat', 'Soil_nutrient', ...
            'N_denit', 'N_Adepo', 'N_Adepo_amm', 'N_Adepo_nit', 'NCanAlMod', 'ParOn',...
            ...%CanopyProcess
            'Sim_species', 'Sim_species_con', 'num_species', 'vanGen', 'RHC',...
            'CO2_Ambient', 'CO2_Elev', 'CO2_Elev_con', 'Temp_Elev', ...
            'Temp_Elev_con', ...
            ...%YearWise
            'Start_Y', 'End_Y', 'YrFilesList', ...
            ...%SoilCNModel
            'Soil_C_pool1', 'N_Uptake_RB', 'opt_root_para', 'N_Fix',...
            'N_Remo',...
            ...%Others
            'root_init', 'root_name', 'nutrient_int',...
            'nutrient_name1', 'nutrient_name2', 'root_init_litter', ...
            'nutrient_int_litter', 'dat_decom', 'dat_decom_litter',...
            'fullpath_forcings', 'working_forcings')
        
        
    end
end

%-------------------------------------------------------------------------%
% OPEN BUTTON ON START SCREEN TO OPEN FILES...                            %
%-------------------------------------------------------------------------%
function button_open_main_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.mat', 'Open MLCan files');
if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
else
    fullpath_filename = [pathname, filename];
    workpath = getappdata(0,'workpath');
    rem_path = strrep(fullpath_filename,workpath,'');
    if strcmp(rem_path,fullpath_filename) == 1
        Msgbox('Please copy project files inside the working folder before loading','MLCan Error');
        return
    else
        working_name = ['.',rem_path];
        config_path = pathname;
        
        try
            %load(working_name, 'crop_name', 'lat_face', 'long_face', 'LAImin_face', 'num_can', 'num_LAD', 'num_root');
            save './Temps/temporary.mat' 'working_name' 'config_path';
            
            copyfile(fullpath_filename,'./Temps/temp_model.mat','f');
            set(handles.pan_start_screen,'Visible','off');
            set(handles.Main,'Color',[0.941 0.941 0.941]);
            set(handles.pan_layer,'Visible','on');
            main_menu_enable(hObject, eventdata, handles);
        catch exception
            msgbox('Incorrect input file','MLCan: Error');
            return
        end
    end
end


%-------------------------------------------------------------------------%
% EXIT BUTTON...                                                          %
%-------------------------------------------------------------------------%
function button_exit_main_Callback(hObject, eventdata, handles)
close


% --------------------------------------------------------------------
function mnu_about_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model_about()


%-------------------------------------------------------------------------%
% Enable some menu in the main screen                                     %
%-------------------------------------------------------------------------%
function main_menu_enable(hObject, eventdata, handles)
set(handles.mnu_model,'Enable','on');
set(handles.mnu_result,'Enable','on');
set(handles.mnu_save,'Enable','on');
set(handles.mnu_saveas,'Enable','on');
set(handles.toolbar_save,'Enable','on');
set(handles.button_setup_model,'Visible','on');
set(handles.txt_model_mlcan,'Visible','on');
set(handles.txt_model_full_name,'Visible','on');


%-------------------------------------------------------------------------%
% Call MODEL SETUP interface from button in main screen                   %
%-------------------------------------------------------------------------%
function button_group_setup_Callback(hObject, eventdata, handles)
% model_setup;
model_setup_App1

%-------------------------------------------------------------------------%
% Call OPTIONS interface from button in main screen                       %
%-------------------------------------------------------------------------%
function button_group_option_Callback(hObject, eventdata, handles)

load './Temps/temp_model.mat' 'DOY_start' 'DOY_end'
if isstr(DOY_start)
    DOY_start=str2num(DOY_start);
end
if isstr(DOY_end)
    DOY_end=str2num(DOY_end);
end

if isempty(DOY_start) || isempty(DOY_end)
    msgbox( 'Please go back to MODEL SETUP - Set Specified DOY',   ...
        'MLCan Error','Error');
    return
else
    %     model_option;
    model_option_App1
    
    load('.\Temps\temporary.mat', 'config_path');
    
    YrPath=[config_path, 'YearFr'];
    
    if ~isdir(YrPath)
        mkdir (YrPath)
    end
    
end

%-------------------------------------------------------------------------%
% Call FORCINGS/INITIAL CONDITONS interface from button in main screen    %
%-------------------------------------------------------------------------%
function button_group_condition_Callback(hObject, eventdata, handles)

load './Temps/temp_model.mat' 'root_init' 'Start_Y'

file_name   = sprintf('YrwdSpec_%sY%s.mat', num2str(1), num2str(Start_Y));

load('.\Temps\temporary.mat', 'config_path');

YrPath=fullfile(config_path, 'YearFr');

if ~isdir(YrPath)
    mkdir (YrPath)
end

fullpath    =...
    fullfile(YrPath, file_name);

load (fullpath, 'num_root', 'dat_root');

% try
dat_length = length(root_init(:,1));
if isstr(num_root) == 1
    root_num = str2num(num_root);
else
    root_num = num_root;
end
if dat_length > root_num
    root_init = root_init(1:root_num,:);
elseif dat_length < root_num
    root_init_add = zeros(root_num-dat_length,3);
    root_init = [root_init; root_init_add];
end
root_init = [dat_root(:,1),root_init(:,2:3)];
model_forcings;
% catch exception
%     msgbox( 'Please go back to MODEL SETUP - Click on Add/Edit Input root profile and/or Leaf area density',   ...
%         'MLCan Error','Error');
%     return
% end

%-------------------------------------------------------------------------%
% Call PARAMETERS interface from button in main screen                    %
%-------------------------------------------------------------------------%
function button_group_parameters_Callback(hObject, eventdata, handles)
% model_parameters;
model_parameters_App1;

load('.\Temps\temporary.mat', 'config_path');

PrPath=fullfile(config_path, 'Par');

if ~isdir(PrPath)
    mkdir (PrPath)
end


%-------------------------------------------------------------------------%
% Call RESULTS interface from button in main screen                       %
%-------------------------------------------------------------------------%
function button_group_result_Callback(hObject, eventdata, handles)
model_results;
%-------------------------------------------------------------------------%
% Loading icon for buttons in main interface                              %

function axes3_CreateFcn(hObject, eventdata, handles)
modelicon = imread('./users/icons/wheat.png');
image(modelicon)
axis off
%imshow('./users/icons/wheat.png');

function axes4_CreateFcn(hObject, eventdata, handles)
optionicon = imread('./users/icons/option.png');
image(optionicon)
axis off
%imshow('./users/icons/option.png');

function axes5_CreateFcn(hObject, eventdata, handles)
forceicon = imread('./users/icons/forcings.png');
image(forceicon)
axis off
%imshow('./users/icons/forcings.png');

function axes6_CreateFcn(hObject, eventdata, handles)
paraicon = imread('./users/icons/parameter.png');
image(paraicon)
axis off
%imshow('./users/icons/parameter.png');

function axes7_CreateFcn(hObject, eventdata, handles)
resulticon = imread('./users/icons/result.png');
image(resulticon)
axis off
%imshow('./users/icons/result.png');

% End of loading icon                                                     %
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Toolbar functions                                                       %
%-------------------------------------------------------------------------%
function toolbar_new_ClickedCallback(hObject, eventdata, handles)
button_new_main_Callback(hObject, eventdata, handles)

function toolbar_save_ClickedCallback(hObject, eventdata, handles)
mnu_save_Callback(hObject, eventdata, handles)

function mnu_setup_model_Callback(hObject, eventdata, handles)
button_group_setup_Callback(hObject, eventdata, handles)

function mnu_options_Callback(hObject, eventdata, handles)
button_group_option_Callback(hObject, eventdata, handles)

function mnu_result_viewer_Callback(hObject, eventdata, handles)
button_group_result_Callback(hObject, eventdata, handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over button_group_parameters.
function button_group_parameters_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to button_group_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_help_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when Main is resized.
function Main_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to Main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
