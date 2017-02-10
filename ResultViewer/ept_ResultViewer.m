%% Result Viewer for ept_TFCE produced result files
% Copyright(C) 2012  Armand Mensen

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Revision History
%
% Version History now tracked in GitHub...
%
% Version 6.1
% 04.10.2013
% - Fixed "ismatrix" call for older versions of matlab
% - Created own ept_fminsearch to avoid the EEGLAB fminsearch overwrite
% - Fixed Factor Selection issues
% 
% Version 6.0
% - Added support for single condition (Ttest) channel*frequency*time data
% - Complete reworked the cluster statistics (and underlying C file) which is now displayed on a separate figure
% - Added an individual topoplots feature to see how consistant the differences are for each person (ie the differences are caused by only one or two participants)
% - Added all analysis and data information to a menu button
% 
%
% Version 5.3
% 15.08.2013
% Added Bar Chart Option and removed EEGLAB topoplot export
% Changed Cluster Table Layout and added the F-Value at peak
% 
% Version 5.21
% 15.02.2013
% Changed the Standard Error calculation of within-group effects to Cousineau method to eliminate individual differences
%  
% Version 5.2
% 08.02.2013
% Added Interaction Plots (no standard error bars, too crowded)
% Changed legend coding to avoid warning
% 
% Version 5.11
% 01.01.2013
% Added TFCE.t support for all factor plots
% 
% Version 5.1
% 12.12.2012
% New button to export the topoplot
% New button which automatically generates original ERP Topoplots + The TFCE-Stat Topoplot (exportable to svg format for publication)
% Cleaned up some random segments of code
% 
% Version 5
% 05.12.2012
% New drop-down menu to select the topoplot colormap
% Ability to set the baseline time in the menu
% New ERP-Plotting using area patches for standard error
% 
% 04.12.2012
% Exports standard error of ERPs instead of deviation
% 
% Version 4
% 16.10.2012
% Updated to streamlined output format
% Generalised ERP Plot making to any number of levels for both factors
% using loops (speed not tested, but shouldn't be a huge issue here
%       Now also easier to change the plot settings
% Included 'publication-ready' correlation/ERP graphs
% 
% 
% 14.10.2012
% Reads the new TFCE output format (Info/Result structure)
%
% 09.01.2011
% Included if statements for different analysis types TFCE or CSize

% Notes
% a=cell2mat(cellfun(@(x) x(:,53,169), Data, 'UniformOutput', false)); %(use to extract data from a particular channel sample for further analysis (e.g. separate ANCOVA).
% Hum = squeeze(mean(cell2maResults_10-Apr-2015t(Data(:,2)),1)); %Use this to extract average of certain group or condition 

%%
function varargout = ept_ResultViewer(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ept_ResultViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @ept_ResultViewer_OutputFcn, ...
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

% --- Executes just before ept_ResultViewer is made visible.
function ept_ResultViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ept_ResultViewer (see VARARGIN)

% Choose default command line output for ept_ResultViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = ept_ResultViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function menu_Load_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
[ResultsFile, ResultsPath] = uigetfile('*.mat', 'Please Select the Results File');
FullFileName = strcat(ResultsPath, ResultsFile);
load (FullFileName);

%% Check whether �frequency sampling rates are specified
if ndims(Data{1}) == 4 && isfield(Info.Parameters, 'fSample') == 0
    
    fSample = inputdlg('Frequency Sampling Rate', 'User Input', 1, {'1'});
    Info.Parameters.fSample = str2num(fSample{1});
    
    % Save the Results File with the new parameter set
    save(strcat(ResultsPath, ResultsFile), 'Info', '-append');
    
elseif isfield(Info.Parameters, 'fSample') == 1
    
else
    Info.Parameters.fSample = [];
    
    % Save the Results File with the new parameter set
    save(strcat(ResultsPath, ResultsFile), 'Info', '-append');
    
end

%% Create handles to the data
handles.Info        = Info;
handles.Data        = Data;
handles.Results     = Results;

% Assign default of factor A which is later changed by the pop menu
% Check to see if this is an ANOVA or single comparision (T, or correlation)
if isstruct(Results.TFCE_Obs)
    handles.Obs         = Results.TFCE_Obs.A;
    handles.nP_Values   = Results.P_Values.A;
    handles.F_Values    = Results.Obs.A;
else
    handles.Obs         = Results.TFCE_Obs;
    handles.nP_Values   = Results.P_Values;
    handles.F_Values    = Results.Obs;
end

%% Check arguments
% Make the pop-up menu invisible if its a single comparison analysis...
if ~isstruct(Results.P_Values)
    set(handles.Pop_FactorSelect,'Value', 1);
    % make popup invisible
    set(handles.Pop_FactorSelect, 'Visible', 'off')
end

%% Plot the Time vs Significance Graph
set(handles.Time_Sig_Plot, 'XTickMode', 'auto');
set(handles.Time_Sig_Plot, 'YTickMode', 'manual');% Reset X axis values
set(handles.Time_Sig_Plot, 'xlim',      [0 Info.Parameters.nSamples]);

if isstruct(Results.P_Values); 
    plot(handles.Time_Sig_Plot, sum(-log(Results.P_Values.A)), 'k', 'LineWidth', 2, 'HitTest', 'off'); 
    
    %Also display main findings...
    %FactorA
    [min_P, idx] = min(Results.P_Values.A(:));
    [Ch, S]      = ind2sub(size(Results.P_Values.A),idx);
    max_Obs      = Results.Obs.A(idx);

    display(['FactorA peak significance found at channel ', num2str(Ch), ' at sample ', num2str(S), ': F(', num2str(size(Data,1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);

    % FactorB
    [min_P, idx] = min(Results.P_Values.B(:));
    [Ch, S]      = ind2sub(size(Results.P_Values.B),idx);
    max_Obs      = Results.Obs.B(idx);

    display(['FactorB peak significance found at channel ', num2str(Ch), ' at sample ', num2str(S), ': F(', num2str(size(Data,2)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);

    % Interaction
    [min_P, idx] = min(Results.P_Values.AB(:));
    [Ch, S]      = ind2sub(size(Results.P_Values.AB),idx);
    max_Obs      = Results.Obs.AB(idx);

    display(['Interaction peak significance found at channel ', num2str(Ch), ' at sample ', num2str(S), ': F(', num2str(size(Data,1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);
    
else
    
    if ndims(Data{1}) == 3
        plot(handles.Time_Sig_Plot, sum(-log(Results.P_Values)), 'k', 'LineWidth', 2, 'HitTest', 'off'); 

        % Also display main findings...
        [min_P, idx] = min(Results.P_Values(:));
        [Ch, S]      = ind2sub(size(Results.P_Values),idx);
        max_Obs      = Results.Obs(idx);

        display(['Peak significance found at channel ', num2str(Ch), ', at sample ', num2str(S), ': T(', num2str(size(Data{1},1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);
    else
        contourf(handles.Time_Sig_Plot, squeeze(sum(-log(Results.P_Values))), 'hittest', 'off', 'linestyle', 'none'); 

        set(handles.Time_Sig_Plot,...
            'NextPlot',     'add'                               ,...
            'YTickMode',    'auto'                              ,...
            'YAxisLocation','right'                             ,...
            'xlim',         [1,size(Results.P_Values,3)]        ,...
            'ylim',         [1,size(Results.P_Values,2)]        );
        
        % Also display main findings...
        [min_P, idx] = min(Results.P_Values(:));
        [Ch, F, S]   = ind2sub(size(Results.P_Values),idx);
        max_Obs      = Results.Obs(idx);

        display(['Peak significance found at channel ', num2str(Ch), ', at frequency bin ', num2str(F), ', at sample ', num2str(S), ': T(', num2str(size(Data{1},1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);
    end
end

set(handles.txt_Name,'String',ResultsFile);

% the colormap to matlab default
colormap(parula);

%% Prepare the interpolation parameters for fast topoplot
PrepareTopoplot(hObject, eventdata, handles)

function menu_FileInfo_Callback(hObject, eventdata, handles)
% hObject    handle to menu_FileInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ndims(handles.Data{1}) == 4
    FileInfo = ['Data Info... \n\n',...
            'Number of Channels: ',   num2str(size(handles.Data{1},2)),             '\n'...
            'Number of Samples: ',    num2str(size(handles.Data{1},4)),             '\n'...
            'Time Sampling Rate: ',   num2str(handles.Info.Parameters.rSample),     '\n'...
            'Number of Freq. Bins: ', num2str(size(handles.Data{1},3)),             '\n'...
            'Freq. Sample Rate: ',    num2str(handles.Info.Parameters.fSample),     '\n'...
            '\nStatistical Info... \n\n',...
            handles.Info.Comments,                                                  '\n'...
            'Statistical Model: ', handles.Info.Parameters.type,                    '\n'...
            'Number of Permutations: ', num2str(handles.Info.Parameters.nPerm),     '\n'...
            'E: ', num2str(handles.Info.Parameters.E_H(1)), ' | H: ', num2str(handles.Info.Parameters.E_H(2)), '\n'...
            ];
else
    FileInfo = ['Data Info... \n\n',...
            'Number of Channels: ',   num2str(size(handles.Data{1},2)),             '\n'...
            'Number of Samples: ',    num2str(size(handles.Data{1},3)),             '\n'...
            'Time Sampling Rate: ',   num2str(handles.Info.Parameters.rSample),     '\n'...
            '\nStatistical Info... \n\n',...
            handles.Info.Comments,                                                  '\n'...
            'Statistical Model: ', handles.Info.Parameters.type,                    '\n'...
            'Number of Permutations: ', num2str(handles.Info.Parameters.nPerm),     '\n'...
            'E: ', num2str(handles.Info.Parameters.E_H(1)), ' | H: ', num2str(handles.Info.Parameters.E_H(2)), '\n'...
            ];    
end

handles.MsgBox = msgbox(sprintf(FileInfo), 'File Info');

%Refresh GUI
guidata(hObject,handles);


% --- Executes on selection change in Pop_FactorSelect.
function Pop_FactorSelect_Callback(hObject, eventdata, handles)

Factor = get(handles.Pop_FactorSelect,'Value');

switch Factor
    case 1
        % Check to see if this is an ANOVA or single comparision (T, or correlation)
        if isstruct(handles.Results.TFCE_Obs)
            handles.Obs         = handles.Results.TFCE_Obs.A;
            handles.nP_Values   = handles.Results.P_Values.A;
            handles.F_Values    = handles.Results.Obs.A;
        else
            handles.Obs         = handles.Results.TFCE_Obs;
            handles.nP_Values   = handles.Results.P_Values;
        end
    case 2
        handles.Obs = handles.Results.TFCE_Obs.B;
        handles.nP_Values = handles.Results.P_Values.B;
        handles.F_Values  =  handles.Results.Obs.B;
    case 3
        handles.Obs = handles.Results.TFCE_Obs.AB;
        handles.nP_Values = handles.Results.P_Values.AB;
        handles.F_Values  = handles.Results.Obs.AB;
end

update_Time_Sig_Plot(hObject, handles)
makeTopoplot(hObject, eventdata, handles)

% Redraw Sigplot with line

%Refresh GUI
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Pop_FactorSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pop_FactorSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on mouse press over axes background.
function Time_Sig_Plot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Time_Sig_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Axes property 'NextPlot' set to 'Add' so that properties are not reset

Curr_Point  = get(handles.Time_Sig_Plot, 'CurrentPoint');
rSample     = handles.Info.Parameters.rSample;
fSample     = handles.Info.Parameters.fSample;

sPoint      = round(Curr_Point(1,1));
fPoint      = round(Curr_Point(1,2));

bL          = str2double(get(handles.txt_Baseline,'String')); %get the specified baseline

%% Set new parameters
set(handles.tx_sPoint,'String',['S ' (num2str(sPoint))]);

if ndims(handles.Data{1}) == 4
    set(handles.tx_fPoint,'String',['F ' (num2str(fPoint))]);
    set(handles.tx_Freq,'String',[num2str(round(fPoint/fSample)) ' Hz']);
end

T_P = round(sPoint/rSample*1000)-bL;
if T_P >= 2000
    set(handles.tx_Time,'String',[num2str(T_P/1000) ' s']);
else
    set(handles.tx_Time,'String',[num2str(T_P) ' ms']);
end

makeTopoplot(hObject, eventdata, handles)

guidata(hObject,handles);

% --- Executes on button press in Left_Button.
function Left_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Left_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rSample = handles.Info.Parameters.rSample;
bL      = str2double(get(handles.txt_Baseline,'String')); %get the specified baseline

sPoint = sscanf(get(handles.tx_sPoint,'String'),'%*s %i');
fPoint  = sscanf(get(handles.tx_fPoint,'String'),'%*s %i');

set(handles.tx_sPoint,'String',['S ' (num2str(sPoint - 1))]);

T_P = round((sPoint-1)/rSample*1000)-bL;
if T_P >= 2000
    set(handles.tx_Time,'String',[num2str(T_P/1000) ' s']);
else
    set(handles.tx_Time,'String',[num2str(T_P) ' ms']);
end

% Reset the Significance Plot and Cursor Position
update_Time_Sig_Plot(hObject, handles)

makeTopoplot(hObject, eventdata, handles)

guidata(hObject,handles);

% --- Executes on button press in Right_Button.
function Right_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Right_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rSample = handles.Info.Parameters.rSample;
bL      = str2double(get(handles.txt_Baseline,'String')); %get the specified baseline

sPoint = sscanf(get(handles.tx_sPoint,'String'),'%*s %i');
fPoint  = sscanf(get(handles.tx_fPoint,'String'),'%*s %i');

set(handles.tx_sPoint,'String',['S ' (num2str(sPoint + 1))]);

T_P = round((sPoint+1)/rSample*1000)-bL;
if T_P >= 2000
    set(handles.tx_Time,'String',[num2str(T_P/1000) ' s']);
else
    set(handles.tx_Time,'String',[num2str(T_P) ' ms']);
end

% Reset the Significance Plot and Cursor Position
update_Time_Sig_Plot(hObject, handles)

makeTopoplot(hObject, eventdata, handles)

guidata(hObject,handles);



%% --- Executes on Button Presses for Extra Tools
function pb_FactorLevels_Callback(hObject, eventdata, handles)

% Get necessary data from the handles
Data        = handles.Data;
Results     = handles.Results;
e_loc       = handles.Info.Electrodes.e_loc;

% Get samples and factor from GUI
sPoint  = sscanf(get(handles.tx_sPoint,'String'),'%*s %i');
fPoint  = sscanf(get(handles.tx_fPoint,'String'),'%*s %i');
Factor  = get(handles.Pop_FactorSelect,'Value');

if ndims(Data{1}) == 4 % For Time Frequency...
    ept_FactorTopoplots(Data, e_loc, Results, Factor, fPoint+size(Data{1},3)*(sPoint-1));
else % For 2D data...
    ept_FactorTopoplots(Data, e_loc, Results, Factor, sPoint);
end


guidata(hObject,handles);

%% --- Executes on button press in pb_IndividualTopoplots.
function pb_IndividualTopoplots_Callback(hObject, eventdata, handles)

% Get necessary data from the handles
Data        = handles.Data;
Results     = handles.Results;
e_loc       = handles.Info.Electrodes.e_loc;

% Get samples and factor from GUI
sPoint  = sscanf(get(handles.tx_sPoint,'String'),'%*s %i');
fPoint  = sscanf(get(handles.tx_fPoint,'String'),'%*s %i');
Factor  = get(handles.Pop_FactorSelect,'Value');

% check for frequency domain and adjust if necessary
if fPoint > 1
    ept_IndividualTopoplots(Data, e_loc, Results, Factor, fPoint+size(Data{1},3)*(sPoint-1));
else
    ept_IndividualTopoplots(Data, e_loc, Results, Factor, sPoint);
end


%% --- Executes on button press in pb_BarChart.
function pb_BarChart_Callback(hObject, eventdata, handles)

% thresh  = str2double(get(handles.Txt_Thresh,'String'));
Data        = handles.Data;
e_loc       = handles.Info.Electrodes.e_loc;

Sample = sscanf(get(handles.tx_sPoint,'String'),'%*s %i');
Freq   = sscanf(get(handles.tx_fPoint,'String'),'%*s %i');

% Get Desired Channel
Channel     = inputdlg('Create Bar Chart', 'Channel Name', 1, {'E129'});
Channel     = Channel{1};
ChannelId   = strcmp(Channel,{e_loc(:).labels});


% Extract Data
if ndims(Data{1})==3
    BarData = cell2mat(cellfun(@(x) mean(x(:,ChannelId,Sample)), Data, 'UniformOutput', false));
    BarSE   = cell2mat(cellfun(@(x) std(x(:,ChannelId,Sample))/size(x,1), Data, 'UniformOutput', false));
else
    BarData = cell2mat(cellfun(@(x) mean(x(:,ChannelId,Freq,Sample)), Data, 'UniformOutput', false));
    BarSE   = cell2mat(cellfun(@(x) std(x(:,ChannelId,Freq,Sample))/size(x,1), Data, 'UniformOutput', false));
end

% Plot Bar Chart
hBC.Figure = figure;
set(hBC.Figure,...
    'Color',            'w'                 ,...
    'Renderer',         'openGL'            );

hBC.Axes = axes;
hBC.Bar  = bar(BarData);

set(hBC.Axes,...
    'NextPlot',         'add'               );

% Plot Error Bars
hBC.Error = zeros(1,size(BarData,2));
for i = 1:size(BarData,2)
    x = get(get(hBC.Bar(i),'children'),'xdata'); % X-Axis locations of that level
    hBC.Error(i) = errorbar(mean(x,1), BarData(:,i), BarSE(:,i), BarSE(:, i), '.k');
end


% Make it look pretty...
if ndims(Data{1})==3
    hBC.Title  = title (['ERP Bar Chart (', Channel, ' at Sample ' num2str(Sample), ')']  );
else
    hBC.Title  = title (['ERP Bar Chart (', Channel, ' at Sample ' num2str(Sample), ' at Frequency Bin ' num2str(Sample) ')']  );
end

hBC.xLabel = xlabel('FactorA Levels');
hBC.ylabel = ylabel('ERP Strength (MicroVolts)');

set(hBC.Axes,...
    'Box',              'off'                   ,...
    'XTickLabel',       ['Level 1'; 'Level 2'; 'Level 3']  );

set([hBC.Title, hBC.xLabel, hBC.ylabel], ...
    'FontName'   , 'AvantGarde'     );
set(hBC.Title, ...
    'FontSize'   , 12               , ...
    'FontWeight' , 'bold'           );

guidata(hObject,handles);

%% --- Executes on button press in Export Button.
function pb_ExportTopoplot_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ExportTopoplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

e_loc       = handles.Info.Electrodes.e_loc;
sPoint         = sscanf(get(handles.tx_sPoint,'String'),'%*s %i');
V           = handles.Obs(:,sPoint);

ept_Topoplot(V, e_loc,...
    'NewFigure', 1, ...
    'PlotChannels', 0);

guidata(hObject,handles);


%% --- Executes on button press in pb_CalculateCluster.
function pb_CalculateCluster_Callback(hObject, eventdata, handles)

%% -- Get necessary parameters -- %%
P_Values    = handles.nP_Values;
TFCE_Obs    = handles.Obs;
Obs_Values  = handles.F_Values;
ChN         = handles.Info.Electrodes.ChannelNeighbours;

thresh      = str2double(get(handles.Txt_Thresh,'String'));

%% -- Process Data -- %%
% For 2D datasets (Channel*Time || Channel*Frequency
if ndims(handles.Data{1}) == 3
    
    % Look for new clusters
    x = TFCE_Obs; % copy TFCE_Obs
    x(P_Values>thresh) = 0; % Threshold the data at alpha
    Cp = ept_ClusRes(x, ChN, 0.01); % Calculate Positive Clusters
    
    xn = x; % copy x
    xn(x>0)=0; % Only show negative x
    xn = abs(xn); % Make negative x's positive
    Cn = ept_ClusRes(xn, ChN, 0.01); % Calculate Negative Clusters
    
    C = Cn-Cp; % Combine the two
    
    % export the cluster data to the global workspace to use channel indices
    fprintf(1, 'Info: Channel cluster results assigned to workspace\n');
    assignin('base', 'current_clusters', C);
    
    b = unique(C); % How many different clusters are there?
    b(b==0)=[]; % Eliminate the 0 from being a unique cluster
    
    if numel(b)==0;
        display (['There are no clusters of significant data at the p = ' num2str(thresh) ' threshold']);
        return
    else
        
        ClusRes= cell(size(b,1),11);
        
        for i = 1:size(b,1);
            x = TFCE_Obs;
            x(C~=b(i))      = 0; %
            idPeak          = find(abs(x)==max(abs(x(:))));
            [PeakC, PeakS]  = ind2sub(size(x),idPeak);
            
            idSize          = find(C== b(i)); % find the rows and columns that are significant
            [SizeC, SizeS]  = ind2sub(size(C),idSize);
            
            ClusRes{i,1}    = PeakC(1); % peak channel (just the first of many possible peak channels (but averaging may result in a channel in between two that is not significant)!
            ClusRes{i,2}    = PeakS(1);
            ClusRes{i,3}    = [];
            ClusRes{i,4}    = Obs_Values(PeakC(1),PeakS(1));
            ClusRes{i,5}    = P_Values(PeakC(1),PeakS(1));
            ClusRes{i,6}    = numel(idSize);
            ClusRes{i,7}    = numel(unique(SizeC));
            ClusRes{i,8}    = numel(unique(SizeS));
            ClusRes{i,9}    = [];
            ClusRes{i,10}   = [num2str(min(SizeS)), ' - ', num2str(max(SizeS))];
            ClusRes{i,11}   = [];
            
        end
    end
    
else % for Time-Freqency Data...
    
    % Calculate clusters above p-value threshold
    data = TFCE_Obs; % copy TFCE_Obs
    data(P_Values>thresh) = 0; % Threshold the data at alpha
    C = ept_ClusRes3D(data, ChN, 0.01); % Calculate Negative Clusters    
    
    b = unique(C); % How many different clusters are there?
    b(b==0)=[]; % Eliminate the 0 from being a unique cluster
    
    if numel(b)==0;
        display (['There are no clusters of significant data at the p = ' num2str(thresh) ' threshold']);
        return
    end
    
    for i = 1:size(b,1);  
        x = TFCE_Obs;
        x(C~=b(i))      = 0; %
        idPeak          = find(abs(x)==max(abs(x(:))));
        [PeakC, PeakF, PeakS] = ind2sub(size(x),idPeak);
        
        idSize          = find(C== b(i)); % find the rows and columns that are significant
        [SizeC, SizeF, SizeS]=ind2sub(size(C),idSize);

        ClusRes{i,1}    = PeakC(1); % peak channel (just the first of many possible peak channels (but averaging may result in a channel in between two that is not significant)!
        ClusRes{i,2}    = PeakS(1);
        ClusRes{i,3}    = PeakF(1);
        ClusRes{i,4}    = Obs_Values(PeakC(1),PeakF(1),PeakS(1));
        ClusRes{i,5}    = P_Values(PeakC(1),PeakF(1),PeakS(1));                
        ClusRes{i,6}    = numel(idSize);
        ClusRes{i,7}    = numel(unique(SizeC));
        ClusRes{i,8}    = numel(unique(SizeS));
        ClusRes{i,9}    = numel(unique(SizeF));
        ClusRes{i,10}   = [num2str(min(SizeS)), ' - ', num2str(max(SizeS))];
        ClusRes{i,11}   = [num2str(min(SizeF)), ' - ', num2str(max(SizeF))];
        
    end    
end     

%% Should make figure positions etc based on pixels because Java is...


%% -- Create external figure -- %%
handles.Fig.ClusterRes.F = figure(...   
    'Name',         'Statistics for Significant Clusters',...
    'Units',        'normalized'            ,...
    'Position',     [0.5 0.6 0.5 0.35]      ,...
    'Color',        'w'                     ,...
    'Toolbar',      'None'                  );

handles.Fig.ClusterRes.T = uitable(...
    'Parent',       handles.Fig.ClusterRes.F,...
    'Units',        'normalized'            ,...
    'Position',     [0.025 0.025 0.95 0.8]  ,...
    'FontSize',     12                      );

set(handles.Fig.ClusterRes.T, 'ColumnName' ,...
    {...
    '<html><center>Peak<br />Channel</center></html>'       ,...
    '<html><center>Peak<br />Sample</center></html>'        ,...
    '<html><center>Peak<br />Frequency</center></html>'     ,...
    '<html><center>Statistic<br />At Peak</center></html>'  ,...
    '<html><center>P-Value<br />At Peak</center></html>'    ,...
    '<html><center>Cluster<br />Size</center></html>'       ,...
    '<html><center>Unique<br />Channels</center></html>'    ,...
    '<html><center>Unique<br />Samples</center></html>'     ,...
    '<html><center>Unique<br />Frequencies</center></html>' ,...
    '<html><center>Sample<br />Range</center></html>'       ,...
    '<html><center>Frequency<br />Range</center></html>'    });

% Get Java handle.
jscroll = findjobj(handles.Fig.ClusterRes.T);
jtable = jscroll.getViewport.getView;    

% Some table settings...
jtable.setRowResizable(true);       % Make row height user resizable...
jtable.setColumnResizable(true);    % Make column width user resizable...
jtable.setRowHeight(30);            % Set specific row height in pixels...

% Make columns sortable...
jtable.setSortable(true);
jtable.setAutoResort(true);
jtable.setMultiColumnSortable(true);
jtable.setPreserveSelectionsAfterSorting(true);

% Make the vertical scroll permananent
% jscroll.java.VERTICAL_SCROLLBAR_ALWAYS;
set(jscroll,'VerticalScrollBarPolicy', 22);

% Set the callback now that we have access to the jtable...
set(handles.Fig.ClusterRes.T, 'CellSelectionCallback', {@Callback_TableSelection, handles, jtable});

% Text box to indicate significant threshold
handles.Fig.ClusterRes.Title = uicontrol(...
    'Parent',       handles.Fig.ClusterRes.F,...
    'Style',        'text'                  ,...
    'Units',        'normalized'            ,...
    'Position',     [0.025 0.85 0.95 0.1]    ,...
    'FontSize',     16                      ,...
    'String',       ['Cluster(s) of Significant Points at Threshold ' num2str(thresh)]);

jtitle = findjobj(handles.Fig.ClusterRes.Title);
jtitle.setVerticalAlignment(0);

%% Put the data in there
set(handles.Fig.ClusterRes.T, 'Data', ClusRes);

guidata(hObject, handles);

%% -- Called when cell is selected in the table --- %%
function Callback_TableSelection(mtable, eventdata, handles, jtable);

Data     = get(mtable, 'data');

if isempty(eventdata.Indices);return;end

RowIdx   = eventdata.Indices(1);
ColIdx   = eventdata.Indices(2);

SelectedValue = Data{RowIdx, ColIdx};

if ColIdx == 1 % Channel selected
    % Highlight the specified channel on the topoplot (if possible)...
    x = get(handles.HeadView, 'children'); %all handles in the HeadView axes (first two are the surface and countours, rest are channels)
    x = x(numel(x)-SelectedValue);
    
    tmpstr = get(x, 'userdata') ;
    set(x, 'userdata', get(x, 'string'));
    set(x, 'string', tmpstr);
    
elseif ColIdx == 2 % Sample selected
    
    set(handles.tx_sPoint,'String',['S ' (num2str(SelectedValue))]);
    
    % Reset the Significance Plot and Cursor Position
    update_Time_Sig_Plot(mtable, handles)
    makeTopoplot(mtable, eventdata, handles)
    
    % Update time point
    T_P = round((SelectedValue)/handles.Info.Parameters.rSample*1000)-str2double(get(handles.txt_Baseline,'String'));
    if T_P >= 2000
        set(handles.tx_Time,'String',[num2str(T_P/1000) ' s']);
    else
        set(handles.tx_Time,'String',[num2str(T_P) ' ms']);
    end
    
elseif ColIdx == 3 % Frequency selected

    set(handles.tx_fPoint,'String',['F ' (num2str(SelectedValue))]);
    set(handles.tx_Freq,'String',[num2str(round(SelectedValue/handles.Info.Parameters.fSample)) ' Hz']);

    
    % Reset the Significance Plot and Cursor Position
    update_Time_Sig_Plot(mtable, handles)
    makeTopoplot(mtable, eventdata, handles)
    
end

jtable.clearSelection();


%% --- Functions not directly accessible through GUI 
function update_Time_Sig_Plot(hObject, handles)

% Get which factor is selected
x       = get(handles.Pop_FactorSelect,'Value');
sPoint  = sscanf(get(handles.tx_sPoint,'String'),'%*s %i');
fPoint  = sscanf(get(handles.tx_fPoint,'String'),'%*s %i');

cla(handles.Time_Sig_Plot)
cla(handles.HeadView)
switch x
    case 1
        % Check to see if this is an ANOVA or single comparision (T, or correlation)
        if isstruct(handles.Results.TFCE_Obs)
            handles.Obs         = handles.Results.TFCE_Obs.A;
            handles.nP_Values   = handles.Results.P_Values.A;
            handles.F_Values    = handles.Results.Obs.A;
        else
            handles.Obs         = handles.Results.TFCE_Obs;
            handles.nP_Values   = handles.Results.P_Values;
        end
    case 2
        handles.Obs = handles.Results.TFCE_Obs.B;
        handles.nP_Values = handles.Results.P_Values.B;
        handles.F_Values  =  handles.Results.Obs.B;
    case 3
        handles.Obs = handles.Results.TFCE_Obs.AB;
        handles.nP_Values = handles.Results.P_Values.AB;
        handles.F_Values  = handles.Results.Obs.AB;
end

if ndims(handles.Data{1}) == 3
    area(handles.Time_Sig_Plot, sum(-log(handles.nP_Values)), 'FaceColor', [0.5 0.5 0.5], 'LineWidth', 2, 'HitTest', 'off')
    YV = get(handles.Time_Sig_Plot, 'YLim');
    axes(handles.Time_Sig_Plot)
    line([sPoint, sPoint], YV, 'color', [0.7, 0.7, 0.7], 'LineWidth', 3, 'Parent', handles.Time_Sig_Plot);
else
    contourf(handles.Time_Sig_Plot, squeeze(sum(-log(handles.nP_Values))), 'hittest', 'off', 'linestyle', 'none');
    axes(handles.Time_Sig_Plot)
    text(sPoint, fPoint, '+',...
        'HorizontalAlignment',  'center'    ,...
        'VerticalAlignment',    'middle'    ,...
        'FontSize',             20          ,...
        'FontWeight',           'bold'      );
end

%Refresh GUI
guidata(hObject,handles);

function PrepareTopoplot(hObject, eventdata, handles)
% hObject    handle to plot_axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles = get(gcf,  'userData');

e_loc = handles.Info.Electrodes.e_loc;
Obs   = handles.Obs;

%% Prepare the interpolation function for ept_topoplot

Th = pi/180*[e_loc.theta];        % Calculate theta values from x,y,z e_loc
Rd = [e_loc.radius];              % Calculate radian values from x,y,z e_loc

x = Rd.*cos(Th);                            % Calculate 2D projected X
y = Rd.*sin(Th);                            % Calculate 2D projected Y

% Squeeze the coordinates into a -0.5 to 0.5 box
intrad = min(1.0,max(abs(Rd))); intrad = max(intrad,0.5); squeezefac = 0.5/intrad;

snx = x'*squeezefac; sny = y'*squeezefac;

gridscale = 100;
Xq = linspace(-0.5,0.5,gridscale);
XYq = Xq(ones(gridscale,1),:);

% Test for TF data type...
if ndims(Obs) == 3
    F  = TriScatteredInterp(snx, sny, Obs(:,1,1), 'natural');
else
    F  = TriScatteredInterp(snx, sny, Obs(:,1), 'natural');
end

%% Set the handles for the topoplot

handles.F   = F;
handles.XYq = XYq;
handles.snx = snx;
handles.sny = sny;

makeTopoplot(hObject, eventdata, handles)

guidata(hObject,handles);


function makeTopoplot(hObject, eventdata, handles)

e_loc       = handles.Info.Electrodes.e_loc;
XYq         = handles.XYq;
gridscale   = 100;           % Topoplot Details (100 is high)

%% Calculate
thresh  = str2double(get(handles.Txt_Thresh,'String'));
sPoint = sscanf(get(handles.tx_sPoint,'String'),'%*s %i');
fPoint = sscanf(get(handles.tx_fPoint,'String'),'%*s %i');

% % Redraw Sigplot with line
update_Time_Sig_Plot(hObject, handles);

if ndims(handles.Data{1})==3
    data = handles.Obs(:, sPoint);
    id{1}    = find(handles.nP_Values(:,sPoint)<thresh); % Significant Channels
    id{2}    = find(handles.nP_Values(:,sPoint)>=thresh); % Non-significant Channels
else
    data = handles.Obs(:,fPoint, sPoint); 
    id{1}    = find(handles.nP_Values(:,fPoint,sPoint)<thresh); % Significant Channels
    id{2}    = find(handles.nP_Values(:,fPoint,sPoint)>=thresh); % Non-significant Channels
end

handles.F.V = data;
Zi = handles.F(XYq', XYq);
unsh = (gridscale+1)/gridscale; % un-shrink the effects of 'interp' SHADING

% Create the Topoplot
cla(handles.HeadView, 'reset')

set(handles.HeadView                        , ...
    'XLim'                , [-0.525 0.525]  , ...
    'YLim'                , [-0.525 0.525]  ,...
    'NextPlot'            , 'add'           ,... % Hold on
    'CameraPosition'      , [0 0 1]         );   % Since surface is a 3D action, change camera to seem 2D
axis(handles.HeadView, 'square', 'off')    

% hSurf       = surf(handles.HeadView, XYq*unsh, XYq'*unsh, zeros(size(Zi)), Zi,...
%     'EdgeColor','none',...
%     'FaceColor','interp');
hContour    = contourf(handles.HeadView, XYq, XYq', Zi*.075, 8,'k');

%Plot channels
labels    = {e_loc.labels};  

axes(handles.HeadView)
for i = 1:size(labels,2)
  hChannels(i) = text(handles.sny(i),handles.snx(i), '.'  ,...
      'userdata'          , char(labels(i)) );
end

if ndims(handles.Data{1})==3
    set(hChannels                               ,...
        'HorizontalAlignment',  'center'        ,...
        'VerticalAlignment',    'middle'        ,...
        'Color',                'k'             ,...
        'FontSize',             10              ,...
        'FontWeight',           'bold'          ,...
        'buttondownfcn',        {@makeERPplot, handles});
    set(hChannels(id{1})                        ,...
        'String',               '+'             );
else
    set(hChannels                               ,...
        'HorizontalAlignment',  'center'        ,...
        'VerticalAlignment',    'middle'        ,...
        'Color',                'k'             ,...
        'FontSize',             10              ,...
        'FontWeight',           'bold'          ,...
        'buttondownfcn',        {@makeWaveletPlot, handles});
    set(hChannels(id{1})                        ,...
        'String',               '+'             );
end

     
guidata(hObject,handles);
 
function makeERPplot(hObject, eventdata, handles)
% hObject    handle to plot_axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

seltype = get(gcf,  'selectiontype');
Ch      = get(gcbo, 'userData');
FactorType = get(handles.Pop_FactorSelect,'Value');

% check whether the name of the electrode or position marker was selected
if any(strcmp(Ch, {'.', '+'}))
   Ch      = get(gcbo, 'string');    
end

switch seltype

    case 'alt'
        thresh      = str2double(get(handles.Txt_Thresh,'String'));
        bL          = str2double(get(handles.txt_Baseline,'String')); %get the specified baseline

        Data        = handles.Data;
        e_loc       = handles.Info.Electrodes.e_loc;
        P_Values    = handles.nP_Values;

        e_labels    = {e_loc(:).labels};
        IdCh        = strcmp(Ch,e_labels);

        if handles.Info.Parameters.type == 'c'
            % Draw the appropriate scatterplot for this channel and sample
            
            sPoint = sscanf(get(handles.tx_sPoint,'String'),'%*s %i'); %Get the current sample
                        
            y = Data{1}(:,IdCh,sPoint);
            x = Data{2};
            
            p = polyfit(x,y,1);   % p returns 2 coefficients fitting r = a_1 * x + a_2
            r = p(1) .* x + p(2); % compute a new vector r that has matching datapoints in x
            
            %% Creating the actual figure; (http://blogs.mathworks.com/loren/2007/12/11/making-pretty-graphs/)
            hScatter = figure('Units', 'pixels', 'Position', [100 100 500 375], 'Color', 'w');
            set(hScatter,'name','ERP and Behavioural Correlation','numbertitle','off')
            hold on;
            
            % Put the raw data
            hData = line(x,y);  % data
            hLine = line(x,r);  % trendline
            
            % Edit the points
            set(hData                           , ...
                'LineStyle'       , 'none'      , ...
                'Marker'          , '.'         , ...
                'Marker'          , 'o'         , ...
                'MarkerSize'      , 6           , ...
                'MarkerEdgeColor' , 'none'      , ...
                'MarkerFaceColor' , [.25 .25 1] );
            
            set(hLine                           , ...
                'LineStyle'       , '-.'        , ...
                'LineWidth'       , 1.5         , ...
                'Color'           , [0 .5 0]    );
            
            % Add Legends and Labels 
            
            hTitle  = title (['ERP (', Ch, ' at Sample ' num2str(sPoint), ') & Behavioural Correlation']  );
            hXLabel = xlabel('Behavioural Measure'            );
            hYLabel = ylabel('ERP Strength (MicroVolts)'      );
                     
            % Get axis limits to adjust text location accordingly
            X_Lim = get(gca, 'XLim'); tx = X_Lim(1)+(X_Lim(2)-X_Lim(1))*.1; % 10% from the left
            Y_Lim = get(gca, 'YLim'); ty = Y_Lim(2)-(Y_Lim(2)-Y_Lim(1))*.1; % 10% from the top
            
            hText   = text(tx, ty, '');
            
            set(hText, ...
                'string', ['\itr = ' num2str(handles.Results.Obs(IdCh, sPoint), '%0.2f') ', p = ' num2str(P_Values(IdCh, sPoint), '%0.3f')]);
            
            % Adjust Font and Axes Properties
            
            set([hTitle, hXLabel, hYLabel, hText], ...
                'FontName'   , 'AvantGarde'     );
            set(hTitle                          , ...
                'FontSize'   , 12               , ...
                'FontWeight' , 'bold'           );
            set(hText                           ,...
                'EdgeColor'  , [0 0 0]          ,...
                'BackgroundColor',[.7 .9 .7]    );
            
            set(gca, ...
                'Box'         , 'off'     , ...
                'TickDir'     , 'out'     , ...
                'TickLength'  , [.02 .02] , ...
                'XMinorTick'  , 'on'      , ...
                'YMinorTick'  , 'on'      , ...
                'YGrid'       , 'on'      , ...
                'XColor'      , [.3 .3 .3], ...
                'YColor'      , [.3 .3 .3], ...
                'LineWidth'   , .2         );
            

        elseif FactorType == 1 || FactorType == 2 % for normal ERP comparison    
            switch FactorType
                case 1
                    
                    if strcmp(handles.Info.Parameters.type, 'i');
                        
                        Data = Data(:);
                        
                        for i = 1:size(Data,1)
                            
                            y(i,:)   = squeeze(mean(cell2mat(cellfun(@(x) x(:,IdCh,:), Data(i,:), 'UniformOutput', false)')));
                            e(i,:)   = squeeze(std(cell2mat(cellfun(@(x) x(:,IdCh,:), Data(i,:), 'UniformOutput', false)')))...
                                / sqrt(handles.Info.Parameters.GroupSizes(1));
                            
                        end
                        
                    elseif strcmp(handles.Info.Parameters.type, 'd');
                        
                        % check for dependent test to 0
                        if max(Data{2}(:)) < 0.05
                            for i = 1:size(Data,1)
                                
                                y(i,:)   = squeeze(mean(cell2mat(cellfun(@(x) x(:,IdCh,:), Data(i, :), 'UniformOutput', false)')));
                                e(i,:)   = squeeze(std(cell2mat(cellfun(@(x) x(:,IdCh,:), Data(i, :), 'UniformOutput', false)')))...
                                    / sqrt(handles.Info.Parameters.GroupSizes(1));
                            end
                            
                        else
                            
                            % make sure Data is in columns
                            Data = reshape(Data, 1, []);
                            
                            % Calculate SE for within group using Cousineau method
                            ParticipantScore = cell2mat(cellfun(@(x) x(:,IdCh,:), Data, 'UniformOutput', false));
                            ParticipantMean  = mean(cell2mat(cellfun(@(x) x(:,IdCh,:), Data, 'UniformOutput', false)),2);
                            ParticipantMeanRep = repmat(ParticipantMean, [1, size(ParticipantScore, 2), 1]);
                            ConditionDifferences = ParticipantScore - ParticipantMeanRep;
                            % single error bar for both (CI = s.e. * 1.96)
                            e = (squeeze(std(ConditionDifferences)) / sqrt(handles.Info.Parameters.GroupSizes(1))) * 1.96;
                            
                            for i = 1:size(Data, 2)
                                y(i,:)   = squeeze(mean(cell2mat(cellfun(@(x) x(:,IdCh,:), Data(:,i), 'UniformOutput', false))));
                            end
                            
                        end
                        
                    else
                        
                        for i = 1:size(Data,1)
                            
                            y(i,:)   = squeeze(mean(cell2mat(cellfun(@(x) x(:,IdCh,:), Data(i,:), 'UniformOutput', false)')));
                            e(i,:)   = squeeze(std(cell2mat(cellfun(@(x) x(:,IdCh,:), Data(i,:), 'UniformOutput', false)')))/sqrt(handles.Info.Parameters.GroupSizes(1));
                            
                        end
                        
                    end
   
                    
                case 2 
                    
%                   Calculate SE for within group using Cousineau method
                    ParticipantScore = cell2mat(cellfun(@(x) x(:,IdCh,:), Data, 'UniformOutput', false));
                    ParticipantMean  = mean(cell2mat(cellfun(@(x) x(:,IdCh,:), Data, 'UniformOutput', false)),2);
                    ParticipantMeanRep = repmat(ParticipantMean, [1,size(ParticipantScore,2),1]);
                    ConditionDifferences = ParticipantScore-ParticipantMeanRep;
                    e = squeeze(std(ConditionDifferences)) / sqrt(handles.Info.Parameters.GroupSizes(1));
                    
                    for i = 1:size(Data, 2)
                    
                        y(i,:)   = squeeze(mean(cell2mat(cellfun(@(x) x(:,IdCh,:), Data(:,i), 'UniformOutput', false))));
                    
                    end

            end 

            time = (1:1:(handles.Info.Parameters.nSamples))*(1000/handles.Info.Parameters.rSample)-bL;

            x = P_Values(IdCh,:);
            x(x>thresh) = 0;

            % Sig_Chan    = abs(P_Values(IdCh,:))/10;
            wSCh        = find(x)*(1000/handles.Info.Parameters.rSample)-bL; % Where Significant Channels (which samples are significant) 
            
            
            nL = size(y,1); % Specify the number of lines in the plot
            x = time; % Rename so I can use the ERPErrorPlot nomenclature :D
            x = x(ones(nL,1),:);
            
            uL = y+e; %Upper Limit
            lL = y-e; %Lower Limit

            yP = [lL, fliplr(uL)];
            xP = [x,fliplr(x)];
             
            %% Actual ERP Plotting

            H.Figure = figure;
            set(H.Figure,...
                'Name',             'ERP Comparison'    ,...
                'NumberTitle',      'off'               ,...
                'Units',            'pixels'            ,...
                'Position',         [100 100 900 500]   ,...
                'Color',            'w'                 ,...
                'Renderer',         'openGL'            ,...
                'DefaultAxesColorOrder',        flipud(hsv(nL))       ,...
                'DefaultAxesLineStyleOrder',    {'-',':', '--', '-.'} );
            H.Axes = axes;
            set(H.Axes,...
                'FontName',         'Helvetica'         ,...
                'XLim',             [x(1,1), x(1,end)]  ,...
                'NextPlot',         'replacechildren'   );
            
            hold on;
            
            H.Line = plot(H.Axes, x',y'); % Create the first line and assign a handle
            set(H.Line,...
                'LineWidth',        3                   );
            
            Y_Lim = get(H.Axes, 'YLim');
            line([0,0],Y_Lim, 'Color', [0 0 0], 'LineStyle', '--')
            
            plot(x',uL',':');
            plot(x',lL',':');

            H.Patch = patch(xP',yP',1);
            set(H.Patch,...
                'EdgeColor',        'none'              ,...
                'FaceColor',        'flat'              ,...
                'FaceVertexCData',  flipud(hsv(nL))     ,...
                'CDataMapping',     'direct'            ,...
                'FaceAlpha',        0.15                );
            
            hTitle  = title (['ERP Comparison for ' Ch]         );
            hXLabel = xlabel('Time (ms)'                        );
            hYLabel = ylabel('ERP Strength (MicroVolts)'        );
            
            % Create Legend Data
            for i = 1:nL
                Legend{i} = ['Level ' int2str(i)];
            end
            hLegend = legend(H.Line, Legend);
            
%             set(H.Axes                       , ...
%                 'FontName'   , 'Helvetica' );
            set([hTitle, hXLabel, hYLabel] , ...
                'FontName'   , 'AvantGarde');
            
            set(hTitle                          , ...
                'FontSize'   , 14               , ...
                'FontWeight' , 'bold'           );
                                   
            set([hXLabel, hYLabel]              , ...
                'FontSize'   , 12               );
            
            set([hLegend, H.Axes]             , ...
                'FontSize'   , 10           );
            
            set(H.Axes, ...
                'Box'         , 'off'     , ...
                'TickDir'     , 'out'     , ...
                'TickLength'  , [.02 .02] , ...
                'XMinorTick'  , 'on'      , ...
                'YMinorTick'  , 'on'      , ...
                'YGrid'       , 'on'      , ...
                'XColor'      , [.3 .3 .3], ...
                'YColor'      , [.3 .3 .3], ...
                'LineWidth'   , .2        );
            
            % Plot significant markers...
            
            for i = 1:length(wSCh)
                text(wSCh(i),(Y_Lim(1)+0.50),'o',...
                     'HorizontalAlignment','left','FontSize',20)
            end
        
        elseif FactorType == 3 %Special for interaction plots    
            
            xData = Data(:);
            
            for i = 1:numel(Data)
                    
                        y(i,:)   = squeeze(mean(cell2mat(cellfun(@(x) x(:,IdCh,:), xData(i,:), 'UniformOutput', false)')));
                        e(i,:)   = squeeze(std(cell2mat(cellfun(@(x) x(:,IdCh,:), xData(i,:), 'UniformOutput', false)')))/sqrt(handles.Info.Parameters.GroupSizes(1));
                    
            end
           
            time = (1:1:(handles.Info.Parameters.nSamples))*(1000/handles.Info.Parameters.rSample)-bL;
            % Assing the ERP data into Matlab for external graph making (E.g. excel)
            
            x = P_Values(IdCh,:);
            x(x>thresh) = 0;
            wSCh        = find(x)*(1000/handles.Info.Parameters.rSample)-bL; % Where Significant Channels (which samples are significant) 
            
            nL = size(Data,1); % Specify the number of lines in the plot
            x = time; % Rename so I can use the ERPErrorPlot nomenclature :D
            x = x(ones(numel(Data),1),:);
            
%             uL = y+e; %Upper Limit
%             lL = y-e; %Lower Limit
% 
%             yP = [lL, fliplr(uL)];
%             xP = [x,fliplr(x)];
            
            %% Actual ERP Plotting for the Interaction
            
            H.Figure = figure;
            set(H.Figure,...
                'Name',                         'ERP Comparison'    ,...
                'NumberTitle',                  'off'               ,...
                'Units',                        'pixels'            ,...
                'Position',                     [100 100 900 500]   ,...
                'Color',                        'w'                 ,...
                'Renderer',                     'openGL'            ,...
                'DefaultAxesColorOrder',        flipud(hsv(nL))     ,...
                'DefaultAxesLineStyleOrder',    {'-',':', '--', '-.'} );
            H.Axes = axes;
            set(H.Axes,...
                'FontName',         'Helvetica'         ,...
                'NextPlot',         'replacechildren'   );
            
%             hold on;

            H.Line = plot(H.Axes, x',y'); % Create the first set of lines
            set(H.Line,...
                'LineWidth',        3                   );
            
            Y_Lim = get(H.Axes, 'YLim');
            line([0,0],Y_Lim, 'Color', [0 0 0], 'LineStyle', '--')
            
%             plot(x',uL',':');
%             plot(x',lL',':');

%             H.Patch = patch(xP',yP',1);
%             set(H.Patch,...
%                 'EdgeColor',        'none'              ,...
%                 'FaceColor',        'flat'              ,...
%                 'FaceVertexCData',  flipud(hsv(nL))     ,...
%                 'CDataMapping',     'direct'            ,...
%                 'FaceAlpha',        0.15                );
%             
            hTitle  = title (['ERP Comparison for ' Ch]         );
            hXLabel = xlabel('Time (ms)'                        );
            hYLabel = ylabel('ERP Strength (MicroVolts)'        );
            
            % Create Legend Data
            for i = 1:size(Data,1)
                for j = 1:size(Data,2)
                    Legend{i,j} = ['Level' int2str(i) '-' int2str(j)];
                end
            end
            hLegend = legend(H.Line, Legend(:));
                      
%           set(H.Axes                       , ...
%                 'FontName'   , 'Helvetica' );
            set([hTitle, hXLabel, hYLabel] , ...
                'FontName'   , 'AvantGarde');
            
            set(hTitle                          , ...
                'FontSize'   , 14               , ...
                'FontWeight' , 'bold'           );
                                   
            set([hXLabel, hYLabel]              , ...
                'FontSize'   , 12               );
            
            set([hLegend, H.Axes]             , ...
                'FontSize'   , 10           );
            
            set(H.Axes, ...
                'Box'         , 'off'     , ...
                'TickDir'     , 'out'     , ...
                'TickLength'  , [.02 .02] , ...
                'XMinorTick'  , 'on'      , ...
                'YMinorTick'  , 'on'      , ...
                'YGrid'       , 'on'      , ...
                'XColor'      , [.3 .3 .3], ...
                'YColor'      , [.3 .3 .3], ...
                'LineWidth'   , .2        );
            
            % Plot significant markers...
            
            for i = 1:length(wSCh)
                text(wSCh(i),(Y_Lim(1)+0.50),'o',...
                     'HorizontalAlignment','left','FontSize',20)
            end
            
        end
        

case 'normal'

        tmpstr = get(gcbo, 'userdata');                 %creates temporary string holding userdata contents (label)
        set(gcbo, 'userdata', get(gcbo, 'string'));     %sets the the userdata to the objects string ('o')
        set(gcbo, 'string', tmpstr);                    %sets the objects string to the old userdata (label
        clear tmpstr;   %next time its called it does the reverse process :D

end

guidata(hObject,handles);

function d = dist_sph(vec,sensloc)
    R = vec(end);
    center = vec(1:end-1);
    % Average distance between the center if mass and the electrodes
    diffvert = bsxfun(@minus, sensloc, center);
    d = mean(abs(sqrt(sum(diffvert.^2,2)) - R));
    
function makeWaveletPlot(hObject, eventdata, handles)

seltype     = get(gcf,  'selectiontype');
Thresh      = str2double(get(handles.Txt_Thresh,'String'));

Ch          = get(gcbo, 'userData');
% e_labels    = {e_loc(:).labels};
IdCh        = find(strcmp(Ch,{handles.Info.Electrodes.e_loc(:).labels}));

switch seltype
    
    % if left-clicked display the channel label/point
    case 'normal'

        tmpstr = get(gcbo, 'userdata');                 %creates temporary string holding userdata contents (label)
        set(gcbo, 'userdata', get(gcbo, 'string'));     %sets the the userdata to the objects string ('o')
        set(gcbo, 'string', tmpstr);                    %sets the objects string to the old userdata (label
        clear tmpstr;   %next time its called it does the reverse process :D
    
    % if right-clicked display the wavelet summary plot
    case 'alt'    
        ept_WaveletChannelOverview(handles.Data, handles.Results, IdCh, Thresh)
        
end
        


function txt_Baseline_Callback(hObject, eventdata, handles)
% hObject    handle to txt_Baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function txt_Baseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_Baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in pop_Colormap.
function pop_Colormap_Callback(hObject, eventdata, handles)
% hObject    handle to pop_Colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Value');
    case 1
        colormap(hsv)
    case 2
        colormap(jet)
    case 3
        colormap(gray)
    case 4
        colormap(hot)
end     

% --- Executes during object creation, after setting all properties.
function pop_Colormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_Colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function menu_File_Callback(hObject, eventdata, handles)
