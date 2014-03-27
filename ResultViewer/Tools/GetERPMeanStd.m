%% Extract ERP Mean and Std from Dataset

% load e_loc file
Ch          = 'E2';
e_labels    = {e_loc(:).labels};
IdCh        = strcmp(Ch,e_labels);

% For multi-factorial studies...
    % Select Factor Type A=1; B=2
    Type = 1; 

    switch Type
        case 1
            AvgA = squeeze(mean(cell2mat(Data(1,:)')));
            AvgB = squeeze(mean(cell2mat(Data(2,:)')));
            StdA = squeeze(std(cell2mat(Data(1,:)')));
            StdB = squeeze(std(cell2mat(Data(2,:)')));
        case 2
            AvgA = squeeze(mean(cell2mat(Data(:,1))));
            AvgB = squeeze(mean(cell2mat(Data(:,2))));
            StdA = squeeze(std(cell2mat(Data(:,1))));
            StdB = squeeze(std(cell2mat(Data(:,2))));
        case 3

    end 

    Export(1,:)  = AvgA(IdCh,:);
    Export(2,:)  = AvgB(IdCh,:);
    Export(3,:)  = StdA(IdCh,:);
    Export(4,:)  = StdB(IdCh,:);

%% Export Log Significance
% Load P_Values

Exportp(1,:) = sum(-log(P_Values.A));
Exportp(2,:) = sum(-log(P_Values.B));
Exportp(3,:) = sum(-log(P_Values.AB));

%% Extract Individual ERPs for comparison

% load e_loc file
Ch          = 'E129';
e_labels    = {e_loc(:).labels};
IdCh        = strcmp(Ch,e_labels);

% For multi-factorial studies...
    % Select Factor Type A=1; B=2
    Type = 1; 
    
    switch Type
        case 1
            A = Data(1,:)';
            AvgA = (A{1}+A{2})./2;
            B = Data(2,:)';
            AvgB = (B{1}+B{2})./2;
        case 2
            A = Data(:,1);
            AvgA = (A{1}+A{2})./2;
            B = Data(:,2);
            AvgB = (B{1}+B{2})./2;
        case 3
            % Interaction
    end
    
% Get channel from average

L1 = squeeze(AvgA(:,IdCh,:));
L2 = squeeze(AvgB(:,IdCh,:));

