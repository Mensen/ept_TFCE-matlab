% Load P_Values, TFCE_Obs, ChN
thresh      = 0.05;


x = TFCE_Obs; % copy TFCE_Obs
% x = CMass_Obs;
x(P_Values>thresh) = 0; % Threshold the data at alpha
Cp = ept_ClusRes(x, ChN, 0.01); % Calculate Positive Clusters

xn = x; % copy x
xn(x>0)=0; % Only show negative x
xn = abs(xn); % Make negative x's positive
Cn = ept_ClusRes(xn, ChN, 0.01); % Calculate Negative Clusters

C = Cn-Cp; % Combine the two
b = unique(C); % How many different clusters are there?
b(b==0)=[]; % Eliminate the 0 from being a unique cluster

nlogP = sum(-log(P_Values));

X = zeros(numel(b), length(nlogP));

for i = 1:numel(b)
    
    a           = P_Values;
    a(C~=b(i))  = 1;
    a           = sum(-log(a));
    
    X(i,:)      = a; 
    
end