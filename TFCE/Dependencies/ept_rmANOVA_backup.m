function [F, Data] = ept_rmANOVA(Data, rand)
%% ANOVA: two-way repeated measures (both within) with reshaping matrices

% This file is part of the program ept_TFCE.
% ept_TFCE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% ept_TFCE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with ept_TFCE.  If not, see <http://www.gnu.org/licenses/>.

if nargin == 0
    display ('Showing test results for random dataset since no input was given');
    Data{1,1} = randi(20,14,125,249);
    Data{1,2} = randi(20,14,125,249)+3;
    Data{1,3} = randi(20,14,125,249);
    Data{2,1} = randi(20,14,125,249);
    Data{2,2} = randi(20,14,125,249);
    Data{2,3} = randi(20,14,125,249)+1;
    
%     Data{1,1} = randi(10,10,6,2);
%     Data{1,2} = randi(10,10,6,2)+3;
%     Data{1,3} = randi(10,10,6,2);
%     Data{2,1} = randi(10,10,6,2);
%     Data{2,2} = randi(10,10,6,2);
%     Data{2,3} = randi(10,10,6,2)+1;

end

if nargin == 1;
    rand = 0;
end

nF1     = size(Data,1);
nF2     = size(Data,2);
nP      = size(Data{1,1},1); % number of participants
nC      = size(Data{1,1},2); % number of channels
nT      = size(Data{1,1},3); % number of time points (samples)

if rand == 0;
    nData   = cell2mat(Data);    % change cells to a matrix for fast calculation
    nData   = reshape (nData,[nP, nF1, nC, nF2, nT]);
    nData   = permute (nData, [2,4,1,3,5]); % re-order the data for appropriate sum calculations
end

if rand == 1;
    nData   = cell2mat(Data(:)); % Vectorise
    nData   = reshape (nData,[nP,nF1, nF2, nC, nT]);

    % generate a new random relabelling for each participant
    randLevels = arrayfun(@(x)randperm(nF1*nF2),1:nP,'UniformOutput',false)';

    randData = zeros(size(nData)); % allocate memory for speed
    for i = 1:nP
        randData(i,:,:) = nData(i,randLevels{i,1},:);
    end
    
    nData = reshape(randData, [nP, nF1, nF2, nC, nT]);
    nData = permute (nData, [2, 3, 1, 4, 5]);
    % reshape to fit to original scheme
    
end
    
% Degrees of Freedom
dfA = nF1-1;
dfB = nF2-1;
dfS = nP-1;
dfAB = dfA*dfB;
dfAS = dfA*dfS;
dfBS = dfB*dfS;
dfABS = dfAB*dfS;


% Row and Column Sums
sBS = squeeze(sum(nData,1));
sAS = squeeze(sum(nData,2));
sAB = squeeze(sum(nData,3));

sS = squeeze(sum(sAS,1));
sA = squeeze(sum(sAB,2));
sB = squeeze(sum(sAB,1));
sT = squeeze(sum(sA,1));

% expected Values
if nT == 1; %if there is only one time point then squeeze shapes the matrix backward
    eA = squeeze(sum(sA.^2)./(nF2*nP));
    eB = squeeze(sum(sB.^2)./(nF1*nP));
    eAB = squeeze(sum(sum(sAB.^2))./nP)';
    eS = squeeze(sum(sS.^2)./(nF1*nF2));
    eAS = squeeze(sum(sum(sAS.^2))./nF2)';
    eBS = squeeze(sum(sum(sBS.^2))./nF1)';
    eY = squeeze(sum(sum(sum(nData.^2))))'; % sum(Y.^2); %that works too... so why the other more complex calculation?
    eT = sT.^2 / (nF1*nF2*nP);  
else % all other cases where there are more than one time point
    eA = squeeze(sum(sA.^2)./(nF2*nP));
    eB = squeeze(sum(sB.^2)./(nF1*nP));
    eAB = squeeze(sum(sum(sAB.^2))./nP);
    eS = squeeze(sum(sS.^2)./(nF1*nF2));
    eAS = squeeze(sum(sum(sAS.^2))./nF2);
    eBS = squeeze(sum(sum(sBS.^2))./nF1);
    eY = squeeze(sum(sum(sum(nData.^2)))); % sum(Y.^2); %that works too... so why the other more complex calculation?
    eT = sT.^2 / (nF1*nF2*nP);
end

% Sum of Squares
ssA = eA - eT;
ssB = eB - eT;
ssAB = eAB - eA - eB + eT;
% ssS = eS - eT;
ssAS = eAS - eA - eS + eT;
ssBS = eBS - eB - eS + eT;
ssABS = eY - eAB - eAS - eBS + eA + eB + eS - eT;
% ssTot = eY - eT;

% Mean Square
msA = ssA / dfA;
msB = ssB / dfB;
msAB = ssAB / dfAB;
% msS = ssS / dfS;
msAS = ssAS / dfAS;
msBS = ssBS / dfBS;
msABS = ssABS / dfABS;

% f statistic
F.A = msA./ msAS;
F.B = msB./ msBS;
F.AB = msAB./ msABS;

% Get values for specific channel and time
% gC = 9;
% gT = 151;
% 
% for j = 1:2;
%     for i = 1:3;
%         
%         nData{j,i} = Data{j,i}(:,gC,gT);
% 
%     end
% end

% Replicate Data
% gC = 2;
% gT = 1;
% 
% for j = 1:2;
%     for i = 1:2;
%         
%         newData{j,i} = repmat(Data{j,i},[1,gC,gT]);
% 
%     end
% end

end
