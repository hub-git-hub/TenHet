clc;
clear all;
addpath('tensor_toolbox');

load('DrugID.mat');
load('DrugChem.mat');
load('DrugSide.mat');

load('TargetID.mat');
load('TargetSW.mat');
load('TargetGO.mat');

load('DiseaseID.mat');
load('DiseaseMim.mat');
load('DiseaseHPO.mat');

fid = fopen('DrugTargetDisease');
rawdata = textscan(fid, '%q %q %q  %*[^\n]','headerlines',1);
fclose(fid);



%------------------------------------------------------------------------------
len = length(rawdata{1,1});
data = zeros(len,3);

% embed multiple side information to A
A{1}{1} = DrugChem;
A{1}{2} = DrugSide;

A{2}{1} = TargetSW;
A{2}{2} = TargetGO;

A{3}{1} = DiseaseMim;
A{3}{2} = DiseaseHPO;



% Mapping Drug/Target/Disease name to their ID
for i = 1:len
    data(i,1) = find(strcmp(DrugID(:,1),rawdata{1,1}{i}));
    data(i,2) = find(strcmp(TargetID(:,1),rawdata{1,2}{i}));
    data(i,3) = find(strcmp(DiseaseID(:,1),rawdata{1,3}{i}));
end

% split the data into 80% training set and 20% test set
folds = 5;
crossval_id = crossvalind('Kfold',(data(:,1)),5);
warning off;
for fold = 1:folds
	fprintf('---------------\nFOLD %d \n', fold)
    trainSample = find(crossval_id ~= fold);
    testSample = find(crossval_id == fold);
    trainData = data(trainSample, :);
    testData = data(testSample, :);
    obs = ones(length(testData),1) ;% observed value encode to 1 for test data

    % construct observed binary tensor X
    X = zeros(length(DrugID), length(TargetID), length(DiseaseID));
    for i = 1:length(trainData)
        X(trainData(i,1),trainData(i,2),trainData(i,3)) = 1;
    end
    
    Omega = logical(X);
    tic
    X(logical(1-Omega)) = mean(X(Omega));
    
   [X_hat, errX_hat] = TenHet(X, Omega, A, 0.001, 0.01, 0.01, 0.5, 50, 200);
   time = toc;
   
   % extract testdataset from completed tensor X_hat and compare with true  value
   pred = zeros(length(testData),1);
   for i = 1:length(testData)
       pred(i) = X_hat(testData(i,1),testData(i,2),testData(i,3));
   end
   
   %Running time and RMSE
   msrErr(fold) = sqrt(mean((pred(:)-obs(:)).^2));
   fprintf('The root mean square error is: %f.\n', msrErr(fold));
   times(fold) = time;
   
end

result.error = msrErr;
result.time = times;
       