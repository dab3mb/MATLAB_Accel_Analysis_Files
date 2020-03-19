%% Set Up Training Data
% Import our models and combine into one data set
pinkyTable = importModel("pinkyTable.csv");             %Poor naming of function, should be "importTable"
ringTable = importModel("ringTable.csv");
middleTable = importModel("middleTable.csv");
indexTable = importModel("indexTable.csv");
trainingTable = [pinkyTable; ringTable; middleTable; indexTable];

% Rn the data looks really bad? Like I don't see any differences in it :/
% scatter(1:30, pinkyModel.StrongestFreqHzZ,'b','*')
% hold on
% scatter(1:30, ringModel.StrongestFreqHzZ,'g','+')
% scatter(1:30, middleModel.StrongestFreqHzZ,'r', 'x')
% scatter(1:30, indexModel.StrongestFreqHzZ,'p', '.')
% hold off

%% MATLAB Create model and find best k and distance 
% https://www.mathworks.com/help/stats/fitcknn.html#namevaluepairarguments
rng(1)
knnModel = fitcknn(trainingTable, 'FingerString','OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));
%% Resubstiution and Cross Evaluation Loss
% Test model for resub loss and cross eval loss
resubLossAmount = resubLoss(knnModel);                         % Roughly 14% Loss <-- Is this okay? 
crossValidatedClassifier = crossval(knnModel);                  
crossValLossValue = kfoldLoss(crossValidatedClassifier);       % Roughly 34% Loss <-- Is this okay?

%% Testing Model
pinkyRingMiddleIndexTest = importModel("TESTFILE1_FixedClasses.csv");
[label, score, cost] = predict(knnModel, pinkyRingMiddleIndexTest);
confusionchart(pinkyRingMiddleIndexTest.FingerString,label);            % 4/16 correct??? 75 % Wrong





%% Find best k
% k=1;
% kValues = [];
% resubLossValues = []; % Fraction of misclassifications from the predictions
% crossValLossValues = []; % Fraction of misclassifications from cross val
% while k < height(trainingTable)
%     trainingModel.NumNeighbors = k;
%     kValues = [kValues ; k];
%     resubLossValues = [resubLossValues; resubLoss(trainingModel)];
%     crossValidatedClassifier = crossval(trainingModel);
%     crossValLossValues = [crossValLossValues; kfoldLoss(crossValidatedClassifier)];
%     k = k + 1;
%     k
% end
% kValuesTable = table(kValues, resubLossValues, crossValLossValues);
% ------ Scrapped bc MATLAB has a built in function to do this LOL ------



