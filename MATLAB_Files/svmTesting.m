%% Set Up Training Data
% Import our models and combine into one data set
pinkyTable = importModel("pinkyTable.csv");             %Poor naming of function, should be "importTable"
ringTable = importModel("ringTable.csv");
middleTable = importModel("middleTable.csv");
indexTable = importModel("indexTable.csv");
trainingTable = [pinkyTable; ringTable; middleTable; indexTable];

%% Create Model
% Creates SVM model! Have to use fitc ecoc bc it's more than 2 classes 
svmModel = fitcecoc(trainingTable, 'FingerString', 'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'));
% Check Cross Evaluation
crossEvalSVM = crossval(svmModel);
crossEvalLoss = kfoldLoss(crossEvalSVM);            % Roughly 33% Loss

%% Test Model
% Test on our model
pinkyRingMiddleIndexTest = importModel("TESTFILE1_FixedClasses.csv");
[label, score, cost] = predict(svmModel, pinkyRingMiddleIndexTest);
confusionchart(pinkyRingMiddleIndexTest.FingerString,label);            % 93% Wrong! 1/16 Correct :(((((
