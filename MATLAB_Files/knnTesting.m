%% Set Up Training Data
% Import our models and combine into one data set

load('accelDataindex30.mat');
load('accelDatamiddle30.mat')
load('accelDataring30.mat')
load('accelDatapinky30.mat')
load('gyroDataindex30.mat');
load('gyroDatamiddle30.mat')
load('gyroDataring30.mat')
load('gyroDatapinky30.mat')


indexTable = [removevars(Final_Model_Accel_INDEX30, {'Finger'}),Final_Model_GYRO_INDEX30];
middleTable = [removevars(Final_Model_Accel_MIDDLE30, {'Finger'}),Final_Model_GYRO_MIDDLE30];
ringTable = [removevars(Final_Model_Accel_RING30, {'Finger'}),Final_Model_GYRO_RING30];
pinkyTable = [removevars(Final_Model_Accel_PINKY30, {'Finger'}),Final_Model_GYRO_PINKY30];

trainingTable = [indexTable; middleTable; ringTable; pinkyTable];


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
knnModel = fitcknn(trainingTable, 'Finger','OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));
%% Resubstiution and Cross Evaluation Loss
% Test model for resub loss and cross eval loss
resubLossAmount = resubLoss(knnModel);                         % Roughly 1.6%
crossValidatedClassifier = crossval(knnModel);                  
crossValLossValue = kfoldLoss(crossValidatedClassifier);       % Roughly 3.3% Loss <-- Is this okay?

%% Testing Model
pinkyRingMiddleIndexTest = [removevars(Final_Model_Accel_TESTING16, {'Finger'}),Final_Model_GYRO_TESTING16];
[label, score, cost] = predict(knnModel, pinkyRingMiddleIndexTest);
confusionchart({'p','p','p','p','r','r','r','r','m','m','m','m','i','i','i','i'}, label);            % 81% Correct!
%     {'p'} p 
%     {'p'} p 
%     {'p'} p 
%     {'p'} p 
%     {'p'} r x
%     {'r'} r 
%     {'r'} r 
%     {'p'} r x
%     {'i'} m x
%     {'m'} m 
%     {'m'} m 
%     {'m'} m 
%     {'i'} i 
%     {'i'} i 
%     {'i'} i 
%     {'i'} i  




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



