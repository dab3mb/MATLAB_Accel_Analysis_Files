%% Set Up Training Data
% Import our models and combine into one data set
pinkyData = importfile("Pinkyx30.csv",[1,inf]);
ringData = importfile("Ringx30.csv",[1,inf]);
middleData = importfile("Middlex30.csv",[1,inf]);
indexData = importfile("Indexx30.csv",[1,inf]);
[pinkyX, pinkyY, pinkyZ, pinkyTime] = setUp(pinkyData);
[ringX, ringY, ringZ, ringTime] = setUp(ringData);
[middleX, middleY, middleZ, middleTime] = setUp(middleData);
[indexX, indexY, indexZ, indexTime] = setUp(indexData);

%% Plot Raw Data
figure(1)
subplot(2,2,1) , plot(pinkyTime, pinkyX,'r')
hold on 
plot(pinkyTime, pinkyY,'g')
plot(pinkyTime, pinkyZ,'b')
hold off
title("Pinky Accel Values and Time")
xlabel("Time in Seconds")
ylabel("Acceleration m/s^2")

subplot(2,2,2), plot(ringTime, ringX,'r')
hold on 
plot(ringTime, ringY,'g')
plot(ringTime, ringZ,'b')
hold off
title("Ring Accel Values and Time")
xlabel("Time in Seconds")
ylabel("Acceleration m/s^2")

subplot(2,2,3) , plot(middleTime, middleX,'r')
hold on 
plot(middleTime, middleY,'g')
plot(middleTime, middleZ,'b')
hold off
title("Middle Accel Values and Time")
xlabel("Time in Seconds")
ylabel("Acceleration m/s^2")

subplot(2,2,4) , plot(indexTime, indexX,'r')
hold on 
plot(indexTime, indexY,'g')
plot(indexTime, indexZ,'b')
hold off
title("Index Accel Values and Time")
xlabel("Time in Seconds")
ylabel("Acceleration m/s^2")


%% Apply and Plot Butterworth filter
% https://www.mathworks.com/help/signal/ref/butter.html#bucse3u-Wn

[pinkyFilteredDataX, pinkyFilteredDataY, pinkyFilteredDataZ] = applyButter(pinkyX,pinkyY, pinkyZ);
[ringFilteredDataX, ringFilteredDataY, ringFilteredDataZ] = applyButter(ringX,ringY, ringZ);
[middleFilteredDataX, middleFilteredDataY, middleFilteredDataZ] = applyButter(middleX,middleY, middleZ);
[indexFilteredDataX, indexFilteredDataY, indexFilteredDataZ] = applyButter(indexX,indexY, indexZ);

%% Plot Filtered Data
figure(2)
subplot(2,2,1) , plot(pinkyTime(3:end),pinkyFilteredDataX, 'r')
hold on
plot(pinkyTime(3:end),pinkyFilteredDataY, 'g')
plot(pinkyTime(3:end),pinkyFilteredDataZ, 'b')
hold off
title("Pinky with Butterworth Filter")
xlabel("Time in Seconds")
ylabel("Y-Direction Acceleration m/s^2")

subplot(2,2,2) , plot(ringTime(3:end),ringFilteredDataX, 'r')
hold on
plot(ringTime(3:end),ringFilteredDataY, 'g')
plot(ringTime(3:end),ringFilteredDataZ, 'b')
hold off
title("Ring with Butterworth Filter")
xlabel("Time in Seconds")
ylabel("Y-Direction Acceleration m/s^2")

subplot(2,2,3) , plot(middleTime(3:end),middleFilteredDataX, 'r')
hold on
plot(middleTime(3:end),middleFilteredDataY, 'g')
plot(middleTime(3:end),middleFilteredDataZ, 'b')
hold off
title("Middle with Butterworth Filter")
xlabel("Time in Seconds")
ylabel("Y-Direction Acceleration m/s^2")

subplot(2,2,4) , plot(indexTime(3:end),indexFilteredDataX, 'r')
hold on
plot(indexTime(3:end),indexFilteredDataY, 'g')
plot(indexTime(3:end),indexFilteredDataZ, 'b')
hold off
title("Index with Butterworth Filter")
xlabel("Time in Seconds")
ylabel("Y-Direction Acceleration m/s^2")

%samplingFrequency = length(y)/max(time); % in Hz

%% Get Taps
% Find the average amplitude for each finger's z
pinkyAverageFilteredZ = sum(abs(pinkyFilteredDataZ))/length(pinkyFilteredDataZ);
ringAverageFilteredZ = sum(abs(ringFilteredDataZ))/length(ringFilteredDataZ);    
middleAverageFilteredZ = sum(abs(middleFilteredDataZ))/length(middleFilteredDataZ);    
indexAverageFilteredZ = sum(abs(indexFilteredDataZ))/length(indexFilteredDataZ);
% Use this and the z data to find the 30 tap events' z data
[pinkyIndexOfTapsZ, pinkyValuesTapsZ] = getTaps(pinkyFilteredDataZ, pinkyAverageFilteredZ, 20, pinkyTime, .8);
[ringIndexOfTapsZ, ringValuesTapsZ] = getTaps(ringFilteredDataZ, ringAverageFilteredZ, 20, ringTime, .8);
[middleIndexOfTapsZ, middleValuesTapsZ] = getTaps(middleFilteredDataZ, middleAverageFilteredZ, 20, middleTime, .8);
[indexIndexOfTapsZ, indexValuesTapsZ] = getTaps(indexFilteredDataZ, indexAverageFilteredZ, 20, indexTime, .8);
% Use z data to get the x and y values
pinkyValuesTapsX = pinkyFilteredDataX(pinkyIndexOfTapsZ);
pinkyValuesTapsY = pinkyFilteredDataY(pinkyIndexOfTapsZ);
ringValuesTapsX = ringFilteredDataX(ringIndexOfTapsZ);
ringValuesTapsY = ringFilteredDataY(ringIndexOfTapsZ);
middleValuesTapsX = middleFilteredDataX(middleIndexOfTapsZ);
middleValuesTapsY = middleFilteredDataY(middleIndexOfTapsZ);
indexValuesTapsX = indexFilteredDataX(indexIndexOfTapsZ);
indexValuesTapsY = indexFilteredDataY(indexIndexOfTapsZ);
% Testing to see if it worked (it does lol)
% i = 1;
% zTimeForGraphing = [];
% yGraphing = [];
% while i <= size(indexIndexOfTapsZ, 1)
%     zTimeForGraphing = [zTimeForGraphing, indexIndexOfTapsZ(i, 1:end)];
%     yGraphing = [yGraphing, indexValuesTapsY(i, 1:end)];
%     i = i + 1;
% end
% plot(indexTime(zTimeForGraphing),yGraphing )
%% Combine Filtered Data
% This step is NOT accounting for chunk size variance. If the chunk size
% changes this block will fail :( 
dataSet = struct;
% smallestDataSize = min([size(pinkyFilteredDataX, 1), size(ringFilteredDataX, 1), size(middleFilteredDataX, 1), size(indexFilteredDataX, 1)]);
dataSet.Labels = [ (repelem(["PinkyX"], 30)'); (repelem(["PinkyY"], 30)'); (repelem(["PinkyZ"], 30)');
    (repelem(["RingX"], 30)');(repelem(["RingY"], 30)'); (repelem(["RingZ"], 30)'); (repelem(["MiddleX"], 30)');  
    (repelem(["MiddleY"], 30)'); (repelem(["MiddleZ"], 30)'); (repelem(["IndexX"], 30)'); (repelem(["IndexY"], 30)');  (repelem(["IndexZ"], 30)');];
dataSet.Data = [pinkyValuesTapsX; pinkyValuesTapsY; pinkyValuesTapsZ; ringValuesTapsX; ringValuesTapsY; ringValuesTapsZ; 
    middleValuesTapsX; middleValuesTapsY; middleValuesTapsZ; indexValuesTapsX; indexValuesTapsY; indexValuesTapsZ;];

%% Make images (using matlab function provided)
parentDir = fullfile(tempdir);
dataDir = 'cnnTestingImages';
helperCreateECGDirectories(dataSet,parentDir,dataDir)
helperCreateRGBfromTF(dataSet,parentDir,dataDir)
% Create an image datastore (image database)
allImages = imageDatastore(fullfile(parentDir,dataDir),...
    'IncludeSubfolders',true,...
    'LabelSource','foldernames');
%% Split image database into 80% Training, 20% Validation
rng default
[imgsTrain,imgsValidation] = splitEachLabel(allImages,0.8,'randomized');
disp(['Number of training images: ',num2str(numel(imgsTrain.Files))]);
disp(['Number of validation images: ',num2str(numel(imgsValidation.Files))]);
%% Loading google model
% Load in the pre-trained NN from GoogLeNet
net = googlenet;
lgraph = layerGraph(net);
numberOfLayers = numel(lgraph.Layers);
% figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
% plot(lgraph)
% title(['GoogLeNet Layer Graph: ',num2str(numberOfLayers),' Layers']);

% Create a dropout layer! Prevents overfitting of model
% Randomly sets inputs to zero 
newDropoutLayer = dropoutLayer(0.3,'Name','new_Dropout');
lgraph = replaceLayer(lgraph,'pool5-drop_7x7_s1',newDropoutLayer);

% Replace the layer 
numClasses = numel(categories(imgsTrain.Labels));
newConnectedLayer = fullyConnectedLayer(numClasses,'Name','new_fc',...
    'WeightLearnRateFactor',5,'BiasLearnRateFactor',5);
lgraph = replaceLayer(lgraph,'loss3-classifier',newConnectedLayer);

% Replace classification layer to one without labels
newClassLayer = classificationLayer('Name','new_classoutput');
lgraph = replaceLayer(lgraph,'output',newClassLayer);

% Set training options
options = trainingOptions('sgdm',...
    'MiniBatchSize',16,... % (16 - 58%)
    'MaxEpochs',35,... % (60 - 50%)
    'InitialLearnRate', .0001,... % (3 - 62.5%) (2 - 51.39%) 
    'ValidationData',imgsValidation,...
    'ValidationFrequency',10,...
    'Verbose',1,...
    'Plots','training-progress');
rng default


trainedGN = trainNetwork(imgsTrain,lgraph,options);
%% Evaluate Network
[YPred,probs] = classify(trainedGN,imgsValidation);
accuracy = mean(YPred==imgsValidation.Labels);
disp(['GoogLeNet Accuracy: ',num2str(100*accuracy),'%'])

% STOPPING HERE! 
% Something is wrong : (((((








%% Helper Functions (From MATLAB Guide + My Own)
function helperCreateRGBfromTF(ECGData,parentFolder,childFolder)
% This function is only intended to support the ECGAndDeepLearningExample.
% It may change or be removed in a future release.

imageRoot = fullfile(parentFolder,childFolder);

data = ECGData.Data;
labels = ECGData.Labels;

[~,signalLength] = size(data); %length of row (Number of columns)

fb = cwtfilterbank('SignalLength',signalLength,'VoicesPerOctave',12);
r = size(data,1); % Row number

for ii = 1:r % For every row
    cfs = abs(fb.wt(data(ii,:)));                                 % Get CWT Transform for row
    im = ind2rgb(im2uint8(rescale(cfs)),jet(128));                % Create imagefile
    
    imgLoc = fullfile(imageRoot,char(labels(ii)));                % Get name of folder where the image is going
    imFileName = strcat(char(labels(ii)),'_',num2str(ii),'.jpg'); % create file name
    imwrite(imresize(im,[224 224]),fullfile(imgLoc,imFileName));  % Write image file
end

end
function helperCreateECGDirectories(ECGData,parentFolder,dataFolder)
% This function is only intended to support the ECGAndDeepLearningExample.
% It may change or be removed in a future release.

rootFolder = parentFolder;
localFolder = dataFolder;
mkdir(fullfile(rootFolder,localFolder))

folderLabels = unique(ECGData.Labels);
for i = 1:numel(folderLabels)
    mkdir(fullfile(rootFolder,localFolder,char(folderLabels(i))));
end
end



function [xData, yData, zData, timeData] = setUp(dataSet)
    xData = dataSet.x(1:end-15);
    yData = dataSet.y(1:end-15);
    zData = dataSet.z(1:end-15);
    zData = zData - 9;                  % Get rid of gravity
    timeData = dataSet.time(1:end-15); 
    timeData = timeData-timeData(1);        % Set starting time to 0
    timeData = timeData/1e+9;           % Change time to seconds
end
function [filteredDataX, filteredDataY, filteredDataZ] = applyButter(datax,datay, dataz)
    % Butterworth Filter Order (Want fast cutoff for cleaner data)
    n = 16;      
    % Cutoff Frequency (0:1 where 1 is Half sampling rate)
    % Can be thought of what percent of our frequencies we're letting through
    Wn = .35;
    % Apply Butterworth filter
    [b,a] = butter(n,Wn,"low");
    filteredDataX = filter(b,a,datax);
    filteredDataX = filteredDataX(3:end);
    filteredDataY = filter(b,a,datay);
    filteredDataY = filteredDataY(3:end);
    filteredDataZ = filter(b,a,dataz);
    filteredDataZ = filteredDataZ(3:end);
end
function [indexList, valuesList] = getTaps(direction, averageDirection, chunkSize, time, significance)
    i = 1;
    indexList = [];
    valuesList = [];
    while i < length(direction)
        % Set Lower and Upper bound of our selection
        begin = i;
        ending = 0;        
        % Check if significantly different
        if abs(direction(i)) > averageDirection + (averageDirection*significance) % Value must be (significance)% more than average  
            % Try to get a little before 
            if (i - 10) >= 1                 
                begin = (i - 10);
            end
            % Check if prefered ending is in bounds
            ending = i + chunkSize;
            if (i + ending > length(direction))
                % If it doesn't, find an ending is in bounds
                while (ending > length(direction))
                    ending = ending - 1;
                end
            end
            % Add to our lists
            indexList = [indexList; begin:ending];
            valuesList = [valuesList, direction(begin:ending)];
%             for j=begin:ending
%                 valuesList = [valuesList, direction(j)];
%             end
        end
        % If we saved a chunk, make sure we don't read the chunk again
        if (ending > 0)
            i = ending;
        else
            i = i + 1;
        end
    end
    valuesList = valuesList';
end

%% Trash
% CNN is mostly used for image classification so we will create images from our FFT models
% https://www.mathworks.com/help/wavelet/examples/classify-time-series-using-wavelet-analysis-and-deep-learning.html
% durationOfSignal = seconds(pinkyTime(3:end));
% signal = timetable(durationOfSignal, pinkyFilteredDataY);
% defineRetimeStep = seconds(0.0097181);
% signalRetimed = retime(signal, 'regular', 'linear', 'TimeStep', defineRetimeStep);
% cwtModel = cwt(signalRetimed);

