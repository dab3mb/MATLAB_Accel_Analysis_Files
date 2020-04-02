%% Set Up Training Data
% Import our models and combine into one data set
pinkyDataAccel = importfile("RecordedDataAccel-pinky30.csv",[1,inf]);
ringDataAccel = importfile("RecordedDataAccel-ring30.csv",[1,inf]);
middleDataAccel = importfile("RecordedDataAccel-middle30.csv",[1,inf]);
indexDataAccel = importfile("RecordedDataAccel-index30.csv",[1,inf]);
pinkyDataGyro = importfile("RecordedDataGyro-pinky30.csv",[1,inf]);
ringDataGyro = importfile("RecordedDataGyro-ring30.csv",[1,inf]);
middleDataGyro = importfile("RecordedDataGyro-middle30.csv",[1,inf]);
indexDataGyro = importfile("RecordedDataGyro-index30.csv",[1,inf]);

[pinkyXAccel, pinkyYAccel, pinkyZAccel, pinkyTimeAccel] = setUpAccel(pinkyDataAccel);
[ringXAccel, ringYAccel, ringZAccel, ringTimeAccel] = setUpAccel(ringDataAccel);
[middleXAccel, middleYAccel, middleZAccel, middleTimeAccel] = setUpAccel(middleDataAccel);
[indexXAccel, indexYAccel, indexZAccel, indexTimeAccel] = setUpAccel(indexDataAccel);

[pinkyXGyro, pinkyYGyro, pinkyZGyro, pinkyTimeGyro] = setUpGyro(pinkyDataGyro);
[ringXGyro, ringYGyro, ringZGyro, ringTimeGyro] = setUpGyro(ringDataGyro);
[middleXGyro, middleYGyro, middleZGyro, middleTimeGyro] = setUpGyro(middleDataGyro);
[indexXGyro,indexYGyro,indexZGyro, indexTimeGyro] = setUpGyro(indexDataGyro);


%% Plot Raw Data
figure(1)
subplot(2,2,1) , plot(pinkyTimeAccel, pinkyXAccel,'r')
hold on 
plot(pinkyTimeAccel, pinkyYAccel,'g')
plot(pinkyTimeAccel, pinkyZAccel,'b')
hold off
title("Pinky Accel Values and Time")
xlabel("Time in Seconds")
ylabel("Acceleration m/s^2")

subplot(2,2,2), plot(ringTimeAccel, ringXAccel,'r')
hold on 
plot(ringTimeAccel, ringYAccel,'g')
plot(ringTimeAccel, ringZAccel,'b')
hold off
title("Ring Accel Values and Time")
xlabel("Time in Seconds")
ylabel("Acceleration m/s^2")

subplot(2,2,3) , plot(middleTimeAccel, middleXAccel,'r')
hold on 
plot(middleTimeAccel, middleYAccel,'g')
plot(middleTimeAccel, middleZAccel,'b')
hold off
title("Middle Accel Values and Time")
xlabel("Time in Seconds")
ylabel("Acceleration m/s^2")

subplot(2,2,4) , plot(indexTimeAccel, indexXAccel,'r')
hold on 
plot(indexTimeAccel, indexYAccel,'g')
plot(indexTimeAccel, indexZAccel,'b')
hold off
title("Index Accel Values and Time")
xlabel("Time in Seconds")
ylabel("Acceleration m/s^2")


%% Apply Butterworth filter
% https://www.mathworks.com/help/signal/ref/butter.html#bucse3u-Wn

[pinkyFilteredDataXAccel, pinkyFilteredDataYAccel, pinkyFilteredDataZAccel] = applyButter(pinkyXAccel,pinkyYAccel, pinkyZAccel);
[ringFilteredDataXAccel, ringFilteredDataYAccel, ringFilteredDataZAccel] = applyButter(ringXAccel,ringYAccel, ringZAccel);
[middleFilteredDataXAccel, middleFilteredDataYAccel, middleFilteredDataZAccel] = applyButter(middleXAccel,middleYAccel, middleZAccel);
[indexFilteredDataXAccel, indexFilteredDataYAccel, indexFilteredDataZAccel] = applyButter(indexXAccel,indexYAccel, indexZAccel);

[pinkyFilteredDataXGyro, pinkyFilteredDataYGyro, pinkyFilteredDataZGyro] = applyButter(pinkyXGyro,pinkyYGyro, pinkyZGyro);
[ringFilteredDataXGyro, ringFilteredDataYGyro, ringFilteredDataZGyro] = applyButter(ringXGyro,ringYGyro, ringZGyro);
[middleFilteredDataXGyro, middleFilteredDataYGyro, middleFilteredDataZGyro] = applyButter(middleXGyro,middleYGyro, middleZGyro);
[indexFilteredDataXGyro, indexFilteredDataYGyro, indexFilteredDataZGyro] = applyButter(indexXGyro,indexYGyro, indexZGyro);


%% Plot Filtered Data
figure(2)
subplot(2,2,1) , plot(pinkyTimeAccel(3:end),pinkyFilteredDataXAccel, 'r')
hold on
plot(pinkyTimeAccel(3:end),pinkyFilteredDataYAccel, 'g')
plot(pinkyTimeAccel(3:end),pinkyFilteredDataZAccel, 'b')
hold off
title("Pinky with Butterworth Filter")
xlabel("Time in Seconds")
ylabel("Y-Direction Acceleration m/s^2")

subplot(2,2,2) , plot(ringTimeAccel(3:end),ringFilteredDataXAccel, 'r')
hold on
plot(ringTimeAccel(3:end),ringFilteredDataYAccel, 'g')
plot(ringTimeAccel(3:end),ringFilteredDataZAccel, 'b')
hold off
title("Ring with Butterworth Filter")
xlabel("Time in Seconds")
ylabel("Y-Direction Acceleration m/s^2")

subplot(2,2,3) , plot(middleTimeAccel(3:end),middleFilteredDataXAccel, 'r')
hold on
plot(middleTimeAccel(3:end),middleFilteredDataYAccel, 'g')
plot(middleTimeAccel(3:end),middleFilteredDataZAccel, 'b')
hold off
title("Middle with Butterworth Filter")
xlabel("Time in Seconds")
ylabel("Y-Direction Acceleration m/s^2")

subplot(2,2,4) , plot(indexTimeAccel(3:end),indexFilteredDataXAccel, 'r')
hold on
plot(indexTimeAccel(3:end),indexFilteredDataY, 'g')
plot(indexTimeAccel(3:end),indexFilteredDataZAccel, 'b')
hold off
title("Index with Butterworth Filter")
xlabel("Time in Seconds")
ylabel("Y-Direction Acceleration m/s^2")

%samplingFrequency = length(y)/max(time); % in Hz

%% Get Taps
% Find the average amplitude for each finger's z
pinkyAverageFilteredZAccel = sum(abs(pinkyFilteredDataZAccel))/length(pinkyFilteredDataZAccel);
ringAverageFilteredZAccel = sum(abs(ringFilteredDataZAccel))/length(ringFilteredDataZAccel);    
middleAverageFilteredZAccel = sum(abs(middleFilteredDataZAccel))/length(middleFilteredDataZAccel);    
indexAverageFilteredZAccel = sum(abs(indexFilteredDataZAccel))/length(indexFilteredDataZAccel);

% Find the average amplitude for each finger's x (gyro)
pinkyAverageFilteredXGyro = sum(abs(pinkyFilteredDataXGyro))/length(pinkyFilteredDataXGyro);
ringAverageFilteredXGyro = sum(abs(ringFilteredDataXGyro))/length(ringFilteredDataXGyro);    
middleAverageFilteredXGyro = sum(abs(middleFilteredDataXGyro))/length(middleFilteredDataXGyro);    
indexAverageFilteredXGyro = sum(abs(indexFilteredDataXGyro))/length(indexFilteredDataXGyro);

% Use this and the z data to find the 30 tap events' z data
[pinkyIndexOfTapsZAccel, pinkyValuesTapsZAccel] = getTaps(pinkyFilteredDataZAccel, pinkyAverageFilteredZAccel, 20, pinkyTimeAccel, 2.35);
[ringIndexOfTapsZAccel, ringValuesTapsZAccel] = getTaps(ringFilteredDataZAccel, ringAverageFilteredZAccel, 20, ringTimeAccel, 2.1);
[middleIndexOfTapsZAccel, middleValuesTapsZAccel] = getTaps(middleFilteredDataZAccel, middleAverageFilteredZAccel, 20, middleTimeAccel, 1.58);
[indexIndexOfTapsZAccel, indexValuesTapsZAccel] = getTaps(indexFilteredDataZAccel, indexAverageFilteredZAccel, 20, indexTimeAccel, 1.58);

% Use this and the z data to find the 30 tap events' x data (gyro)
[pinkyIndexOfTapsXGyro, pinkyValuesTapsXGyro] = getTaps(pinkyFilteredDataXGyro, pinkyAverageFilteredXGyro, 20, pinkyTimeGyro, 2.78);
[ringIndexOfTapsXGyro, ringValuesTapsXGyro] = getTaps(ringFilteredDataXGyro, ringAverageFilteredXGyro, 20, ringTimeGyro, 2.7);
[middleIndexOfTapsXGyro, middleValuesTapsXGyro] = getTaps(middleFilteredDataXGyro, middleAverageFilteredXGyro, 20, middleTimeGyro, 1.78);
[indexIndexOfTapsXGyro, indexValuesTapsXGyro] = getTaps(indexFilteredDataXGyro, indexAverageFilteredXGyro, 20, indexTimeGyro, 2.2);

% Use z data to get the x and y values
pinkyValuesTapsXAccel = pinkyFilteredDataXAccel(pinkyIndexOfTapsZAccel);
pinkyValuesTapsYAccel = pinkyFilteredDataYAccel(pinkyIndexOfTapsZAccel);
ringValuesTapsXAccel = ringFilteredDataXAccel(ringIndexOfTapsZAccel);
ringValuesTapsYAccel = ringFilteredDataYAccel(ringIndexOfTapsZAccel);
middleValuesTapsXAccel = middleFilteredDataXAccel(middleIndexOfTapsZAccel);
middleValuesTapsYAccel = middleFilteredDataYAccel(middleIndexOfTapsZAccel);
indexValuesTapsXAccel = indexFilteredDataXAccel(indexIndexOfTapsZAccel);
indexValuesTapsYAccel = indexFilteredDataYAccel(indexIndexOfTapsZAccel);

% Use x data to get the y and z values (gyro)
pinkyValuesTapsZGyro = pinkyFilteredDataZGyro(pinkyIndexOfTapsXGyro);
pinkyValuesTapsYGyro = pinkyFilteredDataYGyro(pinkyIndexOfTapsXGyro);
ringValuesTapsZGyro = ringFilteredDataZGyro(ringIndexOfTapsXGyro);
ringValuesTapsYGyro = ringFilteredDataYGyro(ringIndexOfTapsXGyro);
middleValuesTapsZGyro = middleFilteredDataZGyro(middleIndexOfTapsXGyro);
middleValuesTapsYGyro = middleFilteredDataYGyro(middleIndexOfTapsXGyro);
indexValuesTapsZGyro = indexFilteredDataZGyro(indexIndexOfTapsXGyro);
indexValuesTapsYGyro = indexFilteredDataYGyro(indexIndexOfTapsXGyro);

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
dataSet.Labels = [ (repelem(["PinkyXAccel"], 30)'); (repelem(["PinkyYAccel"], 30)'); (repelem(["PinkyZAccel"], 30)');
    (repelem(["RingXAccel"], 30)');(repelem(["RingYAccel"], 30)'); (repelem(["RingZAccel"], 30)'); (repelem(["MiddleXAccel"], 30)');  
    (repelem(["MiddleYAccel"], 30)'); (repelem(["MiddleZAccel"], 30)'); (repelem(["IndexXAccel"], 30)'); (repelem(["IndexYAccel"], 30)');  (repelem(["IndexZAccel"], 30)');(repelem(["PinkyXGyro"], 30)'); (repelem(["PinkyYGyro"], 30)'); (repelem(["PinkyZGyro"], 30)');
    (repelem(["RingXGyro"], 30)');(repelem(["RingYGyro"], 30)'); (repelem(["RingZGyro"], 30)'); (repelem(["MiddleXGyro"], 30)');  
    (repelem(["MiddleYGyro"], 30)'); (repelem(["MiddleZGyro"], 30)'); (repelem(["IndexXGyro"], 30)'); (repelem(["IndexYGyro"], 30)');  (repelem(["IndexZGyro"], 30)');];
dataSet.Data = [pinkyValuesTapsXAccel; pinkyValuesTapsYAccel; pinkyValuesTapsZAccel; ringValuesTapsXAccel; ringValuesTapsYAccel; ringValuesTapsZAccel; 
    middleValuesTapsXAccel; middleValuesTapsYAccel; middleValuesTapsZAccel; indexValuesTapsXAccel; indexValuesTapsYAccel; indexValuesTapsZAccel;pinkyValuesTapsXGyro; pinkyValuesTapsYGyro; pinkyValuesTapsZGyro; ringValuesTapsXGyro; ringValuesTapsYGyro; ringValuesTapsZGyro; 
    middleValuesTapsXGyro; middleValuesTapsYGyro; middleValuesTapsZGyro; indexValuesTapsXGyro; indexValuesTapsYGyro; indexValuesTapsZGyro;];

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
    class(labels)
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



function [xData, yData, zData, timeData] = setUpAccel(dataSet)
    xData = dataSet.x(1:end-15);
    yData = dataSet.y(1:end-15);
    zData = dataSet.z(1:end-15);
    zData = zData - 9;                  % Get rid of gravity
    timeData = dataSet.time(1:end-15); 
    timeData = timeData-timeData(1);        % Set starting time to 0
    timeData = timeData/1e+9;           % Change time to seconds
end
function [xData, yData, zData, timeData] = setUpGyro(dataSet)
    xData = dataSet.x(1:end-15);
    yData = dataSet.y(1:end-15);
    zData = dataSet.z(1:end-15);
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

