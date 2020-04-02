%% Set Up Data               % -15 to cut off the touching of stop button
fingerNumber = "index30";
whichFinger  = "i";
data = importfile("RecordedDataGyro-" + fingerNumber +  ".csv",[1,inf]);
x = data.x(15:end-15);
y = data.y(15:end-15);
z = data.z(15:end-15);
time = data.time(15:end-15); 
time = time-time(1);        % Set starting time to 0
time = time/1e+9;           % Change time to seconds

% Plot Original Time
subplot(2,2,1) , plot(time, x,'r')
hold on 
plot(time, y,'g')
plot(time, z,'b')
hold off

title("Rotation Values and Time")
xlabel("Time in Seconds")
ylabel("Rate of rotation around the axes in rad/s")



%% Apply and Plot Butterworth filter
% https://www.mathworks.com/help/signal/ref/butter.html#bucse3u-Wn
% Butterworth Filter Order (Want fast cutoff for cleaner data)
n = 16;      
% Cutoff Frequency (0:1 where 1 is Half sampling rate)
% Can be thought of what percent of our frequencies we're letting through
Wn = .35;

% Apply Butterworth filter
[b,a] = butter(n,Wn,"low");
filteredDataX = filter(b,a,x);
filteredDataX = filteredDataX(3:end);
filteredDataY = filter(b,a,y);
filteredDataY = filteredDataY(3:end);
filteredDataZ = filter(b,a,z);
filteredDataZ = filteredDataZ(3:end);
% Plot data
subplot(2,2,2) , plot(time(3:end),filteredDataX, 'r')
hold on
plot(time(3:end),filteredDataY, 'g')
plot(time(3:end),filteredDataZ, 'b')
hold off

title("Butterworth Filter")
xlabel("Time in Seconds")
ylabel("Rate of rot. rad/s")

samplingFrequency = length(y)/max(time); % in Hz


%% Plotting our FFT WITHOUT Butterworth filter
% https://www.youtube.com/watch?v=dM1y6ZfQkDU&t=389s
fouriedY = fft(y);  
fouriedX = fft(x);
fouriedZ = fft(z);
% Returns a complex number with Phase and Amplitude. We want the amplitude,
% and using absolute value can help do so 
LY = length(y);
twoSidedSpecY = abs(fouriedY/LY);
oneSidedSpecY = twoSidedSpecY(1:LY/2+1);
oneSidedSpecY(2:end-1) = 2*oneSidedSpecY(2:end-1);
frequencyDomainY = samplingFrequency * (0:(LY/2))/LY;

LX = length(x);
twoSidedSpecX = abs(fouriedX/LX);
oneSidedSpecX = twoSidedSpecX(1:LX/2+1);
oneSidedSpecX(2:end-1) = 2*oneSidedSpecX(2:end-1);
frequencyDomainX= samplingFrequency * (0:(LX/2))/LX;

LZ = length(z);
twoSidedSpecZ = abs(fouriedZ/LZ);
oneSidedSpecZ = twoSidedSpecZ(1:LZ/2+1);
oneSidedSpecZ(2:end-1) = 2*oneSidedSpecZ(2:end-1);
frequencyDomainZ = samplingFrequency * (0:(LZ/2))/LZ;

% Get rid of 0-1 in this one, it's making our data hard to read
% We also only need half of it bc the other half is a mirror image!
subplot(2,2,3) , plot(frequencyDomainX(3:end),oneSidedSpecX(3:end),'r')
hold on
plot(frequencyDomainY(3:end),oneSidedSpecY(3:end),'g')
plot(frequencyDomainZ(3:end),oneSidedSpecZ(3:end),'b')
hold off

title("Fourier Transform")
xlabel("Frequency in Hz")
ylabel("Magnitude (decibels)")



%% Fourier WITH Butterworth filter
% https://www.youtube.com/watch?v=dM1y6ZfQkDU&t=389s
fouriedY = fft(filteredDataY);  
fouriedX = fft(filteredDataX);
fouriedZ = fft(filteredDataZ);
% Returns a complex number with Phase and Amplitude. We want the amplitude,
% and using absolute value can help do so 

% % LX = length(x);
% % twoSidedSpecX = abs(fouriedX/LX);
% % oneSidedSpecX = twoSidedSpecX(1:LX/2+1);
% % oneSidedSpecX(2:end-1) = 2*oneSidedSpecX(2:end-1);
% % frequencyDomainX= samplingFrequency * (0:(LX/2))/LX;

LY = length(y);
twoSidedSpecY = abs(fouriedY/LY);
oneSidedSpecY = twoSidedSpecY(1:LY/2+1);
oneSidedSpecY(2:end-1) = 2*oneSidedSpecY(2:end-1);
frequencyDomainY = samplingFrequency * (0:(LY/2))/LY;

LZ = length(z);
twoSidedSpecZ = abs(fouriedZ/LZ);
oneSidedSpecZ = twoSidedSpecZ(1:LZ/2+1);
oneSidedSpecZ(2:end-1) = 2*oneSidedSpecZ(2:end-1);
frequencyDomainZ = samplingFrequency * (0:(LZ/2))/LZ;

% Get rid of 0-1 in this one, it's making our data hard to read
% We also only need half of it bc the other half is a mirror image!
subplot(2,2,4) , plot(frequencyDomainX(3:end/2),oneSidedSpecX(3:end/2),'r')
hold on
plot(frequencyDomainY(3:end/2),oneSidedSpecY(3:end/2),'g')
plot(frequencyDomainZ(3:end/2),oneSidedSpecZ(3:end/2),'b')
hold off

title("Fourier Transform After ButterWorth")
xlabel("Frequency in Hz")
ylabel("Magnitude (decibels)")




%% Get Taps

averageFilteredX = sum(abs(filteredDataX))/length(filteredDataX);
% averageFilteredY = sum(abs(filteredDataY))/length(filteredDataY);    
% averageFilteredZ = sum(abs(filteredDataZ))/length(filteredDataZ);    

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% If I have time later, make it find the significance based on the number
% of samples you say you have
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Use z data to find 'events' or 'taps'
[xTimeIndex,xData,sampleX]= getTaps(filteredDataX, averageFilteredX, 28, time, 1.8, "_gyro_x", whichFinger);                

% Use 'event' or 'tap' times from z (zTime) to get the data for x & y
sampleZ = getZYFromX(filteredDataZ, xTimeIndex, time, "_gyro_z");
sampleY = getZYFromX(filteredDataY, xTimeIndex, time, "_gyro_y");


% Create final Model
Final_Model_GYRO_INDEX30 = [sampleZ(:,1:end), sampleY(:,1:end), sampleX(:, 1:end)];
Final_Model_GYRO_INDEX30
% Get Variables for graphing
i = 1;
xTimeIndexForGraphing = [];
while i <= size(xTimeIndex, 1)
    xTimeIndexForGraphing = [xTimeIndexForGraphing, xTimeIndex(i, 1:end)];
    i = i + 1;
end

%Graph Model
figure(2)
hold on
plot(time(xTimeIndexForGraphing), filteredDataZ(xTimeIndexForGraphing),'r')
plot(time(xTimeIndexForGraphing), filteredDataY(xTimeIndexForGraphing),'g')
plot(time(xTimeIndexForGraphing), xData,'b')
hold off
title("All Directions, Taps Only")
xlabel("Time in Seconds")
ylabel("Rotation Speed ")





% ---- Write our final file -----
fileName = 'gyroData'+ fingerNumber;
save(fileName , 'Final_Model_GYRO_INDEX30')

function [indexList, valuesList, samples] = getTaps(direction, averageDirection, chunkSize, time, significance, fingerAndDirection, justFinger)
    i = 1;
    indexList = [];
    valuesList = [];
    FreqRangeHz = [];
    freqAcrossEachChunk = [];
    LengthSec = [];
    FingerString = [];
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
            % Add to our index list
            indexList = [indexList; begin:ending];
            for j=begin:ending
                valuesList = [valuesList, direction(j)];
            end
            
            % ----- Get frequency data from chunk -----
            fouriedChunk = fft(direction(begin:ending));
            timeChunk = time(begin:ending);
            timeChunk = timeChunk-timeChunk(1);
            lengthChunk = length(fouriedChunk);
            samplingFrequencyChunk = lengthChunk/max(timeChunk); % in Hz
            twoSidedSpecChunk = abs(fouriedChunk/lengthChunk);
            oneSidedSpecChunk = twoSidedSpecChunk(1:lengthChunk/2+1);
            oneSidedSpecChunk(2:end-1) = 2*oneSidedSpecChunk(2:end-1);
            frequencyDomainChunk = samplingFrequencyChunk * (0:(lengthChunk/2))/lengthChunk;
            oneSidedSpecChunk = oneSidedSpecChunk(3:end); % Trim off excess data
            frequencyDomainChunk = frequencyDomainChunk(3:end);% Trim off excess data
            
            
            % ----- Put Chunk Freq Data into a , put that in -----
            freqAcrossEachChunk = [freqAcrossEachChunk; oneSidedSpecChunk'];
            
           
            
             ii = 0;
             while (ii < size(oneSidedSpecChunk,1 ))
                FingerString = [FingerString; fingerAndDirection];
                ii = ii + 1;
             end
        end
        % If we saved a chunk, make sure we don't read the chunk again
        if (ending > 0)
            i = ending;
        else
            i = i + 1;
        end
    end
    % ----- Format Table -----
    tableLabels = [];
    for i = 1:size(freqAcrossEachChunk, 2)
        tableLabels = [tableLabels, fingerAndDirection+int2str(i)];
    end
    tableLabels = [tableLabels, "Finger"];
    samples = array2table(freqAcrossEachChunk);
    samples = addvars(samples, repelem(justFinger, size(freqAcrossEachChunk, 1))');
    samples.Properties.VariableNames = tableLabels(1:end);
    samples
    
end

function samples = getZYFromX(direction, zTimeDataIndex, time, fingerAndDirection)
    freqAcrossEachChunk = [];
    startChunkIndex = 1;
    endChunkIndex = size(zTimeDataIndex, 1);
    row = [];
    % Loop through each row
    rowNumber = 1;
    while rowNumber <= size(zTimeDataIndex, 1)
        row = [zTimeDataIndex(rowNumber, 1:end)];
%         
        % ----- Get frequency data from chunk -----
        fouriedChunk = fft(direction(row));
        timeChunk = time(row);
        timeChunk = timeChunk-timeChunk(1);
        lengthChunk = length(fouriedChunk);
        samplingFrequencyChunk = lengthChunk/max(timeChunk); % in Hz
        twoSidedSpecChunk = abs(fouriedChunk/lengthChunk);
        oneSidedSpecChunk = twoSidedSpecChunk(1:lengthChunk/2+1);
        oneSidedSpecChunk(2:end-1) = 2*oneSidedSpecChunk(2:end-1);
        frequencyDomainChunk = samplingFrequencyChunk * (0:(lengthChunk/2))/lengthChunk;
        oneSidedSpecChunk = oneSidedSpecChunk(3:end);       % Trim off excess data
        frequencyDomainChunk = frequencyDomainChunk(3:end); % Trim off excess data

        %domainAndFrequency = [frequencyDomainChunk',oneSidedSpecChunk];
        freqAcrossEachChunk = [freqAcrossEachChunk; oneSidedSpecChunk'];


        rowNumber = rowNumber + 1;
    end
    % ----- Format Table -----
    tableLabels = [];
    for i = 1:size(freqAcrossEachChunk, 2)
        tableLabels = [tableLabels, fingerAndDirection+int2str(i)];
    end
    samples = array2table(freqAcrossEachChunk);
    samples.Properties.VariableNames = tableLabels(1:end);
    samples
end