clear all; clc; close all
FsFP = 300; %Hz

Path = '.\Dropfoot\';
File = 'DCM\Descalco25_170831_1';
FilePath = [Path,File];
csv = '-Delsys 1.csv';
numSensors = [11,12];

% Cortex data
%     Forces = importdata([FilePath,'.forces']);
%     Fy(:,1) = (Forces.data(:,1))/FsFP;
%     Fy(:,2) = Forces.data(:,10)/2;
    
    Contact1 = ReadDelsys1([FilePath,csv], 'EMG', 'EMG',numSensors);
    Contact2 = ReadDelsys1([FilePath,csv], 'ACC', 'ACC Pitch',numSensors);
    Contact3 = ReadDelsys1([FilePath,csv], 'ACC', 'ACC Roll',numSensors);
    Contact4 = ReadDelsys1([FilePath,csv], 'ACC', 'ACC Yaw',numSensors);
   
    
    plot(Fy(:,1),Fy(:,2))