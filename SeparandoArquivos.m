clear all; clc; close all
FsFP = 300; %Hz

Path = '.\Dropfoot\';
File = 'DCM\Descalco25_170831_1';
FilePath = [Path,File];
csv = '-Delsys 1.csv';
xls = '-Delsys 1.xlsx';
name = 'Descalco25_170831_';
% numSensors = [11,12];

A =  importdata([FilePath,csv]);

Data = [A.textdata, num2cell(A.data)];

Label = Data(1,:);

len = 31512; 
num = 0;
for i = 2: len-1 : length(Data)-len
    num = num +1;
    newData = zeros(len+1,20);
    newData = [Label;Data(i:i+len-1,:)];
    xlswrite([name,num,xls],newData)
end
