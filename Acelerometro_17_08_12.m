%% Gait Analysis (w/ accelerometer)
clear all; clc; close all

%% -- About
% Date: 3Aug2017
% Subject: EGS
% KISTLER Sensor: - 2150551 (Acc1) - Shank L
%                 - 2148496 (Acc2) - Shank R
% DELSYS Sensor:  - Sensor 1  -  Shank L
%                 - Sensor 2  -  Shank R
% --- First support on FP : always with the RIGHT foot

%% -- Load data
% Delsys data
% Name = 'Acc_170803_EGS_';
% Trial = {'1' '2' '3' '4' '5'}
% ShankL = 2;ShankR = 3;
% Name = 'Acc_170731_RNW_';
% csv = '-Delsys 1.csv';
% Trial = {'1' '2' '4' '5'};
Name = 'Acc_170731_DCM_';
csv = '-Delsys 1.csv';
Trial = {'1' '2' '3' '4' };
ShankL = 5;ShankR = 4;

Path = 'C:\Users\desir\Desktop\Acelerometro - Aug17\Coletas\';
ChannelType = 'AUX';
Signal = {'IM ACC Pitch', 'IM ACC Roll', 'IM ACC Yaw'};

p = zeros(149,5,2*length(Trial));
y = zeros(149,2*length(Trial));

for k = 1: length(Trial)
    File = [Name,Trial{k}];
    FilePath = [Path,File];
    
    Pitch = ReadDelsys([FilePath,csv], ChannelType, Signal(1));
    Roll = ReadDelsys([FilePath,csv], ChannelType, Signal(2));
    Yaw = ReadDelsys([FilePath,csv], ChannelType, Signal(3));
    EMG = ReadDelsys([FilePath,csv], 'IMEMG', 'EMG');
    
    % Cortex data
    Forces = importdata([FilePath,'.forces']);
    Fy = (Forces.data(:,1) -1)/300;
    for i =1 : length(Forces.colheaders)
        if strcmp(Forces.colheaders{i}(1:2), 'FY')
            Fy = [Fy, Forces.data(:,i)];
        end
    end
    
    % ANC file
    % ANC = importdata([FileName,'.anc'],'\t',11);
    %
    % Accx1 = ANC.data(:,29);
    % Accy1 = ANC.data(:,30);
    % Accz1 = ANC.data(:,31);
    % Accx2 = ANC.data(:,32);
    % Accy2 = ANC.data(:,33);
    % Accz2 = ANC.data(:,34);
    %
    % Acc1 = ANC.data(:,1);
    
    %% -- Filtering
    
    Fs = 148.39;
    [n,Wn] = buttord(10/(Fs/2),20/(Fs/2),1,60);
    [b,a] = butter(n, Wn);
    
    PitchF = [Pitch(:,1),filtfilt(b,a,Pitch(:,2:end))];
    RollF = [Pitch(:,1),filtfilt(b,a,Roll(:,2:end))];
    YawF = [Pitch(:,1),filtfilt(b,a,Yaw(:,2:end))];
    
    %% -- Resultant
    
    % Delsys
    Acc = [Pitch(:,1), sqrt(Pitch(:,2:end).^2 + Roll(:,2:end).^2 +Yaw(:,2:end).^2)];
    AccF = [Pitch(:,1), sqrt(PitchF(:,2:end).^2 + RollF(:,2:end).^2 +YawF(:,2:end).^2)];
    
    % Kistler
    % Acc1 = [Acc1, sqrt(Accx1.^2 + Accy1.^2 + Accz1.^2)];
    % Acc2 = Acc1;
    % Acc2 = [Acc2, sqrt(Accx2.^2 + Accy2.^2 + Accz2.^2)];
    
    %%
    
    instant = importdata('Instantes_gait.txt',',');
    [HeelContactInstant] = Instants(instant,File);
    HeelContact = floor(HeelContactInstant*Fs);
    
    Var = {'AccF','PitchF','RollF','YawF'};%,'AccF.^2'};
    Sensors = {'ShankR','ShankL'};
    
    for j = 1 : length(Sensors)
        first(j) = HeelContact(j,1);
        last(j) =  HeelContact(j,1)+floor(Fs);
        
        for i = 1: length(Var)
            
                        
            p(:,:,k+(i-1)*length(Trial)) = [p,eval([Var{i},'(first:last,eval(Sensors{j}))'])];
            
            stimulWin = exp(-0.5*((Acc(:,1)-(HeelContactInstant(1,2)-50e-3))/(50e-3/3)).^2);
            y(:,k+(i-1)*length(Trial)) = stimulWin(first:last);
            
%             p(:,:,k+length(Trial)) = [AccF(first:last,ShankR),...
%                 PitchF(first:last,ShankL),...
%                 RollF(HeelContact(1,1):HeelContact(1,1)+floor(Fs),ShankL),...
%                 YawF(HeelContact(1,1):HeelContact(1,1)+floor(Fs),ShankL),...
%                 AccF(HeelContact(1,1):HeelContact(1,1)+floor(Fs),ShankL).^2];
%             
%             stimulWin = exp(-0.5*((Acc(:,1)-(HeelContactInstant(2,2)-50e-3))/(50e-3/3)).^2);
%             y(:,k+length(Trial)) = stimulWin(HeelContact(2,1):HeelContact(2,1)+floor(Fs));
            
        end
    end
 end
% 
% beta = ols(p,y)
% %% -- Plots


figure;
subplot(3,1,1); plot(Fy(:,1), Fy(:,2:end));
legend({'1','2','3','4','5','6','7'})

%Delsys
subplot(3,1,2); plot(Acc(:,1), Acc(:,ShankL),AccF(:,1), AccF(:,ShankL)); ylabel('Shank L')
title('Resultant'); legend('Raw','Filtered')
subplot(3,1,3); plot(Acc(:,1), Acc(:,ShankR),AccF(:,1), AccF(:,ShankR)); ylabel('Shank R')
% Kistler
%subplot(3,1,2); plot(Acc1(:,1),Acc1(:,ShankL)); ylabel('Shank L')
%title('Resultant')
%subplot(3,1,3); plot(Acc2(:,1),Acc2(:,ShankL)); ylabel('Shank R')
ylim([-2 2]);
PlotInstants( instant, File )

%--
figure;
subplot(3,1,1)
plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})

subplot(3,1,2); plot(Acc(:,1), Pitch(:,ShankL),Acc(:,1), PitchF(:,ShankL)); ylabel('Shank L')
title('Pitch'); legend('Raw','Filtered')
subplot(3,1,3); plot(Acc(:,1), Pitch(:,ShankR),AccF(:,1), PitchF(:,ShankR)); ylabel('Shank R')
ylim([-2 2]);
PlotInstants( instant, File )
%--
figure;
subplot(3,1,1)
plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
subplot(3,1,2); plot(Acc(:,1), Roll(:,ShankL), Acc(:,1), RollF(:,ShankL));  ylabel('Shank L')
title('Roll'); legend('Raw','Filtered')
subplot(3,1,3); plot(Acc(:,1), Roll(:,ShankR), AccF(:,1), RollF(:,ShankR)); ylabel('Shank R')
ylim([-2 2]);
PlotInstants( instant, File )
%--
figure;
subplot(3,1,1)
plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
subplot(3,1,2); plot(Acc(:,1), Yaw(:,ShankL), Acc(:,1), YawF(:,ShankL));
ylabel('Shank L'); ylim([-2 2]);
title('Yaw'); legend('Raw','Filtered')
subplot(3,1,3); plot(Acc(:,1), Yaw(:,ShankR), AccF(:,1), YawF(:,ShankR));
ylabel('Shank R');ylim([-2 2]);
PlotInstants( instant, File )

%--
X = Acc; XF = AccF;
X(:,2:end) = Acc(:,2:end).^2 - Acc(:,2:end);
XF(:,2:end) = AccF(:,2:end).^2 - AccF(:,2:end);

figure;
subplot(3,1,1)
plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
subplot(3,1,2); plot(Acc(:,1), X(:,ShankL), Acc(:,1), XF(:,ShankL));
ylabel('Shank L');
title('Yaw'); legend('Raw','Filtered')
subplot(3,1,3); plot(Acc(:,1), X(:,ShankR), Acc(:,1), XF(:,ShankR));
ylabel('Shank R');ylim([-1 10]);
PlotInstants( instant, File )

%%
