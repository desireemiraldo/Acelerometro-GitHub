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

Path = 'C:\Users\BMClab\Downloads\Desiree\Acelerometro\Acelerometro-GitHub\Coletas\';

%Delsys data

Name = 'Acc_170803_EGS_';
csv = '-Delsys 1.csv';
Trial = {'1' '2' '3' '4' '5'};
ShankL = 2;ShankR = 3;

% Name = 'Acc_170731_RNW_';
% csv = '-Delsys 1.csv';
% Trial = {'1' '2' '4' '5'};
%ShankL = 5;ShankR = 4;

% Name = 'Acc_170731_DCM_';
% csv = '-Delsys 1.csv';
% Trial = {'1' '2' '3' '4' };
% ShankL = 5;ShankR = 4;

Fs = 148.39; % Delsys sensor
FsFP = 300;  % force plates


ChannelType = 'AUX';
Signal = {'IM ACC Pitch', 'IM ACC Roll', 'IM ACC Yaw',...
    'IM GYR Pitch', 'IM GYR Roll', 'IM GYR Yaw'};

Var = {'ACCF','ACCPitchF','ACCRollF','ACCYawF',...
    'GYR','GYRPitchF','GYRRollF','GYRYawF'};%,'AccF.^2'}; %Used in linear combination

Sensors = {'ShankR','ShankL'};

% -- Filter
[n,Wn] = buttord(10/(Fs/2),20/(Fs/2),1,60);
[b,a] = butter(n, Wn);
% -- Linear combination
p = zeros(ceil(Fs),length(Var),length(Sensors)*length(Trial));
y = zeros(ceil(Fs),2*length(Trial));


for k = 1: length(Trial)
    File = [Name,Trial{k}];
    FilePath = [Path,File];
    
    for i = 1:length(Signal)
        VarName = strrep(strrep(Signal{i},'IM ',''),' ','');
        eval([VarName '= ReadDelsys([FilePath,csv], ChannelType, Signal(i));']);
        temp = eval(VarName);
        %Filtering
        eval([VarName 'F' ' = [temp(:,1),filtfilt(b,a,temp(:,2:end))];']);
    end
    
    EMG = ReadDelsys([FilePath,csv], 'IMEMG', 'EMG');
    
    cte = ones(size(eval(VarName),1),size(eval(VarName),2));
    
    % Cortex data
    Forces = importdata([FilePath,'.forces']);
    Fy = (Forces.data(:,1) -1)/300;
    for i =1 : length(Forces.colheaders)
        if strcmp(Forces.colheaders{i}(1:2), 'FY')
            Fy = [Fy, Forces.data(:,i)];
        end
    end
    
    
    %% -- Resultant
    
    % Delsys
    ACC = [ACCPitch(:,1), sqrt(ACCPitch(:,2:end).^2 + ACCRoll(:,2:end).^2 + ACCYaw(:,2:end).^2)];
    ACCF = [ACCPitch(:,1), sqrt(ACCPitchF(:,2:end).^2 + ACCRollF(:,2:end).^2 + ACCYawF(:,2:end).^2)];
    
    GYR = [GYRPitch(:,1), sqrt(GYRPitch(:,2:end).^2 + GYRRoll(:,2:end).^2 + GYRYaw(:,2:end).^2)];
    GYRF = [GYRPitch(:,1), sqrt(GYRPitchF(:,2:end).^2 + GYRRollF(:,2:end).^2 + GYRYawF(:,2:end).^2)];
    
    % Kistler
    % Acc1 = [Acc1, sqrt(Accx1.^2 + Accy1.^2 + Accz1.^2)];
    % Acc2 = Acc1;
    % Acc2 = [Acc2, sqrt(Accx2.^2 + Accy2.^2 + Accz2.^2)];
    
    %% Linear combination of different variables apllied in one gait cycle
    %  Gait cycle = HeelContact + 1 sec
    
    instant = importdata('Instantes_gait.txt',',');
    [HeelContactInstant] = Instants(instant,File);
    HeelContact = floor(HeelContactInstant*Fs);
    
    % --- Building inputs for Orthogonal Least Squares Algorithm
    % --- (ols.m) implemented by Renato Naville Watanabe
    for j = 1 : length(Sensors)
        first(j) = HeelContact(j,1);
        last(j) =  HeelContact(j,1)+floor(Fs);
        
        %% --- ISSO NÃO DEVE SER REPETIDO A CADA TRIAL
        for i = 1: length(Var)
            
            Combinations = nchoosek(Var,i);
            
            for kk = 1 : size (Combinations,1)
                
                Features = Combinations(kk,:);
                
         %% --- ISSO NÃO DEVE SER REPETIDO A CADA TRIAL
                
                p(:,i,k+(j-1)*length(Trial)) = eval([Features{i},'(first:last,eval(Sensors{j}))']);
                
                stimulWin = exp(-0.5*((ACC(:,1)-(HeelContactInstant(1,2)-50e-3))/(50e-3/3)).^2);
                
                y(:,k+(j-1)*length(Trial)) = stimulWin(first:last);
            end
            
        end
    end
    
    %     HeelContact = ceil(HeelContactInstant*FsFP);
    %     for j = 1: length(Sensors)
    %         first1(j) = HeelContact(j,1);
    %         last1(j) =  HeelContact(j,1)+floor(FsFP);
    %     end
    %     figure(k);
    %     subplot(2,1,1); plot(Fy(first1(1):last1(1),1), Fy(first1(1):last1(1),2:end));
end

beta = ols(p,y); % Coefs for linear combination

% --- Applying LC
LinearCombination = zeros(ceil(Fs),1);
for i = 1:length(Var)
    for j = 1 : length(Sensors)*length(Trial)
        LinearCombination = LinearCombination + beta(i,j)*p(:,i,j);
    end
end

%%

% for i = 1:length(Trial)
% figure(i)
% subplot(2,1,2); plot(LinearCombination);
% end

%% -- Plots
% % % %
% % figure;
% % subplot(3,1,1); plot(Fy(:,1), Fy(:,2:end));
% % legend({'1','2','3','4','5','6','7'})
% %
% % %Delsys
% % subplot(3,1,2); plot(Acc(:,1), Acc(:,ShankL),AccF(:,1), AccF(:,ShankL)); ylabel('Shank L')
% % title('Resultant'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(Acc(:,1), Acc(:,ShankR),AccF(:,1), AccF(:,ShankR)); ylabel('Shank R')
% % % Kistler
% % %subplot(3,1,2); plot(Acc1(:,1),Acc1(:,ShankL)); ylabel('Shank L')
% % %title('Resultant')
% % %subplot(3,1,3); plot(Acc2(:,1),Acc2(:,ShankL)); ylabel('Shank R')
% % ylim([-2 2]);
% % PlotInstants( instant, File )
% %
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% %
% % subplot(3,1,2); plot(Acc(:,1), Pitch(:,ShankL),Acc(:,1), PitchF(:,ShankL)); ylabel('Shank L')
% % title('Pitch'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(Acc(:,1), Pitch(:,ShankR),AccF(:,1), PitchF(:,ShankR)); ylabel('Shank R')
% % ylim([-2 2]);
% % PlotInstants( instant, File )
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(Acc(:,1), Roll(:,ShankL), Acc(:,1), RollF(:,ShankL));  ylabel('Shank L')
% % title('Roll'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(Acc(:,1), Roll(:,ShankR), AccF(:,1), RollF(:,ShankR)); ylabel('Shank R')
% % ylim([-2 2]);
% % PlotInstants( instant, File )
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(Acc(:,1), Yaw(:,ShankL), Acc(:,1), YawF(:,ShankL));
% % ylabel('Shank L'); ylim([-2 2]);
% % title('Yaw'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(Acc(:,1), Yaw(:,ShankR), AccF(:,1), YawF(:,ShankR));
% % ylabel('Shank R');ylim([-2 2]);
% % PlotInstants( instant, File )
% %
% % %--
% % X = Acc; XF = AccF;
% % X(:,2:end) = Acc(:,2:end).^2 - Acc(:,2:end);
% % XF(:,2:end) = AccF(:,2:end).^2 - AccF(:,2:end);
% %
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(Acc(:,1), X(:,ShankL), Acc(:,1), XF(:,ShankL));
% % ylabel('Shank L');
% % title('Yaw'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(Acc(:,1), X(:,ShankR), Acc(:,1), XF(:,ShankR));
% % ylabel('Shank R');ylim([-1 10]);
% % PlotInstants( instant, File )