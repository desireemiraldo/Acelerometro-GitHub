%% Gait Analysis (w/ accelerometer)
clear all; clc; close all

%% -- About the files
% % Path = '.\Dropfoot\';
% % % -- Settings
% % Name = 'DCM\Descalco25_Dual_170831_1';
% % csv = '-Delsys 1.csv';

paths = {'.\Dropfoot\'};
names = {'DCM\Descalco25_Dual_170831_1'};
ext = {'-Delsys 1.csv'};
shanksL = 2;
shanksR = 1;
sensor = {{'ShankR'}};

Win = {'Gauss','Rect'};

%% -- Initializing some variables

% -- Sample Frequencies
Fs = 148.39; % Delsys sensor
FsFP = 300;  % force plates

% --
ChannelType = 'AUX';
Signal = {'IM ACC Pitch', 'IM ACC Roll', 'IM ACC Yaw',...
    'IM GYR Pitch', 'IM GYR Roll', 'IM GYR Yaw'};

Var = {'ACCF','ACCPitchF','ACCRollF','ACCYawF',...
    'GYRF','GYRPitchF','GYRRollF','GYRYawF','Cte'};%,'AccF.^2'}; %Used in linear combination


% -- Filter
[t,Wn] = buttord(10/(Fs/2),20/(Fs/2),1,60);
[b,a] = butter(t, Wn);

% -- Standard deviation to be used in Linear Combination
sd = 50e-3;

% --
winPts = 202;

for Sub = 1: length(names)
    for w = 1: length(Win)
        
        Path = paths{Sub};
        Name = names {Sub};
        csv = ext{Sub};
        ShankR = shanksR(Sub); %ShankL = shanksL(Sub);
        Sensors = sensor{Sub};
        
        
        FilePath = [Path,Name];
        
        %% Loading Data
        
        [Trials,deltaT,TrLabel] = ReadDelsysWindows([FilePath,csv], winPts);
        
        instant = importdata('Instantes_gait1.txt',';');
        
        [instant] = ReshapeInstants(deltaT, instant,Name);
        
        % Select trials for trainning and test
        [indTr,indTs] = PartData(instant,16);
        
        
        tic
        
        for i = 1:length(Signal)
            VarName = strrep(strrep(Signal{i},'IM ',''),' ','');
            eval([VarName '= SelectVar(Trials, ChannelType, Signal(i));']);
            temp = eval(VarName);
            %Filtering
            eval([VarName 'F = [temp(:,1,:),filtfilt(b,a,temp(:,2:end,:))];']);
        end
        
        %EMG = SelectVar(Trials, 'IMEMG', 'EMG');
        
        Cte = ones(size(eval(VarName),1),size(eval(VarName),2), size(eval(VarName),3));
        
        % Cortex data
        Forces = importdata([FilePath,'.forces']);
        Fy = (Forces.data(:,1) -1)/FsFP;
        for i =1 : length(Forces.colheaders)
            if strcmp(Forces.colheaders{i}(1:2), 'FY')
                Fy = [Fy, Forces.data(:,i)];
            end
        end
        
        
        Fy(:,6:8) = Fy(:,6:8)/10;
        
        %% -- Resultants
        
        % Delsys
        ACC = [ACCPitch(:,1,:), sqrt(ACCPitch(:,2:end,:).^2 + ACCRoll(:,2:end,:).^2 + ACCYaw(:,2:end,:).^2)];
        ACCF = [ACCPitch(:,1,:), sqrt(ACCPitchF(:,2:end,:).^2 + ACCRollF(:,2:end,:).^2 + ACCYawF(:,2:end,:).^2)];
        
        GYR = [GYRPitch(:,1,:), sqrt(GYRPitch(:,2:end,:).^2 + GYRRoll(:,2:end,:).^2 + GYRYaw(:,2:end,:).^2)];
        GYRF = [GYRPitch(:,1,:), sqrt(GYRPitchF(:,2:end,:).^2 + GYRRollF(:,2:end,:).^2 + GYRYawF(:,2:end,:).^2)];
        
        
        
        %% Linear combination of different variables apllied in one gait cycle
        %
        [HeelStrike, ToeOff] = Instants1(instant);
        HeelStrikeInd = floor(HeelStrike*Fs);
        ToeOffInd = floor(ToeOff*Fs);
        
        % --- Building inputs for Orthogonal Least Squares Algorithm
        % --- (ols.m) implemented by Renato Naville Watanabe
        
        % -- Initializing variables for linear combination
        
        p = NaN(length(ACC),length(Var),length(Sensors)*size(Trials,3));
        % y = NaN(ceil(Fs),length(Sensors)*size(Trials,3));
        
        for j = 1 : size(Trials,3)
            for i = 1 : length(Sensors)
                first(i) = HeelStrikeInd(i,1);
                last(i) =  ToeOffInd(i,end);
                
                ind = 1;
                while isnan(last(i))
                    last(i) = ToeOffInd(i,end-1);
                    ind = ind + 1;
                end
                
                if strcmp(Win(w),'Gauss')
                    stimulWin = sum(exp(-0.5*((ACC(:,1)-(ToeOff(j,:)- sd))/(sd/3)).^2),2);
                end
                
                if strcmp(Win(w),'Rect')
                    stimulWin = zeros(length(ACC),1);
                    for kk = 1: size(ToeOff,2)
                        stimulWin(floor((ToeOff(j,kk)-2*sd)*Fs):floor(ToeOff(j,kk)*Fs)) = 1;
                    end
                end
                y(:,j) = stimulWin;
                
                for jj = 1 : length (Var)
                    p(:,jj,:) = eval([Var{jj},'(:,eval(Sensors{i})+1,:)']);
                    
                end
            end
            
            
            %% Plot Fy
            
            firstFP = zeros(length(Sensors),1);
            lastFP = zeros(length(Sensors),1);
            HeelStrikeInd = ceil(HeelStrike*FsFP);
            ToeOffInd = ceil(ToeOff*FsFP);
            
            for i = 1: length(Sensors)
                first(i) = HeelStrikeInd(i,1);
                last(i) =  ToeOffInd(i,end);
                
                ind = 1;
                while isnan(last(i))
                    last(i) = ToeOffInd(i,end-1);
                    ind = ind + 1;
                end
            end
            %     for i = 1 : length(Sensors)
            %         cycle1 = ((firstFP(i):1:lastFP(i))-firstFP(i))/(lastFP(i)-firstFP(i));
            %         figure(j+(i-1)*length(Trial));
            %         subplot(2,1,1);
            %         plot(cycle1, Fy(firstFP(i):lastFP(i),2:end));
            %         title([Name,Trial{j},' (',Sensors{i},')']);
            %
            %     end
        end
        
        %% --- Combinatorial Analysis
        
        t = 0; c = 0;
        for pct = 0: 0.01 : 1
            pct= 0.75; %% APAGAR
            [ResultsTrials,ResultsCombinatorics] = combinatorics1(Var,Name,indTr,Sensors,ToeOffInd,Fs,p,y,pct);
            
            RT(t+1:t+length(ResultsTrials)) = ResultsTrials;
            RC(c+1:c+length(ResultsCombinatorics)) = ResultsCombinatorics;
            
            t = length(RT);
            c = length(RC);
            
        end
        
        save(['RTsquare_',Name,'.mat'],'RT')
        save(['RCsquare_',Name,'.mat'],'RC')
        RT(:) = []; RC(:) = [];
        
        toc
        
    end
end
%% -- Plots

% % figure;
% % subplot(3,1,1); plot(Fy(:,1), Fy(:,2:end));
% % legend({'1','2','3','4','5','6','7'})
% %
% % %Delsys
% % subplot(3,1,2); plot(ACC(:,1), ACC(:,ShankL),ACCF(:,1), ACCF(:,ShankL)); ylabel('Shank L')
% % title('Resultant'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACC(:,ShankR),ACCF(:,1), ACCF(:,ShankR)); ylabel('Shank R')
% % % Kistler
% % %subplot(3,1,2); plot(ACC1(:,1),ACC1(:,ShankL)); ylabel('Shank L')
% % %title('Resultant')
% % %subplot(3,1,3); plot(ACC2(:,1),ACC2(:,ShankL)); ylabel('Shank R')
% % ylim([-2 2]);
% % PlotInstants( instant, File )
% %
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% %
% % subplot(3,1,2); plot(ACC(:,1), ACCPitch(:,ShankL),ACC(:,1), ACCPitchF(:,ShankL)); ylabel('Shank L')
% % title('Pitch'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCPitch(:,ShankR),ACCF(:,1), ACCPitchF(:,ShankR)); ylabel('Shank R')
% % ylim([min(ACCPitch(:,ShankL))-1, max(ACCPitch(:,ShankL))+1]);
% % PlotInstants( instant, File )
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(ACC(:,1), ACCRoll(:,ShankL), ACC(:,1), ACCRollF(:,ShankL));  ylabel('Shank L')
% % title('Roll'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCRoll(:,ShankR), ACCF(:,1), ACCRollF(:,ShankR)); ylabel('Shank R')
% % ylim([min(ACCRoll(:,ShankL))-1, max(ACCRoll(:,ShankL))+1]);
% % PlotInstants( instant, File )
% % %--
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(ACC(:,1), ACCYaw(:,ShankL), ACC(:,1), ACCYawF(:,ShankL)); ylabel('Shank L');
% % title('Yaw'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), ACCYaw(:,ShankR), ACCF(:,1), ACCYawF(:,ShankR));ylabel('Shank R');
% % ylim([min(ACCYaw(:,ShankL))-1, max(ACCYaw(:,ShankL))+1]);
% % PlotInstants( instant, File )
