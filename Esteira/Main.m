%% Gait Analysis (w/ accelerometer)
clear all; clc; close all

%% -- About the files
% % Path = '.\Dropfoot\';
% % % -- Settings
% % Name = 'DCM\Descalco25_Dual_170831_1';
% % csv = '-Delsys 1.csv';

paths = {'.\Dropfoot\'};
folders = {'DCM\'};
names = {'Descalco25_Dual_170831_1'};
ext = {'-Delsys 1.csv'};
shanksL = 2;
shanksR = 1;
sensor = {{'ShankR'}};

Win = {'Gauss'}; %,'Rect'};

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
        Folder = folders{Sub};
        Name = names {Sub};
        csv = ext{Sub};
        ShankR = shanksR(Sub); %ShankL = shanksL(Sub);
        Sensors = sensor{Sub};
        
        
        FilePath = [Path,Folder,Name];
        
        %% Loading Data
        
        [Trials,deltaT,TrLabel] = ReadDelsysWindows([FilePath,csv], winPts);
        
        instant = importdata('Instantes_gait1.txt',';');
        
        %[NewInstant] = ReshapeInstants(deltaT, instant,[Folder,Name]);
        [ToeOff, HeelStrike] = ReshapeInstants1(deltaT, instant,[Folder,Name]);
        
        
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
        
        
        
        %% Linear combination of different variables apllied in one trial
        %
        %[HeelStrike, ToeOff] = Instants1(instant);
        HeelStrikeInd = floor(HeelStrike*Fs);
        ToeOffInd = floor(ToeOff*Fs);
        
        % --- Building inputs for Orthogonal Least Squares Algorithm
        % --- (ols.m) implemented by Renato Naville Watanabe
        
        % -- Initializing variables for linear combination
        
        p = NaN(length(ACC),length(Var),length(Sensors)*size(Trials,3));
        % y = NaN(ceil(Fs),length(Sensors)*size(Trials,3));
        
        for j = 1 : size(Trials,3)
            for i = 1 : length(Sensors)
                first(j) = min([HeelStrike(j,:),ToeOff(j,:)],[],2);
                last(j) = max([HeelStrike(j,:),ToeOff(j,:)],[],2);
                
                tempTO = (ToeOff(j,:));
                tempTO(isnan(tempTO))=[];
                if strcmp(Win(w),'Gauss')
                    stimulWin = sum(exp(-0.5*((ACC(:,1,j) - (tempTO - sd))/(sd/3)).^2),2);
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
        end
        
        %% --- TESTANDO AS JANELAS DE ESTÍMULO
        timeACC(:,:,1) = ACC(:,1,:);
        %         figure;
        %         plot(Fy(:,1),Fy(:,3),'k')
        %         hold on
        %         plot(timeACC,1100*y)
        
        %% --  -- Select trials for trainning and test
        [indTr,indTs] = PartData(ToeOff,16);
        %% --- Combinatorial Analysis
        
        t = 0; c = 0;
        for pct = 0: 0.01 : 1

            [ResultsTrials,ResultsCombinatorics] = combinatorics1(Var,[Folder,Name],Sensors,p,y,pct,Fy,ToeOff,timeACC,indTr,indTs);
            
            RT(t+1:t+length(ResultsTrials)) = ResultsTrials;
            RC(c+1:c+length(ResultsCombinatorics)) = ResultsCombinatorics;
            
            t = length(RT);
            c = length(RC);
            
        end
        
        save(['RT_',Name,'.mat'],'RT')
        save(['RC_',Name,'.mat'],'RC')
        RT(:) = []; RC(:) = [];
        
        toc
        
    end
end


% A = struct('Trial',{0},...
% 'Features',{0},'Threshold',{0},'TP',{0},...
% 'FP',{0},'TN',{0},'FN',{0},'beta',{0});
% repmat(A,51100*16,1)




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
