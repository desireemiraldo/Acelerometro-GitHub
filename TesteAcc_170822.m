%% Gait Analysis (w/ accelerometer)
clear all; clc; close all

%% -- About the files

Path = '.\Coletas\';

% -- Settings

% names = {'Acc_170803_EGS_','Acc_170731_RNW_','Acc_170731_DCM_'};
% ext = {'-Delsys 1.csv','-Delsys 1.csv','-Delsys 1.csv'};
% trials = struct('S1',{{'1' '2' '3' '4' '5'}},'S2',{{'1' '2' '4' '5'}},...
%     'S3',{{'1' '2' '3' '4'}});
% shanksL = [2,5,5];
% shanksR = [3,4,4];

Name = 'Acc_170803_EGS_';
csv = '-Delsys 1.csv';
Trial = {'1' '2' '3' '4' '5'};
ShankL = 2;ShankR = 3;

% Name = 'Acc_170731_RNW_';
% csv = '-Delsys 1.csv';
% Trial = {'1' '2' '4' '5'};
% ShankL = 5;ShankR = 4;

% Name = 'Acc_170731_DCM_';
% csv = '-Delsys 1.csv';
% Trial = {'1' '2' '3' '4' };
% ShankL = 5;ShankR = 4;

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

Sensors = {'ShankR','ShankL'};

% -- Filter
[n,Wn] = buttord(10/(Fs/2),20/(Fs/2),1,60);
[b,a] = butter(n, Wn);

% -- Standard deviation to be used in Linear Combination
sd = 50e-3;


%% Load Data

% for Sub = 1: length(names)
%     
%     Name = names {Sub};
%     csv = ext{Sub};
%     Trial = eval(['trials.S',num2str(Sub)]);
%     ShankL = shanksL(Sub); ShankR = shanksR(Sub);
    
    % -- Initializing variables for linear combination
    p = zeros(ceil(Fs),length(Var),length(Sensors)*length(Trial));
    y = zeros(ceil(Fs),2*length(Trial));
    
    tic
    for j = 1: length(Trial)
        
        File = [Name,Trial{j}];
        FilePath = [Path,File];
        
        for i = 1:length(Signal)
            VarName = strrep(strrep(Signal{i},'IM ',''),' ','');
            eval([VarName '= ReadDelsys([FilePath,csv], ChannelType, Signal(i));']);
            temp = eval(VarName);
            %Filtering
            eval([VarName 'F = [temp(:,1),filtfilt(b,a,temp(:,2:end))];']);
        end
        
        EMG = ReadDelsys([FilePath,csv], 'IMEMG', 'EMG');
        
        Cte = ones(size(eval(VarName),1),size(eval(VarName),2));
        
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
        ACC = [ACCPitch(:,1), sqrt(ACCPitch(:,2:end).^2 + ACCRoll(:,2:end).^2 + ACCYaw(:,2:end).^2)];
        ACCF = [ACCPitch(:,1), sqrt(ACCPitchF(:,2:end).^2 + ACCRollF(:,2:end).^2 + ACCYawF(:,2:end).^2)];
        
        GYR = [GYRPitch(:,1), sqrt(GYRPitch(:,2:end).^2 + GYRRoll(:,2:end).^2 + GYRYaw(:,2:end).^2)];
        GYRF = [GYRPitch(:,1), sqrt(GYRPitchF(:,2:end).^2 + GYRRollF(:,2:end).^2 + GYRYawF(:,2:end).^2)];
        
        %% Linear combination of different variables apllied in one gait cycle
        %
        instant = importdata('Instantes_gait.txt',',');
        [HeelContact, ToeOff] = Instants(instant,File);
        HeelContactInd = floor(HeelContact*Fs);
        
        % --- Building inputs for Orthogonal Least Squares Algorithm
        % --- (ols.m) implemented by Renato Naville Watanabe
        for i = 1 : length(Sensors)
            first(i) = HeelContactInd(i,1);
            last(i) =  HeelContactInd(i,1)+floor(Fs);
           
            
            if first(i)==0 || last(i)==0
               
                
            else
                ToeOffInd(i,j+(i-1)*length(Trial)) = (floor(ToeOff(i)*Fs)-first(i))/(last(i)-first(i));
                for jj = 1 : length (Var)
                    p(:,jj,j+(i-1)*length(Trial)) = eval([Var{jj},'(first(i):last(i),eval(Sensors{i}))']);
                    
                    stimulWin = exp(-0.5*((ACC(:,1)-(ToeOff(i,1)- sd))/(sd/3)).^2);
                    
                    y(:,j+(i-1)*length(Trial)) = stimulWin(first(i):last(i));
                    
                end
            end
        end
%% Plot Fy

%         firstFP = zeros(length(Sensors),1);
%         lastFP = zeros(length(Sensors),1);
%         HeelContactInd = ceil(HeelContact*FsFP);
%         
%         for i = 1: length(Sensors)
%             firstFP(i) = HeelContactInd(i,1);
%             % lastFP(i) =  HeelContactInd(i,2);
%             lastFP(i) =  HeelContactInd(i,1)+ceil(FsFP);
%         end
%         for i = 1 : length(Sensors)
%             cycle1 = ((firstFP(i):1:lastFP(i))-firstFP(i))/(lastFP(i)-firstFP(i));
%             figure(j+(i-1)*length(Trial));
%             subplot(2,1,1); plot(cycle1, Fy(firstFP(i):lastFP(i),2:end));
%             title([Name,Trial{j},' (',Sensors{i},')']);
%         end
    end

    
    n = 0;
for k = 1: length(Var)
    % -- Combinatorial analysis to find the best set of variables
    % to be used in linear combination
    CombinatoricsInd = nchoosek(1:1:length(Var),k);
    
    for kk = 1 : size (CombinatoricsInd,1)
        TPos = 0; TNeg = 0; FPos = 0; FNeg = 0;
        
        Features = CombinatoricsInd(kk,:);        
        pp = p(:,Features,:);
        
        beta = ols(pp,y); % Coefs for linear combination
        betaM = mean(beta(:,1:length(Trial)),2);
        
        % --- Applying Linear Combination
        
        for j = 1 : length(Trial)
            File = [Name,Trial{j}];         
            
            for i = 1 : length(Sensors)
                TP = 0; FP = 0; TN = 0; FN = 0;
                n = n+1; disp(n)
                
                LinearCombination(:,j+(i-1)*length(Trial)) = pp(:,:,j+(i-1)*length(Trial))*betaM;
                
%                 cycle = ((first(i):1:last(i))-first(i))/(last(i)-first(i));
%                 figure(j+(i-1)*length(Trial))
%                 subplot(2,1,2); plot(cycle,LinearCombination(:,j+(i-1)*length(Trial)));
%                 ylim([min(LinearCombination(:,j+(i-1)*length(Trial)))*1.1 max(LinearCombination(:,j+(i-1)*length(Trial)))*1.1]);
                
                % --- Checking the combination's quality              
                threshold = (max(LinearCombination(:,j+(i-1)*length(Trial))))*0.75;
                [pks,locs] = findpeaks(LinearCombination(:,j+(i-1)*length(Trial)),Fs,'MinPeakHeight',threshold);
                
%                 Line = line([locs locs], [-1 100],'Linewidth',1,'Linestyle','--','Color',[0 0 0]);
%                 set(Line,'Clipping','off')
                
                
                A = rangesearch(locs,ToeOffInd(i,j+(i-1)*length(Trial))-0.05,0.05);
                if isempty(A{1})
                    FN = FN+1;
                else
                    TP = TP+1;
                end
                A = rangesearch(locs,(ToeOffInd(i,j+(i-1)*length(Trial))-0.1)/2,(ToeOffInd(i,j+(i-1)*length(Trial))-0.1)/2);
                B = rangesearch(locs,(1+ToeOffInd(i,j+(i-1)*length(Trial)))/2,(1-ToeOffInd(i,j+(i-1)*length(Trial)))/2);
                if isempty(A{1}) || isempty(B{1})
                    TN = TN+1;
                else
                    FP = FP+1;
                end
                
                %---
                TPos = TPos + TP;
                TNeg = TNeg + TN; 
                FPos = FPos + FP; 
                FNeg = FNeg + FN;
              
                
                % --- save
                ResultsTrials(n) = struct('Trial',{File},'Sensor',{Sensors(i)},...
                    'k',{k},'Features',{Var{Features}},'Locs',{locs},'TP',{TP},...
                    'FP',{FP},'TN',{TN},'FN',{FN},'beta',{beta});
                % structSave = struct('teste',{File,k,Features,TP,FN,TN,FP,beta});
                
                % keyboard %breakpoint
            end 
        end
        ResultsCombinatorics(n) = struct('Trials',{Name},'k',{k},'Features',...
            {Var{Features}},'Locs',{locs},'TP',{TPos},'FP',{FPos},'TN',{TNeg},...
            'FN',{FNeg},'beta',{beta});
    end
end
    
    save(['ResultTrials_',Name,'.mat'],'ResultsTrials')
    save(['ResultCombinatorics',Name,'.mat'],'ResultsTrials')
    toc
% end
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

%% --
% % X = ACC; XF = ACCF;
% % X(:,2:end) = ACC(:,2:end).^2 - ACC(:,2:end);
% % XF(:,2:end) = ACCF(:,2:end).^2 - ACCF(:,2:end);
% %
% % figure;
% % subplot(3,1,1)
% % plot(Fy(:,1), Fy(:,2:end)); legend({'1','2','3','4','5','6','7'})
% % subplot(3,1,2); plot(ACC(:,1), X(:,ShankL), ACC(:,1), XF(:,ShankL));
% % ylabel('Shank L');
% % title('Yaw'); legend('Raw','Filtered')
% % subplot(3,1,3); plot(ACC(:,1), X(:,ShankR), ACC(:,1), XF(:,ShankR));
% % ylabel('Shank R');ylim([-1 10]);
% % PlotInstants( instant, File )
