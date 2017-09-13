function [ResultsindTrs,ResultsCombinatorics] = combinatorics1(Var,Name,Sensors,p,y,pct,ForceY,ToeOff,timeACC,indTr,indTs)


%% -- Dividing forces
trialsNumber = length(ToeOff);
ind = 0;
while rem(size(ForceY(1:end-ind,1)),trialsNumber)~=0
    ind = ind+1;
end

[row,col] = size(ForceY(1:end-ind,:));
len = row/trialsNumber;

Fy = zeros(len,col,trialsNumber);

for i = 1: trialsNumber
    Fy(:,:,i) = ForceY((i-1)*len +1 : i*len,:);
end


%% -- Select trials for trainning and test
% [indTr,indTs] = PartData(ToeOff,16);

%% --

n = 0;
for k = 1: length(Var)
    % -- Combinatorial analysis to find the best set of variables
    % to be used in linear combination
    CombinatoricsInd = nchoosek(1:1:length(Var),k);
    
    for kk = 1 : size (CombinatoricsInd,1)
        TPos = 0; TNeg = 0; FPos = 0; FNeg = 0;
        
        Features = CombinatoricsInd(kk,:);
        pTr = p(:,Features,indTr);
        yTr = y(:,indTr);
        
        pTs = p(:,Features,indTs);
        
        beta = ols(pTr,yTr); % Coefs for linear combination
        betaM = mean(beta,2);
        
        % --- Applying Linear Combination
        
        for j = 1:length(indTs)
            for i = 1 : length(Sensors)
                
                index = indTs(j)+(i-1)*length(indTs);
                
                TP = 0; FP = 0; TN = 0; FN = 0;
                n = n+1; disp(n)
                
                LinearCombination(:,index) = pTs(:,:,(j)+(i-1)*length(indTs))*betaM;
                
                %% Plot Fy
                
%                 figure(index)
%                 subplot(2,1,1); plot(Fy(:,1,index), Fy(:,2:end,index),'k');
%                 hold on;
%                 plot(timeACC(:,index),1100*y(:,index),'r')
%                 title([' Trial ', num2str(index)]);
%                 
%                 subplot(2,1,2); plot(timeACC(:,index),LinearCombination(:,index));
%                 ylim([min(LinearCombination(:,index))*1.1 max(LinearCombination(:,index))*1.1]);
                
                %%  --- Checking the combination's quality
                threshold = (max(LinearCombination(:,index)))*pct;
                [pks,locs] = findpeaks(LinearCombination(:,index),timeACC(:,index),'MinPeakHeight',threshold);
                
%                 for z = 1:length(locs)
%                     Line = line([locs(z) locs(z)], [-1 100],'Linewidth',1,'Linestyle','--','Color',[0 0 0]);
%                     set(Line,'Clipping','off')
%                 end
                
                ind = [Fy(1,1,index), ToeOff(index,:),Fy(end,1,index)+0.1];
                ind(isnan(ind))=[];
                
                for w = 2: length(ind)
                    center = ((ind(w)-0.1) + ind(w-1))/2;
                    radius = abs(ind(w)-0.1 - center);
                    
                    if w~= length(ind)
                        A = rangesearch(locs,ind(w)-0.05,0.05);
                        if isempty(A{1})
                            FN = FN+1;
                        else
                            TP = TP+1;
                        end
                    end
                    
                    if center >= ind(1)
                        A = rangesearch(locs,center,radius);
                        if isempty(A{1})
                            TN = TN+1;
                        else
                            FP = FP+1;
                        end
                    end
                end
                
                %---
                TPos = TPos + TP;
                TNeg = TNeg + TN;
                FPos = FPos + FP;
                FNeg = FNeg + FN;
                
                % --- save
                ResultsindTrs(n) = struct('Trial',{indTs(j)},'Sensor',{Sensors(i)},...
                    'Features',{Var(Features)},'Locs',{locs},'Threshold',...
                    {pct},'TP',{TP},'FP',{FP},'TN',{TN},'FN',{FN},'beta',{beta});
                % structSave = struct('teste',{File,k,Features,TP,FN,TN,FP,beta});
                
                % keyboard %breakpoint
            end
        end
        ResultsCombinatorics(n/length(indTs)) = struct('Trial',{Name},...
            'Features',{Var(Features)},'Threshold',{pct},'TP',{TPos},...
            'FP',{FPos},'TN',{TNeg},'FN',{FNeg},'beta',{beta});
    end
end


end

