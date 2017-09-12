function [ResultsindTrs,ResultsCombinatorics, indTr] = combinatorics1_teste(Var,Name,Sensors,ToeOffInd,p,y,pct, winPts,Fy,HeelStrike,ToeOff,FsFP)

% -- Select trials for trainning and test
[indTr,indTr] = PartData(ToeOffInd,16);

% -- For plotting Fy
HeelStrikeFP = ceil(HeelStrike*FsFP);
ToeOffFP = ceil(ToeOff*FsFP);

firstFP = HeelStrikeFP(indTr,1);
lastFP =  ToeOffInd(indTr,end);

% --

ToeOffInd = ToeOffInd./max(ToeOffInd,[],2);


n = 0;
for k = 9: length(Var)
    % -- Combinatorial analysis to find the best set of variables
    % to be used in linear combination
    CombinatoricsInd = nchoosek(1:1:length(Var),k);
    
    for kk = 1 : size (CombinatoricsInd,1)
        TPos = 0; TNeg = 0; FPos = 0; FNeg = 0;
        numTrials = length(indTr)*length(Sensors);
        
        Features = CombinatoricsInd(kk,:);
        pTr = p(:,Features,indTr);
        yTr = y(:,indTr);
        
        pTs = p(:,Features,indTr);
        
        beta = ols(pTr,yTr); % Coefs for linear combination
        while sum(isnan(beta(1,:)))~= 0
            for i = 1:size(beta,2)-sum(isnan(beta(1,:)))
                if isnan(beta(1,i))
                    beta(:,i)=[];
                end
            end
        end
        betaM = mean(beta,2);
        
        % --- Applying Linear Combination
        
        for j = 1 : length(indTr)
            % File = [Name,indTr{j}];
            
            for i = 1 : length(Sensors)
                
                TP = 0; FP = 0; TN = 0; FN = 0;
                n = n+1; disp(n)
                
                j+(i-1)*length(indTr) = 1;
                
                LinearCombination(:,j+(i-1)*length(indTr)) = pTs(:,:,j+(i-1)*length(indTr))*betaM;
                
                %% Plot Fy
                
                ind = 1;
                while isnan(lastFP(j))
                    lastFP(j) = ToeOffFP(j,end-1);
                    ind = ind + 1;
                end
                cycle1 = ((firstFP(j):1:lastFP(j))-firstFP(j))/(lastFP(j)-firstFP(j));
                cycle = 0:1/(2*winPts - 1):1;
                
                figure(j+(i-1)*length(indTr))
                subplot(2,1,1); plot(cycle1, Fy(firstFP(j):lastFP(j),2:end));
                %title([Name,' Trial ', i]);
                subplot(2,1,2); plot(cycle,LinearCombination(:,j+(i-1)*length(indTr)));
                ylim([min(LinearCombination(:,j+(i-1)*length(indTr)))*1.1 max(LinearCombination(:,j+(i-1)*length(indTr)))*1.1]);
                
                % --- Checking the combination's quality
                threshold = (max(LinearCombination(:,j+(i-1)*length(indTr))))*pct;
                [pks,locs] = findpeaks(LinearCombination(:,j+(i-1)*length(indTr)),cycle,'MinPeakHeight',threshold);
                
                %                 x1(1:length(locs))=-1;
                %                 x2(1:length(locs))=100;
                %
                %                 Line = line([locs locs], [x1 x2],'Linewidth',1,'Linestyle','--','Color',[0 0 0]);
                %                 set(Line,'Clipping','off')
                
                for z = 1:length(locs)
                    Line = line([locs(z) locs(z)], [-1 100],'Linewidth',1,'Linestyle','--','Color',[0 0 0]);
                    set(Line,'Clipping','off')
                    
                    for w = 1 : size(ToeOffInd,2)
                        ind = ToeOffInd(indTr(j),w);
                        
                        if isnan(ind)
                            
                        else
                            
                            A = rangesearch(locs(z),ind-0.05,0.05);
                            if isempty(A{1})
                                FN = FN+1;
                            else
                                TP = TP+1;
                            end
                            A = rangesearch(locs(z),(ind-0.1)/2,(ind-0.1)/2);
                            B = rangesearch(locs(z),(1+ind)/2,(1-ind)/2);
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
                        end
                    end
                end
                % --- save
                ResultsindTrs(n) = struct('Trial',{indTr(j)},'Sensor',{Sensors(i)},...
                    'Features',{Var(Features)},'Locs',{locs},'Threshold',...
                    {pct},'TP',{TP},'FP',{FP},'TN',{TN},'FN',{FN},'beta',{beta});
                % structSave = struct('teste',{File,k,Features,TP,FN,TN,FP,beta});
                
                % keyboard %breakpoint
            end
        end
        ResultsCombinatorics(n/numTrials) = struct('Trial',{Name},...
            'Features',{Var(Features)},'Threshold',{pct},'TP',{TPos},...
            'FP',{FPos},'TN',{TNeg},'FN',{FNeg},'beta',{beta});
    end
end


end

