function [ResultsindTrs,ResultsCombinatorics] = combinatorics(Var,Name,indTr,Sensors,ToeOffInd,Fs,p,y,pct)

ToeOffInd = ToeOffInd./max(ToeOffInd,[],2);

n = 0;
for k = 1: length(Var)
    % -- Combinatorial analysis to find the best set of variables
    % to be used in linear combination
    CombinatoricsInd = nchoosek(1:1:length(Var),k);
    
    for kk = 1 : size (CombinatoricsInd,1)
        TPos = 0; TNeg = 0; FPos = 0; FNeg = 0;
        numTrials = length(indTr)*length(Sensors);
        
        Features = CombinatoricsInd(kk,:);
        pp = p(:,Features,:);
        
        beta = ols(pp,y); % Coefs for linear combination
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
                
                
                LinearCombination(:,j+(i-1)*length(indTr)) = pp(:,:,j+(i-1)*length(indTr))*betaM;
                
                cycle = 0:1/403:1;
                figure(j+(i-1)*length(indTr))
                subplot(2,1,2); plot(cycle,LinearCombination(:,j+(i-1)*length(indTr)));
                ylim([min(LinearCombination(:,j+(i-1)*length(indTr)))*1.1 max(LinearCombination(:,j+(i-1)*length(indTr)))*1.1]);
                
                % --- Checking the combination's quality
                threshold = (max(LinearCombination(:,j+(i-1)*length(indTr))))*pct;
                [pks,locs] = findpeaks(LinearCombination(:,j+(i-1)*length(indTr)),Fs,'MinPeakHeight',threshold);
                
                Line = line([locs locs], [-1 100],'Linewidth',1,'Linestyle','--','Color',[0 0 0]);
                set(Line,'Clipping','off')
                
                for w = 1 : size(ToeOffInd,2)
                    ind = ToeOffInd(indTr(j),w);
                    
                    if isnan(ind)
                        
                    else
                        
                        A = rangesearch(locs,ind-0.05,0.05);
                        if isempty(A{1})
                            FN = FN+1;
                        else
                            TP = TP+1;
                        end
                        A = rangesearch(locs,(ind-0.1)/2,(ind-0.1)/2);
                        B = rangesearch(locs,(1+ind)/2,(1-ind)/2);
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

