function [ResultsTrials,ResultsCombinatorics] = combinatorics(Var,Name,Trial,Sensors,ToeOffInd,Fs,p,y)

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


end

