function [] = SavingCombinatorics(Name, Trial, Sensors, ResultsStruct)

n = length(Sensors)*length(Trial);

for i = 1:(length(ResultsStruct)/n)
    TPos = 0; FPos = 0; TNeg = 0; FNeg = 0;
    
    for j = (i-1)*n+1: i*n
        TPos = TPos + ResultsStruct(j).TP;
        FPos = FPos + ResultsStruct(j).FP;
        TNeg = TNeg + ResultsStruct(j).TN;
        FNeg = FNeg + ResultsStruct(j).FN;
    end
    
    ResultsCombinatorics(i) = struct('Trials',{ResultsStruct(i*n).Trial},'k',...
        {ResultsStruct(i*n).k},'Features',{ResultsStruct(i*n).Features},'TP',...
        {TPos},'FP',{FPos},'TN',{TNeg},'FN',{FNeg},'beta',{ResultsStruct(i).beta});
end

save(['ResultCombinatorics_',Name,'.mat'],'ResultsCombinatorics')
end