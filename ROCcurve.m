function [Sensib,Specif] = ROCcurve(ResultsStruct)

Sensib = NaN (length(ResultsStruct),1);
Specif = NaN (length(ResultsStruct),1);

for i = 1: length(ResultsStruct)
    
    Sensib(i) = ResultsStruct(i).TP/...
        (ResultsStruct(i).TP + ResultsStruct(i).FN);
    
    Specif(i) = ResultsStruct(i).TN/...
        (ResultsStruct(i).FP + ResultsStruct(i).TN);
end


plot(1-Specif,Sensib,'.','MarkerSize',12)
line([0 max(1-Specif)], [0 max(Sensib)],'LineStyle','--','color','k')
xlabel('1-Specif');
ylabel('Sensib');

end