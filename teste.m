clc;clear all; close all

% load('ResultCombinatorics_Acc_170803_EGS_.mat');
load('ResultCombinatorics_Acc_170731_RNW_.mat')
% load('ResultCombinatorics_Acc_170731_DCM_.mat')

[Sensib,Specif] = ROCcurve(ResultsCombinatorics);


%%---
D = sqrt((1-Sensib).^2 + (1-Specif).^2);
k=1;
for i = 1: length(D)
    if D(i)== min(D)
        Features{k,1} = ResultsCombinatorics(i).Features;
        k = k+1;
    end
end

