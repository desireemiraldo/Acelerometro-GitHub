function [indTr,indTs] = PartData(instant,Ntest)

N = length(instant.data)/2;


[Train, Test] = crossvalind('LeaveMOut', N, Ntest);

[indTr val] = find(Train==1);

[indTs val] = find(Test==1);

end


