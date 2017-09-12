function [ ] = TestandoJanelas( HeelStrike, ToeOff, FsFP, Fy, y, Fs, winPts )

y1 = reshape(1100*y,[],1);
plot(0:1/FsFP:59.9967,Fy(:,3),'b'); hold on
plot(0:1/Fs:59.8895,y1,'r')


instant = importdata('Instantes_gait1.txt',';');
[NewInstant] = ReshapeInstants(deltaT, instant,Name);
NewInstant = NewInstant + deltaT(:,1);

first = min(NewInstant,[],2);
last = max(NewInstant,[],2);

HeelStrikeFP = ceil(HeelStrike*FsFP);
ToeOffFP = ceil(ToeOff*FsFP);

firstFP = min([HeelStrikeFP,ToeOffFP],[],2);
lastFP = max([HeelStrikeFP,ToeOffFP],[],2);

ToeOffInd = floor(ToeOff*FsFP);
ToeOffInd = ToeOffInd./715;

for j = 1: size(y,2)
cycle1 = ((firstFP(j):1:lastFP(j))-firstFP(j))/(lastFP(j)-firstFP(j));
cycle = 0:1/(2*winPts - 1):1;
figure(j)
subplot(2,1,1); plot(cycle1, Fy(firstFP(j):lastFP(j),2:end))
subplot(2,1,2); plot(y(:,j)); xlim([0 404])
end

end

