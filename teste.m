Path = '.\Coletas\';
Name = 'Acc_170731_RNW_';
csv = '-Delsys 1.csv';
Trial = {'1' '2' '3' '4' '5'};
ShankL = 5;ShankR = 4;
Sensors = {'ShankR','ShankL'};
FsFP = 300;
Fs = 148.39;

f = [0 0 181 174 0 0 165 240];
l = [148 148 329 322 148 148 313 388];

instant = importdata('Instantes_gait.txt',',');
[HeelContact, ToeOff] = Instants(instant,File);
HeelContactInd = floor(HeelContact*Fs);

for j = 1: length(Trial)
    
    File = [Name,Trial{j}];
    FilePath = [Path,File];
    
    % Cortex data
    Forces = importdata([FilePath,'.forces']);
    Fy = (Forces.data(:,1) -1)/FsFP;
    for i =1 : length(Forces.colheaders)
        if strcmp(Forces.colheaders{i}(1:2), 'FY')
            Fy = [Fy, Forces.data(:,i)];
        end
    end
    
    % Plot Fy
    
    firstFP = zeros(length(Sensors),1);
    lastFP = zeros(length(Sensors),1);
    HeelContactInd = ceil(HeelContact*FsFP);
    
    for i = 1: length(Sensors)
        firstFP(i) = HeelContactInd(i,1);
        % lastFP(i) =  HeelContactInd(i,2);
        lastFP(i) =  HeelContactInd(i,1)+ceil(FsFP);
    end
    for i = 1 : length(Sensors)
        cycle1 = ((firstFP(i):1:lastFP(i))-firstFP(i))/(lastFP(i)-firstFP(i));
        figure(j+(i-1)*length(Trial));
        %subplot(2,1,1);
        plot(cycle1, Fy(firstFP(i):lastFP(i),2:end));
        title([Name,Trial{j},' (',Sensors{i},')']);
        ylim([min(Fy)-1 max(Fy)+1]);
    end
    instant = importdata('Instantes_gait.txt',',');
    [HeelContact, ToeOff] = Instants(instant,File);
    
    TO(:,j) = ToeOff(:,1);
end
for j = 1: length(Trial)
    for i = 1: length(Sensors)
        TOInd(i,j) = (floor(TO(i,j)*148.39)-f(i))/(l(i)-f(i))
    end
end