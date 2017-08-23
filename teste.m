Name = 'Acc_170731_RNW_';
csv = '-Delsys 1.csv';
Trial = {'1' '2' '3' '4' '5'};
ShankL = 5;ShankR = 4;
Sensors = {'ShankR','ShankL'};

f = [0 0 181 174 0 0 165 240];
l = [148 148 329 322 148 148 313 388];

for j = 1: length(Trial)
    
    File = [Name,Trial{j}];
    instant = importdata('Instantes_gait.txt',',');
    [HeelContact, ToeOff] = Instants(instant,File);
    
    TO(:,j) = ToeOff(:,1);
end
for j = 1: length(Trial)
    for i = 1: length(Sensors)
        TOInd(i,j) = (floor(TO(i,j)*148.39)-f(i))/(l(i)-f(i))
    end
end