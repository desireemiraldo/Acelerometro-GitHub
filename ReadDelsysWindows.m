function Variable = ReadDelsysWindows(FilePath, ChannelType, Signal, winPts)

A =  importdata([FilePath,csv]);
Data = [A.textdata, num2cell(A.data)];
Label = Data(1,:);
Data(1,:) = [];

Variable = []; t = [];

linesNumber = 78*winPts;
BlocksNumber = (length(Data)-1)/linesNumber;

index = 1;
while index <= BlocksNumber 
    Trials(:,:,index) = [Label;Data(linesNumber*(index-1) + 1:linesNumber*index);
    index = index +1;
end


%%corrigir


for i = 1: height(Data)
   if  strcmp(Data.Var1(i),ChannelType) && strcmp(Data.Var2(i),Signal)
       Variable = [Variable; Data(i,5:end-1)];
       t = [t; str2num(cell2mat(Data.Var4(i)))];
   end
end
Variable = table2array(Variable);
% t = str2num(cell2mat(t.Var4));

Variable = [t,Variable];
end