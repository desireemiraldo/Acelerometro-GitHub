function [newInstant] = ReshapeInstants(Data, instant,Name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

A = strcmp(instant.textdata(:,1),Name);
ind = find(A==1)-1; %first row is label

newInstant = NaN;

for k = 1 : size(Data,3)
    
    first = str2double(Data(2,4,k));
    last = str2double(Data(end,4,k));
    
    i = find(instant.data(ind(1),:)>=first & instant.data(ind(1),:)<=last); %edit when using both legs
    
    New = instant.data(ind(1),i);
       
    ColNumber = size(newInstant,2);
    
    if length(New)> ColNumber
        newInstant(:,ColNumber + 1 :length(New)) = NaN;
    end
    if length(New)< ColNumber
        newInstant(k,:) = NaN;
    end
    newInstant(k,1:length(New))  = New;

end



end

