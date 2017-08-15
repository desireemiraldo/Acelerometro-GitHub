function [HeelContact] = Instants(instant,File)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
HeelContact = [];
k=1;
for i=1:2:length(instant.textdata)
   if strcmp(instant.textdata(i,1),File)
       for j = 1 : size(instant.data,2)
           if instant.data (i,j)> 0
               HeelContact(1,k) = instant.data(i,j);
               HeelContact(2,k) = instant.data(i+1,j);
               k=k+1;
           end
       end
   end
end
               