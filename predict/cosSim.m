function [ result ] = cosSim( data )
%COSSIM Summary of this function goes here
%   Detailed explanation goes here
   
rows=size(data,1);
result=zeros(rows,rows);

    for i=1:rows
        
        for j=1:rows
            
            result(i,j)=dot(data(i,:),data(j,:))/(norm(data(i,:))*norm(data(j,:)));
            
        end
        
        
    end

end

