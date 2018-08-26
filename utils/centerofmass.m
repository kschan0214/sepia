function coord=centerofmass(data);
data=abs(data);
dims=size(data);
for k=1:length(dims)
%     datatemp=permute(data,[k ]);
dimsvect=ones([1, length(dims)]);
dimsvect(k)=dims(k);
temp=bsxfun(@times,(data),reshape(1:dims(k),dimsvect));
coord(k)=sum(temp(:))./sum(data(:));
end;
