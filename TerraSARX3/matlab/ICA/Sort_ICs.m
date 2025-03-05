function [ICs,SRs] = Sort_ICs(IC,SR)%(Z,W)
%%归一化处�?
%IC:时间序列；SR:是空间响�?
%sica的时候，[ICs,SRs] = Sort_ICs(W',Z')
%tica的时候，[ICs,SRs] = Sort_ICs(Z,W)
%Sorting the ICs using Ratio
%IC concludes the ICs in each row
%SR concludes the SRs in each column
X=SR*IC;
[p,q]=size(IC);
[m,n]=size(X);
%m stations
%p ICs

sumX=0;
for i=1:m
    x=X(i,:);
    sumx=sum(x.^2);   
    sumX=sumX+sumx;        
end


for i=1:p
    sr=SR(:,i);
    r=max(abs(sr));
    IC(i,:)=IC(i,:)*r;
    SR(:,i)=SR(:,i)/r;
    GPSICi=SR(:,i)*IC(i,:);
    sumICi=0;
    for j=1:m
        GPSICj=GPSICi(j,:);
        sumGPSICj=sum(GPSICj.^2); 
        sumICi=sumICi+sumGPSICj;
    end
    Ratio(i)=10*log10(sumX/sumICi);
   
end


[Ratio1 I]=sort(Ratio);

for i=1:p
    ICs(i,:)=IC(I(i),:);
    SRs(:,i)=SR(:,I(i));
end

end



