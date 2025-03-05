 clear;clc;
load K:\GAMMAUS1\TerraSARX3\output\ICA\Phase.mat; 
S=Phase;
[P,D] = pca(S', 10, 'eig');
%
for u=1:length(D)
   D1(u)=D(u)/sum(D);%D(38*1)
end
D2=zeros(1,length(D));
for u=1:length(D)
    D2(u)=sum(D1(1,1:u));
end




% %贡献率画单坐标轴折线�?
% x=1:1:10;
% y=[0.179721928941714,0.0835900187123608,0.0551681793485837,0.0499787741246863,0.0419684781326573,0.0372262942391094,0.0358279860263477,0.0331868427683266,0.0303749966731335,0.0298179078870903]; %a数据y�?
% plot(x,y,'-*b'); 
% axis([0,11,0,0.2]) 
% set(gca,'XTick',[0:1:11]) 
% set(gca,'YTick',[0:0.1:0.2])
% for i=1:10
%    text(x(i),y(i),num2str(y(i)));
% end




%
x=1:1:10;
y1=D1(1:10);
y2=D2(1:10);
figure;
[AX,h1,h2]=plotyy(x,y1,x,y2);
set(AX(1),'ycolor','r')
set(AX(2),'ycolor','b')
set(h1,'linestyle','-','marker','diamond','color','r')
set(h2,'linestyle','-','marker','pentagram','color','b')
for u=1:10
   text(AX(1),x(u)+0.1,y1(u)+0.002,num2str(y1(u),3));
   text(AX(2),x(u)+0.12,y2(u)-0.004,num2str(y2(u),3));
end
legend('Variance ','Cumulative');
title('Principal Component contribution rate');