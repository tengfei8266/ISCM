% 对于小矩阵可采用循环，对于大矩阵干涉对分开拟合
% 干涉对1 对于小矩阵可采用循环，对于大矩阵

% clear;clc;
% load J:\GAMMAUS1\TerraSARX3\output\ICDEM\Phase.mat Phase;
% load J:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat row_col_ct; 
% for i=1:36
% chr=[num2str(i),'y.mat'];
% b=Phase(:,i);
% save(chr,'b');
% end
% PS_N=size(row_col_ct,1);
% a = randperm(PS_N); 
% x =(sort(a))';
% save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x','x');


% 干涉对1 
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\1y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y1=b;
y_baseline1=y_baseline;
y_baseline_correction1=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y1','y1');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline1','y_baseline1');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction1','y_baseline_correction1');
figure(1)
h=plot(x,y1,'bla',x,y_baseline_correction1,'r',x,y_baseline1,'c');

%干涉对2
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\2y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y2=b;
y_baseline2=y_baseline;
y_baseline_correction2=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y2','y2');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline2','y_baseline2');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction2','y_baseline_correction2');

figure(1)
h=plot(x,y2,'bla',x,y_baseline_correction2,'r',x,y_baseline2,'c');

%干涉对3
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\3y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y3=b;
y_baseline3=y_baseline;
y_baseline_correction3=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y3','y3');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline3','y_baseline3');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction3','y_baseline_correction3');

figure(1)
h=plot(x,y3,'bla',x,y_baseline_correction3,'r',x,y_baseline3,'c');


%干涉对4
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\4y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y4=b;
y_baseline4=y_baseline;
y_baseline_correction4=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y4','y4');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline4','y_baseline4');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction4','y_baseline_correction4');

figure(1)
h=plot(x,y4,'bla',x,y_baseline_correction4,'r',x,y_baseline4,'c');

%干涉对5
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\5y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y5=b;
y_baseline5=y_baseline;
y_baseline_correction5=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y5','y5');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline5','y_baseline5');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction5','y_baseline_correction5');

figure(1)
h=plot(x,y5,'bla',x,y_baseline_correction5,'r',x,y_baseline5,'c');

%干涉对6
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\6y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y6=b;
y_baseline6=y_baseline;
y_baseline_correction6=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y6','y6');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline6','y_baseline6');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction6','y_baseline_correction6');

figure(1)
h=plot(x,y6,'bla',x,y_baseline_correction6,'r',x,y_baseline6,'c');

%干涉对7
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\7y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y7=b;
y_baseline7=y_baseline;
y_baseline_correction7=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y7','y7');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline7','y_baseline7');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction7','y_baseline_correction7');

figure(1)
h=plot(x,y7,'bla',x,y_baseline_correction7,'r',x,y_baseline7,'c');

%干涉对8
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\8y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y8=b;
y_baseline8=y_baseline;
y_baseline_correction8=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y8','y8');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline8','y_baseline8');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction8','y_baseline_correction8');

figure(1)
h=plot(x,y8,'bla',x,y_baseline_correction8,'r',x,y_baseline8,'c');


%干涉对9
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\9y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y9=b;
y_baseline9=y_baseline;
y_baseline_correction9=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y9','y9');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline9','y_baseline9');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction9','y_baseline_correction9');

figure(1)
h=plot(x,y9,'bla',x,y_baseline_correction9,'r',x,y_baseline9,'c');

%干涉对10
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\10y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y10=b;
y_baseline10=y_baseline;
y_baseline_correction10=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y10','y10');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline10','y_baseline10');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction10','y_baseline_correction10');

figure(1)
h=plot(x,y10,'bla',x,y_baseline_correction10,'r',x,y_baseline10,'c');


clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\11y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y11=b;
y_baseline11=y_baseline;
y_baseline_correction11=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y11','y11');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline11','y_baseline11');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction11','y_baseline_correction11');
figure(1)
h=plot(x,y11,'bla',x,y_baseline_correction11,'r',x,y_baseline11,'c');

%干涉对12
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\12y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y12=b;
y_baseline12=y_baseline;
y_baseline_correction12=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y12','y12');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline12','y_baseline12');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction12','y_baseline_correction12');

figure(1)
h=plot(x,y12,'bla',x,y_baseline_correction12,'r',x,y_baseline12,'c');

%干涉对13
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\13y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y13=b;
y_baseline13=y_baseline;
y_baseline_correction13=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y13','y13');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline13','y_baseline13');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction13','y_baseline_correction13');

figure(1)
h=plot(x,y13,'bla',x,y_baseline_correction13,'r',x,y_baseline13,'c');


%干涉对14
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\14y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y14=b;
y_baseline14=y_baseline;
y_baseline_correction14=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y14','y14');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline14','y_baseline14');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction14','y_baseline_correction14');

figure(1)
h=plot(x,y14,'bla',x,y_baseline_correction14,'r',x,y_baseline14,'c');

%干涉对15
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\15y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y15=b;
y_baseline15=y_baseline;
y_baseline_correction15=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y15','y15');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline15','y_baseline15');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction15','y_baseline_correction15');

figure(1)
h=plot(x,y15,'bla',x,y_baseline_correction15,'r',x,y_baseline15,'c');

%干涉对16
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\16y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y16=b;
y_baseline16=y_baseline;
y_baseline_correction16=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y16','y16');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline16','y_baseline16');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction16','y_baseline_correction16');

figure(1)
h=plot(x,y16,'bla',x,y_baseline_correction16,'r',x,y_baseline16,'c');

%干涉对17
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\17y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y17=b;
y_baseline17=y_baseline;
y_baseline_correction17=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y17','y17');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline17','y_baseline17');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction17','y_baseline_correction17');

figure(1)
h=plot(x,y17,'bla',x,y_baseline_correction17,'r',x,y_baseline17,'c');

%干涉对18
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\18y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y18=b;
y_baseline18=y_baseline;
y_baseline_correction18=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y18','y18');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline18','y_baseline18');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction18','y_baseline_correction18');

figure(1)
h=plot(x,y18,'bla',x,y_baseline_correction18,'r',x,y_baseline18,'c');


%干涉对19
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\19y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y19=b;
y_baseline19=y_baseline;
y_baseline_correction19=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y19','y19');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline19','y_baseline19');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction19','y_baseline_correction19');

figure(1)
h=plot(x,y19,'bla',x,y_baseline_correction19,'r',x,y_baseline19,'c');

%干涉对20
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\20y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y20=b;
y_baseline20=y_baseline;
y_baseline_correction20=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y20','y20');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline20','y_baseline20');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction20','y_baseline_correction20');

figure(1)
h=plot(x,y20,'bla',x,y_baseline_correction20,'r',x,y_baseline20,'c');


% 干涉对21 
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\21y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y21=b;
y_baseline21=y_baseline;
y_baseline_correction21=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y21','y21');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline21','y_baseline21');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction21','y_baseline_correction21');
figure(1)
h=plot(x,y21,'bla',x,y_baseline_correction21,'r',x,y_baseline21,'c');

%干涉对22
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\22y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y22=b;
y_baseline22=y_baseline;
y_baseline_correction22=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y22','y22');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline22','y_baseline22');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction22','y_baseline_correction22');

figure(1)
h=plot(x,y22,'bla',x,y_baseline_correction22,'r',x,y_baseline22,'c');

%干涉对23
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\23y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y23=b;
y_baseline23=y_baseline;
y_baseline_correction23=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y23','y23');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline23','y_baseline23');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction23','y_baseline_correction23');

figure(1)
h=plot(x,y23,'bla',x,y_baseline_correction23,'r',x,y_baseline23,'c');


%干涉对24
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\24y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y24=b;
y_baseline24=y_baseline;
y_baseline_correction24=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y24','y24');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline24','y_baseline24');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction24','y_baseline_correction24');

figure(1)
h=plot(x,y24,'bla',x,y_baseline_correction24,'r',x,y_baseline24,'c');

%干涉对25
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\25y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y25=b;
y_baseline25=y_baseline;
y_baseline_correction25=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y25','y25');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline25','y_baseline25');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction25','y_baseline_correction25');

figure(1)
h=plot(x,y25,'bla',x,y_baseline_correction25,'r',x,y_baseline25,'c');

%干涉对26
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\26y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y26=b;
y_baseline26=y_baseline;
y_baseline_correction26=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y26','y26');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline26','y_baseline26');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction26','y_baseline_correction26');

figure(1)
h=plot(x,y26,'bla',x,y_baseline_correction26,'r',x,y_baseline26,'c');

%干涉对27
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\27y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y27=b;
y_baseline27=y_baseline;
y_baseline_correction27=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y27','y27');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline27','y_baseline27');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction27','y_baseline_correction27');

figure(1)
h=plot(x,y27,'bla',x,y_baseline_correction27,'r',x,y_baseline27,'c');

%干涉对28
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\28y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y28=b;
y_baseline28=y_baseline;
y_baseline_correction28=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y28','y28');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline28','y_baseline28');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction28','y_baseline_correction28');

figure(1)
h=plot(x,y28,'bla',x,y_baseline_correction28,'r',x,y_baseline28,'c');


%干涉对29
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\29y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y29=b;
y_baseline29=y_baseline;
y_baseline_correction29=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y29','y29');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline29','y_baseline29');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction29','y_baseline_correction29');

figure(1)
h=plot(x,y29,'bla',x,y_baseline_correction29,'r',x,y_baseline29,'c');

%干涉对30
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\30y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y30=b;
y_baseline30=y_baseline;
y_baseline_correction30=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y30','y30');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline30','y_baseline30');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction30','y_baseline_correction30');

figure(1)
h=plot(x,y30,'bla',x,y_baseline_correction30,'r',x,y_baseline30,'c');

% 干涉对31 
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\31y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y31=b;
y_baseline31=y_baseline;
y_baseline_correction31=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y31','y31');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline31','y_baseline31');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction31','y_baseline_correction31');
figure(1)
h=plot(x,y31,'bla',x,y_baseline_correction31,'r',x,y_baseline31,'c');

%干涉对32
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\32y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y32=b;
y_baseline32=y_baseline;
y_baseline_correction32=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y32','y32');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline32','y_baseline32');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction32','y_baseline_correction32');

figure(1)
h=plot(x,y32,'bla',x,y_baseline_correction32,'r',x,y_baseline32,'c');

%干涉对33
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\33y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y33=b;
y_baseline33=y_baseline;
y_baseline_correction33=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y33','y33');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline33','y_baseline33');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction33','y_baseline_correction33');

figure(1)
h=plot(x,y33,'bla',x,y_baseline_correction33,'r',x,y_baseline33,'c');


%干涉对34
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\34y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y34=b;
y_baseline34=y_baseline;
y_baseline_correction34=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y34','y34');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline34','y_baseline34');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction34','y_baseline_correction34');

figure(1)
h=plot(x,y34,'bla',x,y_baseline_correction34,'r',x,y_baseline34,'c');

%干涉对35
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\35y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y35=b;
y_baseline35=y_baseline;
y_baseline_correction35=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y35','y35');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline35','y_baseline35');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction35','y_baseline_correction35');

figure(1)
h=plot(x,y35,'bla',x,y_baseline_correction35,'r',x,y_baseline35,'c');

%干涉对36
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\x.mat x; 
load J:\GAMMAUS1\TerraSARX3\matlab\BCDEM\36y.mat b;

PS_N=size(x,1);
y = sgolayfilt(b,6,11);

n=2;%多项式阶数

[p0,s0,mu0]=polyfit(x,y,n);%多项式拟合
y_fit0=polyval(p0,x,[],mu0);%计算拟合值
r0=y-y_fit0;
dev0=sqrt(sum((r0-mean(r0)).^2)/length(r0));%计算残差
y_remove0=y(find(y<=y_fit0));%峰值消除
x_remove0=x(find(y<=y_fit0));%峰值消除
i=1;
judge=1;
while(judge)
[p1,s1,mu1]=polyfit(x_remove0,y_remove0,n);%多项式拟合
y_fit1=polyval(p1,x_remove0,[],mu1);%计算拟合值
r1=y_remove0-y_fit1;
dev(i)=sqrt(sum((r1-mean(r1)).^2)/length(r1));%计算残差
if i==1
    judge=abs(dev(i)-dev0)/dev(i)>0.05;
else
    judge=abs((dev(i)-dev(i-1))/dev(i))>0.05;%残差判断条件
end
index=find(y_remove0<=y_fit1);
y_remove0(index)=y_fit1(index);%光谱重建,大于拟合值的数据用拟合值代替，小于拟合值的数组采用原始数据
i=i+1;
end
y_baseline=polyval(p1,x,[],mu1);%基线
find(y>0);

%基线去除条件，根据相位信号的正负值决定
 y_baseline_correction=zeros(PS_N,1);
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline>0))=y(find(y>y_baseline & y>0 & y_baseline>0))-y_baseline(find(y>y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y<y_baseline & y>0 & y_baseline>0))=-y(find(y<y_baseline & y>0 & y_baseline>0))+y_baseline(find(y<y_baseline & y>0 & y_baseline>0));
 
 y_baseline_correction(find(y>y_baseline & y>0 & y_baseline<0))=y(find(y>y_baseline & y>0 & y_baseline<0))+y_baseline(find(y>y_baseline & y>0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline<0))=y(find(y<y_baseline & y<0 & y_baseline<0))-y_baseline(find(y<y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y>y_baseline & y<0 & y_baseline<0))=-y(find(y>y_baseline & y<0 & y_baseline<0))+y_baseline(find(y>y_baseline & y<0 & y_baseline<0));
 
 y_baseline_correction(find(y<y_baseline & y<0 & y_baseline>0))=y(find(y<y_baseline & y<0 & y_baseline>0))+y_baseline(find(y<y_baseline & y<0 & y_baseline>0));
hold on
hold off
y36=b;
y_baseline36=y_baseline;
y_baseline_correction36=y_baseline_correction;
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y36','y36');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline36','y_baseline36');
save('J:\GAMMAUS1\TerraSARX3\output\BCDEM\y_baseline_correction36','y_baseline_correction36');

figure(1)
h=plot(x,y36,'bla',x,y_baseline_correction36,'r',x,y_baseline36,'c');

