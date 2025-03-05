% 计算沉降速率和高程改正值

% 基于沉降速率的DInSAR模型，（笔记3，或者尹宏杰的论文）
% 对于第j行，位于主辅影像获取时间之间的列B(j,k)=t（k+1)-t(k)，否则B=0.

clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\ICA\SICA\SS.mat;
load J:\GAMMAUS1\TerraSARX3\output\ICA\SICA\S1.mat;
load J:\GAMMAUS1\TerraSARX3\output\ICA\SICA\S2.mat;
load J:\GAMMAUS1\TerraSARX3\output\ICA\SICA\S3.mat;
load J:\GAMMAUS1\TerraSARX3\output\ICA\SICA\S4.mat;

Phase=S1+S2+S3+S4;
save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\Phase','Phase');

clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\ICDEM\Phase.mat;
Phase=[Phase(:,2) Phase(:,3)-Phase(:,1) Phase(:,4)-Phase(:,1) Phase(:,4)-Phase(:,2) Phase(:,6)-Phase(:,2) Phase(:,4)-Phase(:,3) Phase(:,6)-Phase(:,3) Phase(:,7)-Phase(:,3) Phase(:,6)-Phase(:,4) Phase(:,9)-Phase(:,4) Phase(:,8)-Phase(:,5) Phase(:,9)-Phase(:,5) Phase(:,7)-Phase(:,6) Phase(:,9)-Phase(:,6) Phase(:,8)-Phase(:,7) Phase(:,10)-Phase(:,7) Phase(:,10)-Phase(:,8) Phase(:,11)-Phase(:,8) Phase(:,11)-Phase(:,9) Phase(:,15)-Phase(:,9) Phase(:,11)-Phase(:,10) Phase(:,15)-Phase(:,10) Phase(:,13)-Phase(:,11) Phase(:,15)-Phase(:,11) Phase(:,17)-Phase(:,11) Phase(:,13)-Phase(:,12) Phase(:,18)-Phase(:,12) Phase(:,14)-Phase(:,13) Phase(:,18)-Phase(:,13) Phase(:,15)-Phase(:,14) Phase(:,17)-Phase(:,14) Phase(:,17)-Phase(:,15) Phase(:,18)-Phase(:,15) Phase(:,17)-Phase(:,16) Phase(:,18)-Phase(:,16) Phase(:,18)-Phase(:,17)];
save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\Phase','Phase');

clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\ICDEM\Phase.mat;
load J:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat row_col_ct; 
load J:\GAMMAUS1\TerraSARX3\output\dem\dem30_1.mat; 
load J:\GAMMAUS1\TerraSARX3\output\dem\finaldem1.mat;
PS_n=size(row_col_ct,1); M=size(Phase,2); n=18; % n为影像个数减一，共获取n+1景SAR影像
% 雷达参数
Lambda=0.032 ; % center_range_slc  地理编码的文件中可以找到
theta=39.2655;  %  incident angle
R=645683.1875; % 轨道半径 center_range_slc
t1=cos(theta*pi./180); % cosd计算度数，而cos计算弧度
t2=sind(theta);
q=-4*pi./Lambda;


%DEM高程改正项
B_p=[-157.2305 35.1051 69.3624 -137.2989 44.0192...
        34.2578 215.5732 260.3242 181.3243 186.1814...
        -6.4580 -72.4147 44.7433 4.8601 26.0727...
        -38.9418 -65.0142 -111.0810  -45.1226  90.8878... 
        -46.0639   89.9407  -194.8841 136.0074 -24.5105...
         8.4395 126.2530 245.0317 117.8157  85.8539...
        -74.6641 -160.5214 -213.0784 155.1250 102.5684 -52.5601]'; % 垂直基线的长度
 q1=Lambda.*R.*sind(theta);
 B_perp=B_p*((-4*pi)/q1); % 高程改正项系数
 

%多速率形变项
TT=[0 22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
T=[0 22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264]./365;
% TT=[780 768 756 744 720 708 696 684 672 660 648 636 624 612 588 576 564 552 528 516 504 492 480 468 432 384 360 336 312 288 264 240 216 192 120 96 72 48 24 0];

% B矩阵的建立可以参考"穆家雷"的《基于SBAAS_InSAR的深部硬石膏开采引起地表形变规律研究》
DT=zeros(1,n);
for i=1:n
    DT(i)=(TT(i+1)-TT(i))./365;     % ./365; % for multi-v model, the 15th differenced T
end
    B1=zeros(M,n);
%0-1 
B1(1,1)=[DT(1)];
%1-3 1-4
B1(2,2:3)=[DT(2) DT(3)];
B1(3,2:4)=[DT(2) DT(3) DT(4)];
%2-4;2-6
B1(4,3:4)=[DT(3) DT(4)];
B1(5,3:6)=[DT(3) DT(4) DT(5) DT(6)];
%3-5,3-6,3-7
B1(6,4:5)=[DT(4) DT(5)];
B1(7,4:6)=[DT(4) DT(5) DT(6)];
B1(8,4:7)=[DT(4) DT(5) DT(6) DT(7)];
%4-6，4-9
B1(9,5:6)=[DT(5) DT(6)];
B1(10,5:9)=[DT(5) DT(6) DT(7) DT(8) DT(9)];
%5-8，5-9
B1(11,6:8)=[DT(6) DT(7) DT(8)];
B1(12,6:9)=[DT(6) DT(7) DT(8) DT(9)];
%6-7，6-9
B1(13,7)=[DT(7)];
B1(14,7:9)=[DT(7) DT(8) DT(9)];
%7-8，7-10
B1(15,8)=[DT(8)];
B1(16,8:10)=[DT(8) DT(9) DT(10)];
%8-10，8-11
B1(17,9:10)=[DT(9) DT(10)];
B1(18,9:11)=[DT(9) DT(10) DT(11)];
%9-11，9-15
B1(19,10:11)=[DT(10) DT(11)];
B1(20,10:15)=[DT(10) DT(11) DT(12) DT(13) DT(14)  DT(15)];
%10-11，10-15
B1(21,11)=[DT(11)];
B1(22,11:15)=[DT(11) DT(12) DT(13) DT(14)  DT(15)];
%11-13，11-15，11-17
B1(23,12:13)=[DT(12) DT(13)];
B1(24,12:15)=[DT(12) DT(13) DT(14) DT(15)];
B1(25,12:17)=[DT(12) DT(13) DT(14) DT(15) DT(16) DT(17)];
%12-13，12-18
B1(26,13)=[DT(13)];
B1(27,13:18)=[DT(13) DT(14) DT(15) DT(16) DT(17) DT(18)];
%13-14，13-18
B1(28,14)=[DT(14)];
B1(29,14:18)=[DT(14) DT(15) DT(16) DT(17) DT(18)];
%14-15，14-17
B1(30,15)=[DT(15)];
B1(31,15:17)=[DT(15) DT(16) DT(17)];
%15-17，15-18
B1(32,16:17)=[DT(16) DT(17)];
B1(33,16:18)=[DT(16) DT(17) DT(18)];
%16-17，16-18
B1(34,17)=[DT(17)];
B1(35,17:18)=[DT(17) DT(18)];
%17-18
B1(36,18)=[DT(18)];

B=B1*((-4*pi)/Lambda); % 正负号问题可以参考廖明生的书第八节
B13=[B,B_perp];
% Low frequency linear deformation
% for multi-v model
% 将相位沿LOS方向转化为平均相位速度V=[ V1=(phi1-phi0)/(t1-t0),V2=(phi2-ph1)/(t2-t1),...V(N-1)=(phi(n-1)-phi(N-2))/(t(N-1)-t(N-2))
Def_LP_X=zeros(n+1,PS_n); 

for i=1:PS_n
    L3=Phase(i,:)'; % High coherence's phase
    Def_LP_X(:,i)=SVD(B13,L3);     % BX_p=L；X_p=INV(B)L；在这里L相当于方程的常数矩阵 
end
save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\Def_LP_X','Def_LP_X');

clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\ICDEM\Def_LP_X.mat; 
n=18;
V_LP=Def_LP_X(1:n,:); % deformation rate，unti:mm/year
SigmaH1=Def_LP_X(n+1,:);
V_LP=V_LP'; 
save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\SigmaH1','SigmaH1');
save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\V_LP','V_LP');

clear;clc;
% 形变求解,当内存不足时，可先清空文件，然后重新打开成果文件。
load J:\GAMMAUS1\TerraSARX3\output\ICDEM\SigmaH1.mat; 
load J:\GAMMAUS1\TerraSARX3\output\ICDEM\V_LP.mat; 
TT=[0 22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
n=18;theta=39.2655;  %  incident angle
DT=zeros(1,n);
for i=1:n
    DT(i)=(TT(i+1)-TT(i))./365;     % ./365; % for multi-v model, the 15th differenced T
end
S_Los(:,1)=DT(1)*V_LP(:,1);  % LOS deformation value, v*t
for i=2:n
    S_Los(:,i)=S_Los(:,i-1)+DT(i)*V_LP(:,i);
end
S=1000.*S_Los./(cosd(theta)); % vetical deformation
for x=1:n
    eval(['V2_',num2str(x),'=V_LP(:,',num2str(x),');'])
end
save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\S','S');
save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\SigmaH1','SigmaH1');

% for x=1:n
%     eval(['save(''J:\GAMMAUS1\TerraSARX3\output\ICDEM\V2\V2_',num2str(x),'.mat'',''V2_',num2str(x),''');']);
% end
% V2=zeros(PS_n,n);
% for i=1:n
%     V2(:,i)=eval(['V2_',num2str(i),'(:)']); 
% end
% save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\V2','V2');
% V2_mean=mean(V2,2);
% save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\V2','V2_mean');


% % 稀疏DEM反演
% HighlyDEM1=dem30_1+SigmaH1';
% % 生成稀疏矩阵
% H=zeros(7000,5000);x1=row_col_ct(:,1);y1=row_col_ct(:,2);
% H(sub2ind(size(H),x1,y1))=SigmaH1'; 
% % 整体DEM反演
% HighDEM1=finaldem1+H;

% save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\HighlyDEM1','HighlyDEM1');
% save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\HighDEM1','HighDEM1');
% save('J:\GAMMAUS1\TerraSARX3\output\ICDEM\H','H');
