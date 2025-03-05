% 选择基线校正后结果的文件夹
% clear;clc;
% n=10;
% filepath1=uigetdir();%
% for i=1:n
%     load([filepath1 '\' 'y_baseline_correction' num2str(i) '.mat'])
% end
% Phase1=[y_baseline_correction1 y_baseline_correction2 y_baseline_correction3 y_baseline_correction4 y_baseline_correction5 y_baseline_correction6 y_baseline_correction7 y_baseline_correction8 y_baseline_correction9 y_baseline_correction10];
% save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase1','Phase1')
% 
% clear;clc;
% n=20;
% filepath1=uigetdir();%
% for i=11:n
%     load([filepath1 '\' 'y_baseline_correction' num2str(i) '.mat'])
% end
% Phase2=[y_baseline_correction11 y_baseline_correction12 y_baseline_correction13 y_baseline_correction14 y_baseline_correction15 y_baseline_correction16 y_baseline_correction17 y_baseline_correction18 y_baseline_correction19 y_baseline_correction20];
% save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase2','Phase2')
% 
% clear;clc;
% n=30;
% filepath1=uigetdir();%
% for i=21:n
%     load([filepath1 '\' 'y_baseline_correction' num2str(i) '.mat'])
% end
% Phase3=[y_baseline_correction21 y_baseline_correction22 y_baseline_correction23 y_baseline_correction24 y_baseline_correction25 y_baseline_correction26 y_baseline_correction27 y_baseline_correction28 y_baseline_correction29 y_baseline_correction30];
% save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase3','Phase3')
% 
% clear;clc;
% n=36;
% filepath1=uigetdir();%
% for i=31:n
%     load([filepath1 '\' 'y_baseline_correction' num2str(i) '.mat'])
% end
% Phase4=[y_baseline_correction31 y_baseline_correction32 y_baseline_correction33 y_baseline_correction34 y_baseline_correction35 y_baseline_correction36];
% save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase4','Phase4')
% clear;clc;
% load F:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase1.mat Phase1;
% load F:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase2.mat Phase2; 
% load F:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase3.mat Phase3;
% load F:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase4.mat Phase4;
% Phase=[Phase1 Phase2 Phase3 Phase4];
% save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase','Phase') 

clear;clc;
load F:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase.mat;
load F:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat; 
load F:\GAMMAUS1\TerraSARX3\output\dem\dem30_1.mat; 
load F:\GAMMAUS1\TerraSARX3\output\dem\finaldem1.mat;
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
save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\Def_LP_X','Def_LP_X');

clear;clc;
load F:\GAMMAUS1\TerraSARX3\output\BCDEM\Def_LP_X.mat; 
n=18;
V_LP=Def_LP_X(1:n,:); % deformation rate，unti:mm/year
SigmaH1=Def_LP_X(n+1,:);
V_LP=V_LP'; 
save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\SigmaH1','SigmaH1');
save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\V_LP','V_LP');

clear;clc;
% 形变求解,当内存不足时，可先清空文件，然后重新打开成果文件。
load F:\GAMMAUS1\TerraSARX3\output\BCDEM\SigmaH1.mat; 
load F:\GAMMAUS1\TerraSARX3\output\BCDEM\V_LP.mat; 
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
save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\S','S');
save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\SigmaH1','SigmaH1');

% for x=1:n
%     eval(['save(''F:\GAMMAUS1\TerraSARX3\output\BCDEM\V2\V2_',num2str(x),'.mat'',''V2_',num2str(x),''');']);
% end
% V2=zeros(PS_n,n);
% for i=1:n
%     V2(:,i)=eval(['V2_',num2str(i),'(:)']); 
% end
% save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\V2','V2');
% V2_mean=mean(V2,2);
% save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\V2','V2_mean');


% % 稀疏DEM反演
% HighlyDEM1=dem30_1+SigmaH1';
% % 生成稀疏矩阵
% H=zeros(7000,5000);x1=row_col_ct(:,1);y1=row_col_ct(:,2);
% H(sub2ind(size(H),x1,y1))=SigmaH1'; 
% % 整体DEM反演
% HighDEM1=finaldem1+H;

% save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\HighlyDEM1','HighlyDEM1');
% save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\HighDEM1','HighDEM1');
% save('F:\GAMMAUS1\TerraSARX3\output\BCDEM\H','H');
