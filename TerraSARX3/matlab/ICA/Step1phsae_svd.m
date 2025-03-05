%将干涉对上的相位用SVD分解得到时序上的相位
clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\highcp\Num_mask_lungui.mat Num_mask_lungui; 
load J:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat row_col_ct; 
a=isnan(Num_mask_lungui);
Num_mask_lungui(a)=0;
PS_n=size(row_col_ct,1); M=size(Num_mask_lungui,2); n=18; % n为影像个数减一，共获取n+1景SAR影像
HCPPhase=Num_mask_lungui;

Lambda=0.032 ; % center_range_slc  地理编码的文件中可以找到
theta=39.2655;  %  incident angle
R=645683.1875; % 轨道半径 center_range_slc
t1=cos(theta*pi./180); % cosd计算度数，而cos计算弧度
t2=sind(theta);
q=-4*pi./Lambda;

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


% Low frequency linear deformation
% for multi-v model
% 将相位沿LOS方向转化为平均相位速度V=[ V1=(phi1-phi0)/(t1-t0),V2=(phi2-ph1)/(t2-t1),...V(N-1)=(phi(n-1)-phi(N-2))/(t(N-1)-t(N-2))
Def_LP_X=zeros(n,PS_n); 
B13=B1; % 正负号问题可以参考廖明生的书第八节

for i=1:PS_n
    L3=Num_mask_lungui(i,:)'; % High coherence's phase
    Phase_V(:,i)=SVD(B13,L3);     % BX_p=L；X_p=INV(B)L；在这里L相当于方程的常数矩阵 
end
   
save J:\GAMMAUS1\TerraSARX3\output\ICA\Phase_V.mat Phase_V; 
save J:\GAMMAUS1\TerraSARX3\output\ICA\row_col_ct.mat row_col_ct 

clear;clc;
load J:\GAMMAUS1\TerraSARX3\output\ICA\Phase_V.mat Phase_V;
Lambda=0.032 ; % center_range_slc  地理编码的文件中可以找到
theta=39.2655;  %  incident angle
TT=[0 22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
T=[0 22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264]./365;

n=18;
DT=zeros(1,n);
for i=1:n
    DT(i)=(TT(i+1)-TT(i))./365;     % ./365; % for multi-v model, the 15th differenced T
end

Phase_V=Phase_V';
Phase(:,1)=DT(1)*Phase_V(:,1); 

for i=2:n
    Phase(:,i)=Phase(:,i-1)+DT(i)*Phase_V(:,i);
end
S=1000.*Phase*(Lambda/(4*pi))./(cosd(theta)); % vetical deformation-
save J:\GAMMAUS1\TerraSARX3\output\ICA\Phase.mat Phase; 
