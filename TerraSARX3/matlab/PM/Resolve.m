% 计算沉降速率和高程改正值

% 基于沉降速率的DInSAR模型，（笔记3，或者尹宏杰的论文）
% 对于第j行，位于主辅影像获取时间之间的列B(j,k)=t（k+1)-t(k)，否则B=0.

clear;clc;
load K:\GAMMAUS1\TerraSARX3\output\highcp\Num_mask_lungui.mat Num_mask_lungui;
load K:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat row_col_ct; 

PS_n=size(row_col_ct,1); M=size(Num_mask_lungui,2); n=3; % n为影像个数减一，共获取n+1景SAR影像
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
T=[22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
TT=[22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264].^2;
TTT=[22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264].^3;
T1=[T(:,1) T(:,3)-T(:,1)  T(:,4)-T(:,1) T(:,4)-T(:,2) T(:,6)-T(:,2) T(:,5)-T(:,3) T(:,6)-T(:,3) T(:,7)-T(:,3) T(:,6)-T(:,4) T(:,9)-T(:,4) T(:,8)-T(:,5) T(:,9)-T(:,5) T(:,7)-T(:,6) T(:,9)-T(:,6) T(:,8)-T(:,7) T(:,10)-T(:,7) T(:,10)-T(:,8) T(:,11)-T(:,8) T(:,11)-T(:,9) T(:,15)-T(:,9) T(:,11)-T(:,10) T(:,15)-T(:,10) T(:,13)-T(:,11) T(:,15)-T(:,11) T(:,17)-T(:,11) T(:,13)-T(:,12) T(:,18)-T(:,12) T(:,14)-T(:,13) T(:,18)-T(:,13) T(:,15)-T(:,14) T(:,17)-T(:,11) T(:,17)-T(:,15) T(:,18)-T(:,15) T(:,17)-T(:,16) T(:,18)-T(:,16) T(:,18)-T(:,17)];
T2=T1.^2;
T3=T1.^3;

B1=((-4*pi)/Lambda).*T1;
B2=(((-4*pi)/Lambda).*T2)/2;
B3=(((-4*pi)/Lambda).*T3)/6;
% TT=[780 768 756 744 720 708 696 684 672 660 648 636 624 612 588 576 564 552 528 516 504 492 480 468 432 384 360 336 312 288 264 240 216 192 120 96 72 48 24 0];

B13=[B1',B2',B3',B_perp];
% Low frequency linear deformation
% for multi-v model
% 将相位沿LOS方向转化为平均相位速度V=[ V1=(phi1-phi0)/(t1-t0),V2=(phi2-ph1)/(t2-t1),...V(N-1)=(phi(n-1)-phi(N-2))/(t(N-1)-t(N-2))
Def_LP_X=zeros(n+1,PS_n); 

for i=1:PS_n
    L3=Num_mask_lungui(i,:)'; % High coherence's phase
    Def_LP_X(:,i)=SVD(B13,L3);     % BX_p=L；X_p=INV(B)L；在这里L相当于方程的常数矩阵 
end

V1=Def_LP_X(1,:); % deformation rate，unti:mm/year
V2=Def_LP_X(2,:); % deformation rate，unti:mm/year
V3=Def_LP_X(3,:); % deformation rate，unti:mm/year
S=1000.*(V1'.*T+(V2'.*TT)/2+(V3'.*TTT)/6)/cosd(theta);
SigmaH1=Def_LP_X(n+1,:);

save('K:\GAMMAUS1\TerraSARX3\output\CTMDEM\S','S');
save('K:\GAMMAUS1\TerraSARX3\output\CTMDEM\SigmaH1','SigmaH1');

