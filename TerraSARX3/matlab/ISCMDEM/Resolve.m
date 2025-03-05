clear;clc;
load K:\GAMMAUS1\TerraSARX3\output\BCDEM\Phase\Phase.mat Phase;
load K:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat row_col_ct; 


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
        -46.0639   89.9407  -194.8841 136.0074 -24.5105...SigmaH-4
         8.4395 126.2530 245.0317 117.8157  85.8539...
        -74.6641 -160.5214 -213.0784 155.1250 102.5684 -52.5601]'; % 垂直基线的长度
Lambda=0.028;
q1=Lambda.*R.*sind(theta);
B_perp=B_p*((-4*pi)/q1); % 高程改正项系数
 

%多速率形变项
TT=[0 22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
T=[22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
% TT=[780 768 756 744 720 708 696 684 672 660 648 636 624 612 588 576 564 552 528 516 504 492 480 468 432 384 360 336 312 288 264 240 216 192 120 96 72 48 24 0];

% B矩阵的建立可以参考"穆家雷"的《基于SBAAS_InSAR的深部硬石膏开采引起地表形变规律研究》
DT=zeros(1,n);
for i=1:n
    DT(i)=(TT(i+1)-TT(i))./365;     % ./365; % for multi-v model, the 15th differenced T
end
    B1=zeros(M,n);
%0-2 
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

% 周期形变的系数
TT1=[T(:,1) T(:,3)-T(:,1)  T(:,4)-T(:,1) T(:,4)-T(:,2) T(:,6)-T(:,2) T(:,5)-T(:,3) T(:,6)-T(:,3) T(:,7)-T(:,3) T(:,6)-T(:,4) T(:,9)-T(:,4) T(:,8)-T(:,5) T(:,9)-T(:,5) T(:,7)-T(:,6) T(:,9)-T(:,6) T(:,8)-T(:,7) T(:,10)-T(:,7) T(:,10)-T(:,8) T(:,11)-T(:,8) T(:,11)-T(:,9) T(:,15)-T(:,9) T(:,11)-T(:,10) T(:,15)-T(:,10) T(:,13)-T(:,11) T(:,15)-T(:,11) T(:,17)-T(:,11) T(:,13)-T(:,12) T(:,18)-T(:,12) T(:,14)-T(:,13) T(:,18)-T(:,13) T(:,15)-T(:,14) T(:,17)-T(:,11) T(:,17)-T(:,15) T(:,18)-T(:,15) T(:,17)-T(:,16) T(:,18)-T(:,16) T(:,18)-T(:,17)];
T2=[T(:,1) T(:,3)  T(:,4) T(:,4) T(:,6) T(:,5) T(:,6) T(:,7) T(:,6) T(:,9) T(:,8) T(:,9) T(:,7) T(:,9) T(:,8)  T(:,10) T(:,10) T(:,11) T(:,11) T(:,15) T(:,11) T(:,15) T(:,13) T(:,15) T(:,17) T(:,13) T(:,18) T(:,14) T(:,18) T(:,15) T(:,17) T(:,17) T(:,18) T(:,17) T(:,18) T(:,18)];
T1=[0      T(:,1)  T(:,1) T(:,2) T(:,2) T(:,3) T(:,3) T(:,3) T(:,4) T(:,4) T(:,5) T(:,5) T(:,6) T(:,6) T(:,7)  T(:,7)  T(:,8)  T(:,8)  T(:,9)  T(:,9)  T(:,10) T(:,10) T(:,11) T(:,11) T(:,11) T(:,12) T(:,12) T(:,13) T(:,13) T(:,14) T(:,14) T(:,15) T(:,15) T(:,16) T(:,16) T(:,17)];
% BB1=((-4*pi)/Lambda).*TT1;
B11=(sind(T2.*2*pi/365)-sind(T1.*2*pi/365))*((-4*pi)/Lambda);
B22=(cosd(T2.*2*pi/365)-cosd(T1.*2*pi/365))*((-4*pi)/Lambda);

% 环境形变的系数
Tem=[23.8375  23.625	21.375	19.125	18.125  17.125	15 	    12.875	13.437	14	    13.75	13.5	13.70	13.91	14.125	14.33	14.54  14.75]/30;
Pre=[0.8      0.92	   	5.7	     5.7    6.25	6.25	72	    72	    138.5	138.5	57.5	 57.5	21.5	21.5	21.5	34      34     34]/30;
Te=[Tem(:,1)-23.56 Tem(:,3)-Tem(:,1)  Tem(:,4)-Tem(:,1) Tem(:,4)-Tem(:,2) Tem(:,6)-Tem(:,2) Tem(:,5)-Tem(:,3) Tem(:,6)-Tem(:,3) Tem(:,7)-Tem(:,3) Tem(:,6)-Tem(:,4) Tem(:,9)-Tem(:,4) Tem(:,8)-Tem(:,5) Tem(:,9)-Tem(:,5) Tem(:,7)-Tem(:,6) Tem(:,9)-Tem(:,6) Tem(:,8)-Tem(:,7) Tem(:,10)-Tem(:,7) Tem(:,10)-Tem(:,8) Tem(:,11)-Tem(:,8) Tem(:,11)-Tem(:,9) Tem(:,15)-Tem(:,9) Tem(:,11)-Tem(:,10) Tem(:,15)-Tem(:,10) Tem(:,13)-Tem(:,11) Tem(:,15)-Tem(:,11) Tem(:,17)-Tem(:,11) Tem(:,13)-Tem(:,12) Tem(:,18)-Tem(:,12) Tem(:,14)-Tem(:,13) Tem(:,18)-Tem(:,13) Tem(:,15)-Tem(:,14) Tem(:,17)-Tem(:,11) Tem(:,17)-Tem(:,15) Tem(:,18)-Tem(:,15) Tem(:,17)-Tem(:,16) Tem(:,18)-Tem(:,16) Tem(:,18)-Tem(:,17)];
Pr=[Pre(:,1)-0.800 Pre(:,3)-Pre(:,1)  Pre(:,4)-Pre(:,1) Pre(:,4)-Pre(:,2) Pre(:,6)-Pre(:,2) Pre(:,5)-Pre(:,3) Pre(:,6)-Pre(:,3) Pre(:,7)-Pre(:,3) Pre(:,6)-Pre(:,4) Pre(:,9)-Pre(:,4) Pre(:,8)-Pre(:,5) Pre(:,9)-Pre(:,5) Pre(:,7)-Pre(:,6) Pre(:,9)-Pre(:,6) Pre(:,8)-Pre(:,7) Pre(:,10)-Pre(:,7) Pre(:,10)-Pre(:,8) Pre(:,11)-Pre(:,8) Pre(:,11)-Pre(:,9) Pre(:,15)-Pre(:,9) Pre(:,11)-Pre(:,10) Pre(:,15)-Pre(:,10) Pre(:,13)-Pre(:,11) Pre(:,15)-Pre(:,11) Pre(:,17)-Pre(:,11) Pre(:,13)-Pre(:,12) Pre(:,18)-Pre(:,12) Pre(:,14)-Pre(:,13) Pre(:,18)-Pre(:,13) Pre(:,15)-Pre(:,14) Pre(:,17)-Pre(:,11) Pre(:,17)-Pre(:,15) Pre(:,18)-Pre(:,15) Pre(:,17)-Pre(:,16) Pre(:,18)-Pre(:,16) Pre(:,18)-Pre(:,17)];

B33=Te*((-4*pi)/Lambda);
B44=Pr*((-4*pi)/Lambda);

B13=[B,B11',B22',B33',B44',B_perp];
% Low frequency linear deformation
% for multi-v model
% 将相位沿LOS方向转化为平均相位速度V=[ V1=(phi1-phi0)/(t1-t0),V2=(phi2-ph1)/(t2-t1),...V(N-1)=(phi(n-1)-phi(N-2))/(t(N-1)-t(N-2))
Def_LP_X=zeros(n+5,PS_n); 

for i=1:PS_n
    L3=Phase(i,:)'; % High coherence's T
    Def_LP_X(:,i)=SVD(B13,L3);     % BX_p=L；X_p=INV(B)L；在这里L相当于方程的常数矩阵 
end
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\Def_LP_X','Def_LP_X');

clear;clc;
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\Def_LP_X.mat; 
n=18;
V_LP=Def_LP_X(1:n,:); % deformation rate，unti:mm/year
% VV1=Def_LP_X(n+1,:);
V11=Def_LP_X(n+1,:);
V22=Def_LP_X(n+2,:);
V_LP=V_LP';
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\V_LP','V_LP');
% save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\VV1','VV1');
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\V11','V11');
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\V22','V22');

clear;clc;
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\Def_LP_X.mat; 
n=18;
V33=Def_LP_X(n+3,:);
V44=Def_LP_X(n+4,:);
SigmaH1=Def_LP_X(n+5,:);
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\SigmaH1','SigmaH1');
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\V33','V33');
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\V44','V44');

clear;clc;
% 线性形变求解,当内存不足时，可先清空文件，然后重新打开成果文件。
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\V_LP.mat; 
TT=[0 22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
n=18;theta=39.2655;   %  incident angle
DT=zeros(1,n);
for i=1:n
    DT(i)=(TT(i+1)-TT(i))./365;     % ./365; % for multi-v model, the 15th differenced T
end
S_Los(:,1)=DT(1)*V_LP(:,1);  % LOS deformation value, v*t
for i=2:n
    S_Los(:,i)=S_Los(:,i-1)+DT(i)*V_LP(:,i);
end
S1=1000.*S_Los./(cosd(theta)); % vetical deformation
for x=1:n
    eval(['V2_',num2str(x),'=V_LP(:,',num2str(x),');'])
end
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\S1','S1');

clear;clc;
% 周期形变求解,当内存不足时，可先清空文件，然后重新打开成果文件。
% load K:\GAMMAUS1\TerraSARX3\output\EFDEM\VV1.mat VV1;
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\V11.mat V11; 
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\V22.mat V22; 
theta=39.2655;   %  incident angle
T=[22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
S2=1000.*(V11'.*sind(T.*2*pi/365)+V22'.*cosd(T.*2*pi/365))/(cosd(theta)); % vetical deformation
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\S2','S2');

clear;clc;
% 环境形变求解,当内存不足时，可先清空文件，然后重新打开成果文件。
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\V33.mat V33; 
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\V44.mat V44; 
theta=39.2655;   %  incident angle
Tem=[23.8375  23.625	21.375	19.125	18.125  17.125	15 	    12.875	13.437	14	    13.75	13.5	13.70	13.91	14.125	14.33	14.54  14.75]/30;
Pre=[0.8      0.92	   	5.7	     5.7    6.25	6.25	72	    72	    138.5	138.5	57.5	 57.5	21.5	21.5	21.5	34      34     34]/30;
S3=1000.*(V33'.*Tem+V44'.*Pre)/(cosd(theta));
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\S3','S3');

clear;clc;
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\S1.mat S1; 
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\S2.mat S2; 
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\S3.mat S3; 
load K:\GAMMAUS1\TerraSARX3\output\EFDEM\SigmaH1.mat SigmaH1; 
S=S1+S2+S3;
save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\S','S');
% for x=1:n
%     eval(['save(''K:\GAMMAUS1\TerraSARX3\output\EFDEM\V2\V2_',num2str(x),'.mat'',''V2_',num2str(x),''');']);
% end
% V2=zeros(PS_n,n);
% for i=1:n
%     V2(:,i)=eval(['V2_',num2str(i),'(:)']); 
% end
% save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\V2','V2');
% V2_mean=mean(V2,2);
% save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\V2','V2_mean');


% % 稀疏DEM反演
% HighlyDEM1=dem30_1+SigmaH1';
% % 生成稀疏矩阵
% H=zeros(7000,5000);x1=row_col_ct(:,1);y1=row_col_ct(:,2);
% H(sub2ind(size(H),x1,y1))=SigmaH1'; 
% % 整体DEM反演
% HighDEM1=finaldem1+H;

% save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\HighlyDEM1','HighlyDEM1');
% save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\HighDEM1','HighDEM1');
% save('K:\GAMMAUS1\TerraSARX3\output\EFDEM\H','H');
