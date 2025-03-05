clear;clc;
load F:\GAMMAUS1\TerraSARX3\output\ICA\Phase.mat; 
load F:\GAMMAUS1\TerraSARX3\output\ICA\row_col_ct; 
Lambda=0.032 ; % center_range_slc  地理编码的文件中可以找到
%多速率形变项
TT=[0 22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
T=[0 22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264]./365;
%空间ICA
S=Phase';  %i*t时 为时间ICA   t*i时 为空间ICA
i=5;       %分离的个数;
[B,Z]=mFastICA(S,5);
W=pinv(B);
[ICs,SRs]=Sort_ICs(W',Z');
save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\W.mat W
save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\Z.mat Z
save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\ICs.mat ICs
save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\SRs.mat SRs

%原始数据检查
clear;clc;
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\W.mat W
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\Z.mat Z

% figure(4);plot(T,ICs(4,:))
  S11=(Z(1,:))'.*(W(:,1))';
  S22=(Z(2,:))'.*(W(:,2))';
  S33=(Z(3,:))'.*(W(:,3))';
  S44=(Z(4,:))'.*(W(:,4))';
  SS11=S11+S22+S33+S44;
  save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\SS11.mat SS11
%原始数据检查(内存不足时分开）
clear;clc;
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\W.mat W
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\Z.mat Z
  S55=(Z(5,:))'.*(W(:,5))';
%   S66=(Z(6,:))'.*(W(:,6))';
%   S77=(Z(7,:))'.*(W(:,7))';
%   S88=(Z(8,:))'.*(W(:,8))';
SS22=S55;
save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\SS22.mat SS22

%原始数据检查(内存不足时分开）
clear;clc;
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\SS11.mat SS11
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\SS22.mat SS22
SS1=[SS11 SS22];
save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\SS1.mat SS1

%分离后数据
clear;clc;
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\ICs.mat ICs
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\SRs.mat SRs
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\SS1.mat SS1
TT=[22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
  figure(1);plot(TT,ICs(1,:))
  figure(2);plot(TT,ICs(2,:))
  figure(3);plot(TT,ICs(3,:))
  figure(4);plot(TT,ICs(4,:))
%   figure(5);plot(TT,ICs(5,:))
%   figure(6);plot(TT,ICs(6,:))
%   figure(7);plot(TT,ICs(7,:))
%   figure(8);plot(TT,ICs(8,:))
  S1=SRs(:,1).*ICs(1,:);
  S2=SRs(:,2).*ICs(2,:);
  S3=SRs(:,3).*ICs(3,:);
  S4=SRs(:,4).*ICs(4,:);
  S11=S1+S2+S3+S4;
  save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S11.mat S11
  save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S1.mat S1
  save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S2.mat S2
  save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S3.mat S3
  save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S4.mat S4
%分离后数据
clear;clc;
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\ICs.mat ICs
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\SRs.mat SRs
  S5=SRs(:,5).*ICs(5,:);
%   S6=SRs(:,6).*ICs(6,:);
%   S7=SRs(:,7).*ICs(7,:);
%   S8=SRs(:,8).*ICs(8,:);
  S22=S5;
save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S22.mat S22
save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S5.mat S5
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S6.mat S6
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S7.mat S7
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S8.mat S8
clear;clc;
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S11.mat S11
load G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\S22.mat S22
SS=[S11 S22];
save G:\GAMMAUS1\TerraSARX3\output\ICA5\SICA\SS.mat SS
















% %时间ICA
% S=Phase;  %i*t时 为时间ICA   t*i时 为空间ICA
% i=4;       %分离的个数;
% [B,Z]=mFastICA(S,4);W=pinv(B);
% [ICs,SRs] = Sort_ICs(Z,W);
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\TICA\W.mat W
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\TICA\Z.mat Z
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\TICA\ICs.mat ICs
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\TICA\SRs.mat SRs
%   figure(1);plot(TT,ICs(1,:)')
%   figure(2);plot(TT,ICs(2,:)')
%   figure(3);plot(TT,ICs(3,:)')
%   figure(4);plot(TT,ICs(4,:)')
%   S11=W(:,1).*Z(1,:);
%   S22=W(:,2).*Z(2,:);
%   S33=W(:,3).*Z(3,:);
%   S44=W(:,4).*Z(4,:);
% SS1=S11+S22+S33+S44;
%   S1=SRs(:,1).*ICs(1,:);
%   S2=SRs(:,2).*ICs(2,:);
%   S3=SRs(:,3).*ICs(3,:);
%   S4=SRs(:,4).*ICs(4,:);
% SS=S1+S2+S3+S4;
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\TICA\S1.mat S1
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\TICA\S2.mat S2
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\TICA\S3.mat S3
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\TICA\S4.mat S4
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\TICA\SS1.mat SS1
% save G:\GAMMAUS1\TerraSARX3\output\ICA5\TICA\SS.mat SS
