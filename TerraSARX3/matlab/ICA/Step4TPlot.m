clear,clc;
load K:\GAMMAchangsha\TerraSARX3\SBAS\output\ICA\TICA\ICs; 
load K:\GAMMAchangsha\TerraSARX3\SBAS\output\ICA\TICA\SRs;
load K:\GAMMAchangsha\TerraSARX3\SBAS\output\high\row_col_ct;
% 时间序列绘图
  TT=[22 44 66 110 132 154];
  figure(1);plot(TT,ICs(1,:)')
  figure(2);plot(TT,ICs(2,:)')
  figure(3);plot(TT,ICs(3,:)')
  figure(4);plot(TT,ICs(4,:)')
% 空间响应绘图
x1=row_col_ct(:,1);y1=row_col_ct(:,2);
%
T1=SRs(:,1);
H1=zeros(4000,2400);
H1(sub2ind(size(H1),x1,y1))=T1';
figure(11);imagesc(H1);axis image off;
%
T2=SRs(:,2);
H2=zeros(4000,2400);
H2(sub2ind(size(H2),x1,y1))=T2';
figure(22);imagesc(H2);axis image off;
%  
T3=SRs(:,3);
H3=zeros(4000,2400);
H3(sub2ind(size(H3),x1,y1))=T3';
figure(33);imagesc(H3);axis image off;
%
T4=SRs(:,4);
H4=zeros(4000,2400);
H4(sub2ind(size(H4),x1,y1))=T4';
figure(44);imagesc(H4);axis image off;
  
