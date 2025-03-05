clear;clc;
load F:\GAMMAUS1\TerraSARX3\output\ICA\SICA\SS.mat;
load F:\GAMMAUS1\TerraSARX3\output\ICA\SICA\S1.mat;
load F:\GAMMAUS1\TerraSARX3\output\ICA\SICA\S2.mat;
load F:\GAMMAUS1\TerraSARX3\output\ICA\SICA\S3.mat;
load F:\GAMMAUS1\TerraSARX3\output\ICA\SICA\S4.mat;
load F:\GAMMAUS1\TerraSARX3\output\ICA\SICA\S5.mat;
load F:\GAMMAUS1\TerraSARX3\output\ICA\SICA\S6.mat;
PhaseICA=S1+S2+S3+S4+S5+S6;
save('F:\GAMMAUS1\TerraSARX3\output\ICDEM\PhaseICA','PhaseICA');

clear;clc;
load F:\GAMMAUS1\TerraSARX3\output\ICDEM\PhaseICA.mat;
load F:\GAMMAUS1\TerraSARX3\output\highcp\Num_mask_lungui.mat Num_mask_lungui; 
a=isnan(Num_mask_lungui);
Num_mask_lungui(a)=0;
HCPPhase=Num_mask_lungui;
PhaseICA=[PhaseICA(:,2) PhaseICA(:,3)-PhaseICA(:,1) PhaseICA(:,4)-PhaseICA(:,1) PhaseICA(:,4)-PhaseICA(:,2) PhaseICA(:,6)-PhaseICA(:,2) PhaseICA(:,4)-PhaseICA(:,3) PhaseICA(:,6)-PhaseICA(:,3) PhaseICA(:,7)-PhaseICA(:,3) PhaseICA(:,6)-PhaseICA(:,4) PhaseICA(:,9)-PhaseICA(:,4) PhaseICA(:,8)-PhaseICA(:,5) PhaseICA(:,9)-PhaseICA(:,5) PhaseICA(:,7)-PhaseICA(:,6) PhaseICA(:,9)-PhaseICA(:,6) PhaseICA(:,8)-PhaseICA(:,7) PhaseICA(:,10)-PhaseICA(:,7) PhaseICA(:,10)-PhaseICA(:,8) PhaseICA(:,11)-PhaseICA(:,8) PhaseICA(:,11)-PhaseICA(:,9) PhaseICA(:,15)-PhaseICA(:,9) PhaseICA(:,11)-PhaseICA(:,10) PhaseICA(:,15)-PhaseICA(:,10) PhaseICA(:,13)-PhaseICA(:,11) PhaseICA(:,15)-PhaseICA(:,11) PhaseICA(:,17)-PhaseICA(:,11) PhaseICA(:,13)-PhaseICA(:,12) PhaseICA(:,18)-PhaseICA(:,12) PhaseICA(:,14)-PhaseICA(:,13) PhaseICA(:,18)-PhaseICA(:,13) PhaseICA(:,15)-PhaseICA(:,14) PhaseICA(:,17)-PhaseICA(:,14) PhaseICA(:,17)-PhaseICA(:,15) PhaseICA(:,18)-PhaseICA(:,15) PhaseICA(:,17)-PhaseICA(:,16) PhaseICA(:,18)-PhaseICA(:,16) PhaseICA(:,18)-PhaseICA(:,17)];
ResPhase=HCPPhase-PhaseICA;
Res=mean(ResPhase,2);
save('F:\GAMMAUS1\TerraSARX3\output\ICDEM\PhaseICA','PhaseICA');
save('F:\GAMMAUS1\TerraSARX3\output\ICDEM\ResPhase','ResPhase');

load F:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat row_col_ct;
x1=row_col_ct(:,1);y1=row_col_ct(:,2);
T1=Res;
H1=zeros(7200,6000);   % 生成0值矩阵
H1(H1==0)=NaN;         % 改变矩阵为空值
H1(sub2ind(size(H1),x1,y1))=T1';  % 矩阵元素对应赋值
H1=fliplr(H1);          % 讲矩阵进行上下颠倒编码
Figure1=imagesc(H1);
set(Figure1,'alphadata',~isnan(H1));
axis image off;colormap(jet);caxis([-4 4]);