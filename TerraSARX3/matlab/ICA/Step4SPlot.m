clear,clc;
load F:\GAMMAUS1\TerraSARX3\output\ICA7\SICA\ICs.mat ICs; 
load F:\GAMMAUS1\TerraSARX3\output\ICA7\SICA\SRs.mat SRs;
load F:\GAMMAUS1\TerraSARX3\output\ICA7\SICA\SS.mat SS;
load F:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat row_col_ct;
load F:\GAMMAUS1\TerraSARX3\output\dem\finaldem1.mat finaldem1;

x1=row_col_ct(:,1);y1=row_col_ct(:,2);
% T11=SS(:,1);
% H11=zeros(7200,5000);
% H11(sub2ind(size(H11),x1,y1))=T11';
% figure(111);imagesc(H11);axis image off;
% 时间序列绘图
TT=[22 44 55 77 88 99 121 132 154 165 187 198 209 220 231 242 253 264];
  figure(1);plot(TT,ICs(1,:)')
  figure(2);plot(TT,ICs(2,:)')
  figure(3);plot(TT,ICs(3,:)')
%   figure(4);plot(TT,ICs(4,:)')
%   figure(5);plot(TT,ICs(5,:)')
%   figure(6);plot(TT,ICs(6,:)')
%   figure(7);plot(TT,ICs(7,:)')
%   figure(8);plot(TT,ICs(8,:)')
% 空间响应绘图

%
T1=SRs(:,1);
H1=zeros(7200,6000);   % 生成0值矩阵
H1(H1==0)=NaN;         % 改变矩阵为空值
H1(sub2ind(size(H1),x1,y1))=T1';  % 矩阵元素对应赋值
H1=fliplr(H1);          % 讲矩阵进行上下颠倒编码
Figure1=imagesc(H1);
set(Figure1,'alphadata',~isnan(H1));
axis image off;colormap(jet);caxis([-1 1]);

%
T2=SRs(:,2);
H2=zeros(7200,6000);
H2(H2==0)=NaN; 
H2(sub2ind(size(H2),x1,y1))=T2';
H2=fliplr(H2); 
Figure2=imagesc(H2); 
set(Figure2,'alphadata',~isnan(H2));
axis image off;colormap(jet);caxis([-0.7 1]);

%  
% T3=SRs(:,3);
% H3=zeros(7200,6000);
% H3(H3==0)=NaN; 
% H3(sub2ind(size(H3),x1,y1))=T3';
% H3=fliplr(H3); 
% Figure3=imagesc(H3); 
% set(Figure3,'alphadata',~isnan(H3));
% axis image off;colormap(jet);caxis([-1 1]);

 
% T4=SRs(:,4);
% H4=zeros(7200,6000);
% H4(H4==0)=NaN; 
% H4(sub2ind(size(H4),x1,y1))=T4';
% H4=fliplr(H4); 
% Figure4=imagesc(H4); 
% set(Figure4,'alphadata',~isnan(H4));
% axis image off;colormap(jet);caxis([-1 1]);

%  
%  
% T5=SRs(:,5);
% H5=zeros(7200,6000);
% H5(H5==0)=NaN; 
% H5(sub2ind(size(H5),x1,y1))=T5';
% H5=fliplr(H5); 
% Figure5=imagesc(H5); 
% set(Figure5,'alphadata',~isnan(H5));
% axis image off;colormap(jet);caxis([-1 1]);
% % 
% % %  
% T6=SRs(:,6);
% H6=zeros(7200,6000);
% H6(H6==0)=NaN; 
% H6(sub2ind(size(H6),x1,y1))=T6';
% H6=fliplr(H6); 
% Figure6=imagesc(H6); 
% set(Figure6,'alphadata',~isnan(H6));
% axis image off;colormap(jet);caxis([-1 1]);
% 
%  
% T7=SRs(:,7);
% H7=zeros(7200,6000);
% H7(H7==0)=NaN; 
% H7(sub2ind(size(H7),x1,y1))=T7';
% H7=fliplr(H7); 
% Figure7=imagesc(H7); 
% set(Figure7,'alphadata',~isnan(H7));
% axis image off;colormap(jet);caxis([-1 1]);

% %  
% T8=SRs(:,8);
% H8=zeros(7200,6000);
% H8(H8==0)=NaN; 
% H8(sub2ind(size(H8),x1,y1))=T8';
% H8=fliplr(H8); 
% Figure8=imagesc(H8); 
% set(Figure8,'alphadata',~isnan(H8));
% axis image off;colormap(jet);caxis([-1 1]);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % 空间响应的单点绘图
% PS_n=size(row_col_ct,1);
% for a=1
%     eval(['Subsidence',num2str(a),'=','SRs(:,a);']);
% 
% % x1=row_col_ps(:,2);y1=row_col_ps(:,1);
% % PS_n=size(row_col_ps,1);
% figure(a);
% eval(['Subsidence=Subsidence',num2str(a),';'])
% Subsidence_max=max(Subsidence);
% Subsidence_min=min(Subsidence);
%  
% % 导入底图:强度图
% 
% % AC=CCaver(1200:2000,700:1050);
% % M=log10(AC);clear Da;
% % figure(a);imagesc(M);axis image off;colormap(gray);
% 
% hold on;


% %%
% %间隔设置
% step=0.02;
% n=100;
% J=jet(n); 
% 
% %%
%     for i=1:PS_n 
%     if Subsidence(i)>1
%         plot(x1(i),y1(i),'.','MarkerEdgeColor',[J(1,1) J(1,2) J(1,3)],'MarkerSize',1);
%     end
%     if Subsidence(i)<=-1
%         plot(x1(i),y1(i),'.','MarkerEdgeColor',[J(n,1) J(n,2) J(n,3)],'MarkerSize',1);
%     end
%     max1=1;
%     for  j1=2:(n-1)
%         if Subsidence(i)>max1-step && Subsidence(i)<=max1
%           plot(x1(i),y1(i),'.','MarkerEdgeColor',[J(j1,1) J(j1,2) J(j1,3)],'MarkerSize',1);
%         end
%         max1=max1-step;
%     end      
% 
%    hold on;
% 
%     end
% end
% 
% axis off




