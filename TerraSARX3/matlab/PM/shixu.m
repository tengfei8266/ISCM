clear,clc;
load L:\GAMMAUS1\TerraSARX3\output\high\row_col_ct row_col_ct;
load L:\GAMMAUS1\TerraSARX3\output\CTMDEM\SigmaH1 SigmaH1;
load L:\GAMMAUS1\TerraSARX3\output\dem\dem30_1.mat dem30_1; 
load L:\GAMMAUS1\TerraSARX3\output\dem\finaldem1.mat finaldem1;
load L:\GAMMAUS1\TerraSARX3\output\CTMDEM\S S;

% mat格式转为tiff,先上下翻转矩阵；
finaldem11=fliplr(finaldem1);
myimage=finaldem11;
myimage=flipud(myimage);
latlim=[32.59,32.71];
lonlim=[-116.99,-116.89];
R=georefcells(latlim,lonlim,size(myimage));
geotiffwrite('finaldem11.tif',myimage,R);

D=S(:,1);
HH=zeros(7200,6000);x1=row_col_ct(:,1);y1=row_col_ct(:,2);
HH(sub2ind(size(HH),x1,y1))=D';
HH=fliplr(HH);
figure(11);imagesc(HH);axis image off;colormap(jet);

% 稀疏DEM反演
CTMHighlyDEM1=dem30_1-SigmaH1';

save('L:\GAMMAUS1\TerraSARX3\output\CTMDEM\CTMHighlyDEM1','CTMHighlyDEM1');

% 生成稀疏矩阵
H=zeros(7200,6000);x1=row_col_ct(:,1);y1=row_col_ct(:,2);
H(sub2ind(size(H),x1,y1))=-SigmaH1';
H=fliplr(H);
figure(1);imagesc(H);axis image off;colormap(jet);

%DEM误差平滑插值拟合，去除森林和建筑物噪点
[x,y]=size(H);
%DEM误差平滑拟合
 % 给DEM矩阵外围加一圈，方便边界位置的处理
H=[zeros(x,1),H];  % 最左侧加一列 
H=[H,zeros(x,1)];  % 最左侧加一列 
H=[zeros(1,y+2);H];  % 最上侧加一列 
H=[H;zeros(1,y+2)];  % 最下侧加一列 
H_rigine=H; % 存储此时的矩阵，方便检查，查看是否有误。

% % 平滑处理
for t=1:40    % 平滑处理的重复次数，最多4O次
    for i=2:x+1
        for j=2:y+1
            if H(i,j)==0  % 意味着背景值或最外围一圈，跳过
                continue
            else
                a=sum([H(i-1,j),H(i+1,j),H(i,j-1),H(i,j+1)]==0);  % 计算目标栅格周围为0的栅格数
                % 计算中心点的高程,4/(4-a)会根据周围0的个数而变化，从而改变权重，默认不出现周围全是0的情况，a在0-3之间
                H(i,j)=0.5*H(i,j)+0.125*4/(4-a)*(H(i-1,j)+H(i+1,j)+H(i,j-1)+H(i,j+1));

            end
        end
    end
    if t==8||t==16||t==24||t==32||t==40     % 重复平滑处理8，16等次数后输出图像
        H_out=H;
        H_out(:,y+2)=[];
        H_out(:,1)=[];
        H_out(x+2,:)=[];
        H_out(1,:)=[];   %删除外围一圈
    end
end
z=isnan(H_out);
H_out(z)=0;
figure(3);imagesc(H_out);axis image off;colormap(jet);caxis([-50 50]);

HighDEM1=finaldem11+H_out;

% save('L:\GAMMAUS1\TerraSARX3\output\CTMDEM\HighDEM1','HighDEM1');
% save('L:\GAMMAUS1\TerraSARX3\output\CTMDEM\H','H');
% mat格式转为tiff,先上下翻转矩阵；

PMHighDEM11=HighDEM1;
figure(2);imagesc(PMHighDEM11);axis image off;colormap(jet);
myimage=PMHighDEM11;
myimage=flipud(myimage);
latlim=[32.59,32.71];
lonlim=[-116.99,-116.89];
R=georefcells(latlim,lonlim,size(myimage));
geotiffwrite('PMHighDEM11.tif',myimage,R);
