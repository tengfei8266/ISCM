 clear;clc;

load F:\GAMMAUS1\TerraSARX3\output\Data\Da.mat;
load F:\GAMMAUS1\TerraSARX3\output\CCaver\CCaver.mat;
load F:\GAMMAUS1\TerraSARX3\output\Data\mli.mat;

% figure(2);imagesc(mli);axis image off;
% 此处的范围值与DEM读取的值保持一致。
cc=CCaver(401:7600,101:6100);  
Da=Da(401:7600,101:6100);
mli=mli(401:7600,101:6100);
% % 掩膜赋值。
% cc(1:7000,1:1325)=1; 
% cc(4830:7000,1326:3006)=1;

figure(1);imagesc(Da);axis image off;colormap(gray)
figure(2);imagesc(mli);axis image off;colormap(gray)
figure(3);imagesc(cc);axis image off;colormap(gray)

for i=1:7200
    for j=1:6000
        if mli(i,j)<0.35||Da(i,j)>0.2||cc(i,j)<0.60
            
            cc(i,j)=0;
        end
    end
 end
% 
cc(isnan(cc)) = 1;
[x,y]=find(cc~=0);
row_col_ct=[x,y];
CT_n=size(row_col_ct,1);
hold on;
plot(row_col_ct(:,2),row_col_ct(:,1),'.','MarkerSize',1,'color','r');
% % 
save('F:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat');

clear;clc;
load F:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat row_col_ct;
n=size(row_col_ct,1);
x=row_col_ct(:,1);y=row_col_ct(:,2);

filepath='F:\GAMMAUS1\TerraSARX3\output\unw\';%文件夹的路径
   for i=1:36 %n是要读入的文件的个数
       load([filepath 'unw' num2str(i) '.mat'])
   end


 for j=1:36
     load([filepath 'unw' num2str(j) '.mat'])
     temp=temp(401:7600,101:6100);
     
     for i=1:n
        Num_mask_lungui(i,j)=temp(x(i),y(i));
     end
 end
 
save('F:\GAMMAUS1\TerraSARX3\output\highcp\Num_mask_lungui.mat');


load F:\GAMMAUS1\TerraSARX3\output\dem\finaldem1.mat;

final1=finaldem1;
 for j=1
     eval(['final',num2str(j),'=final',num2str(j),'(1:7200,1:6000);'])
     for ii=1:n
         eval(['dem30_1(ii,',num2str(j),')=final',num2str(j),'(x(ii),y(ii));'])
     end
 end

 save('F:\GAMMAUS1\TerraSARX3\output\dem\dem30_1.mat');
 

 