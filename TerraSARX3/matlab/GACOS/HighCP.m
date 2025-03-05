clear;clc;
load M:\GAMMAUS1\TerraSARX3\output\high\row_col_ct.mat row_col_ct;
n=size(row_col_ct,1);
x=row_col_ct(:,1);y=row_col_ct(:,2);

filepath='M:\GAMMAUS1\TerraSARX3\output\UTM\';%文件夹的路径
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
 
save('M:\GAMMAUS1\TerraSARX3\output\UTM\Num_mask_lungui.mat');


 

 