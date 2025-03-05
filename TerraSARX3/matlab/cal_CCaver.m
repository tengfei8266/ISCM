clear;clc;
filepath=uigetdir();%选择cc文件夹的路径
cc_sum=0;
n=36;
for i=1:n %n是要读入的文件的个数
    load([filepath '\' 'cc' num2str(i) '.mat'])
    cc_sum=cc_sum+temp;
    display(num2str(i))
end
CCaver=cc_sum./n;
figure(1);imagesc(CCaver);axis image off;colormap(gray);
filepath2=uigetdir();%选择残差CCaver.mat文件的保存路径
save([filepath2 '\' 'CCaver.mat'],'CCaver');
figure(1);imagesc(CCaver);axis image off;
% exportgraphics(figure(1),[filepath2 '\' 'CCaver' '.tif'],'Resolution',300)