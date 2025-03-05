clc;clear;
final1 = freadbkB('F:\GAMMAUS1\TerraSARX2\file\DEM_GEO\20090808.dem.final',8000,'float32');
figure(1);imagesc(final1);axis image off;

finaldem1=final1(401:7600,101:6100);

figure(2);imagesc(finaldem1);axis image off;

% finaldem2=final1(4001:8000,101:2500);
% figure(3);imagesc(finaldem2);axis image off;
% 
% finaldem3=final1(1:4000,2501:5000);
% figure(4);imagesc(finaldem3);axis image off;
% 
% finaldem4=final1(4001:8000,2501:5000);
% figure(5);imagesc(finaldem4);axis image off;

% % mat格式转为tiff,先左右翻转矩阵；
% finaldem11=fliplr(finaldem1);
% imwrite(uint16(finaldem11),'finaldem11.tif','tif');
% 
% % imwrite(uint16(finaldem1),'finaldem11.tif','tif');
save('F:\GAMMAUS1\TerraSARX3\output\dem\final1','final1');
save('F:\GAMMAUS1\TerraSARX3\output\dem\finaldem1','finaldem1');
% save('G:\GAMMAchangsha\TerraSARX3\SBAS1F\output\dem\finaldem2','finaldem2');
% save('G:\GAMMAchangsha\TerraSARX3\SBAS1F\output\dem\finaldem3','finaldem3');
% save('G:\GAMMAchangsha\TerraSARX3\SBAS1F\output\dem\finaldem4','finaldem4');




