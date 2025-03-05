
%数据预处理
[Fname,T,B,L,H,dN,dE,dU] = readdata();%读取原始数据
[T]=transdate_o(T);%日期转换
%提取2006.1.1-2015.12.31之间的U向时间序列数据 3652days
N1=datenum('2010 1 1');
N2=datenum('2019 12 31');
X=zeros(3652,length(T))*nan;
for i=1:length(T)
    u=dU{i};t=T{i};
    u=u(t>=N1&t<=N2);
    t=t(t>=N1&t<=N2);
    p=t-N1+1;
    x=zeros(3652,1)*nan;
    x(p)=u;
    X(:,i)=x;
    clear u t p x;
end
%每个台站数据缺失统计
p=isnan(X);
q=abs(1-p);
A=cell(1,358);
for i=1:length(T)
    A{i}=X(:,i);
    A{i}(p(:,i))=[];
    D(i)=length(A{i})/3652;
    q(:,i)=q(:,i)*i;
end
D=D';    
q(q==0)=nan;    
%画数据统计图（358个站点）
plot(q(:,1:358),'k','linewidth',3);
axis([1,3890,0,359])
set(gca,'XTick',[1,365,731,1096,1461,1826,2192,2557,2922,3287,3652],'ytick',[])
set(gca,'XTickLabel',{'2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016'})
set(gca,'fontweight','bold');
xlabel('Time/year','fontweight','bold');
box on;
for i=1:358
    text(3656,i,[num2str(D(i)*100,2),'%'],'fontweight','bold');
    text(0,i,num2str(Fname{i}),'horizontalalignment','right','fontweight','bold');
end

ylabel('Station Numbers','fontweight','bold');
%筛选大于数据长度90%的数据
D1=D;
for i = 1:358
    if D1(i) <0.90;
        D1(i) = 0;
    elseif D1(i)>=0.90;
        D1(i) = D(i);
    end
end
D1(D1==0)=nan; p1=isnan(D1);q1=abs(1-p1);

X1=X;Sname=Fname;
 for i=1:358
    if pp1(i,1)==1;
        X1(:,i)=X1(:,i);
        Sname{i}=Fname{i};
    elseif pp1(i,1)==0;
        X1(:,i)=zeros(3652,1);
        Sname{i}=[];
    end
 end
Sname(cellfun('isempty',Sname))=[];%清除Sname中包含的[]
X2=X1;X2(:,all(X2==0,1))=[];%清除X1中包含的0列


X0=X;clear X;
%画45各站点数据长度图
plot(q(:,1:45),'k','linewidth',3);
axis([1,3890,0,46])
set(gca,'XTick',[1,365,731,1096,1461,1826,2192,2557,2922,3287,3652],'ytick',[])
set(gca,'XTickLabel',{'2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016'})
set(gca,'fontweight','bold');
xlabel('Time/year','fontweight','bold');
box on;
for i=1:45
    text(3656,i,[num2str(D(i)*100,2),'%'],'fontweight','bold');
    text(0,i,num2str(Fname{i}),'horizontalalignment','right','fontweight','bold');
end


for i=1:13
    text(2400,i,[num2str(D(i)*100,2),'%'],'fontweight','bold');
    text(0,i,num2str(Sname_XZ{i}),'horizontalalignment','right','fontweight','bold');
end

% 坐标转换 XYZ用来计算ATML
fid = fopen('45BLH.txt');
data = textscan(fid,'%s%f%f%f'); 
B=data{2};
L=data{3};
H=data{4};
Sname=data{1};
for i=1:length(B)
    [X(i),Y(i),Z(i)]=geo2xyz(B(i),L(i),H(i));
end
fclose(fid);  

Time=2006+1/365.25:1/365.25:2006+3652/365.25;
T_offset=cell(1,62);
T_offset{4}(1,1)=2006+2845/365.25;
T_offset{4}(2,1)=2006+2949/365.25;
T_offset{8}=2006+3189/365.25;
T_offset{13}=2006+1460/365.25;
T_offset{14}(1,1)=2006+2860/365.25;
T_offset{14}(2,1)=2006+2922/365.25;
T_offset{17}=2006+2158/365.25;
T_offset{19}=2006+2734/365.25;
T_offset{21}=2006+3185/365.25;
T_offset{22}=2006+2045/365.25;
T_offset{24}=2006+1468/365.25;
T_offset{31}=2006+2747/365.25;
T_offset{32}=2006+2763/365.25;
T_offset{34}(1,1)=2006+2408/365.25;
T_offset{34}(2,1)=2006+2861/365.25;
T_offset{34}(3,1)=2006+3480/365.25;
T_offset{39}=2006+293/365.25;
T_offset{56}(1,1)=2006+3149/365.25;
T_offset{56}(2,1)=2006+3194/365.25;
T_offset{59}(1,1)=2006+1910/365.25;
T_offset{59}(2,1)=2006+2732/365.25;
T_offset{60}=2006+1459/365.25;

T_offset{2}=2006+1946/365.25;
T_offset{3}(1,1)=2006+3149/365.25;
T_offset{3}(2,1)=2006+3193/365.25;
T_offset{5}(1,1)=2006+91/365.25;
T_offset{5}(2,1)=2006+2407/365.25;
T_offset{8}(1,1)=2006+1472/365.25;
T_offset{8}(2,1)=2006+2032/365.25;
T_offset{10}(1,1)=2006+1565/365.25;
T_offset{10}(2,1)=2006+2464/365.25;
T_offset{13}(1,1)=2006+2408/365.25;
T_offset{13}(2,1)=2006+2861/365.25;
T_offset{14}=2006+2761/365.25;
T_offset{15}=2006+2746/365.25;
T_offset{17}=2006+360/365.25;
T_offset{19}=2006+1539/365.25;
T_offset{21}=2006+258/365.25;
T_offset{22}=2006+1462/365.25;
T_offset{24}=2006+1290/365.25;
T_offset{27}=2006+2045/365.25;
T_offset{28}=2006+3178/365.25;
T_offset{29}(1,1)=2006+2642/365.25;
T_offset{29}(2,1)=2006+2678/365.25;
T_offset{29}(3,1)=2006+2709/365.25;
T_offset{33}=2006+2761/365.25;
T_offset{34}=2006+72/365.25;
T_offset{35}(1,1)=2006+2425/365.25;
T_offset{35}(2,1)=2006+3266/365.25;
T_offset{40}=2006+1011/365.25;
T_offset{43}=2006+3128/365.25;
T_offset{41}=2006+2412/365.25;
T_offset{45}(1,1)=2006+3149/365.25;
T_offset{45}(2,1)=2006+3168/365.25;

T_offset(cellfun('isempty',T_offset))=num2cell(0);%T_offset中的【】变成0
%去除线性项和阶跃项
for i=1:45
    [X_detrend(:,i) fitline(:,i)]=fitmodel(Time,X_(:,i),T_offset{i});
end





%画简易版cGPS台站点位图
load('/Users/maxiaojun/Desktop/TW/TW_mat/SRs.mat')
for i=1:5
    mapshow(readL); hold on;
    set(gca,'fontweight','bold','color',[0.83,0.82,0.78]);
    axis([119.67,122.5,21.7,25.8]);
    scatter(L,B);box on;grid on;
end
for i=1:44
    D1(i)=D(i)/sum(D);
end


%m_map画cGPS台站点位图(代码运行两次)
figname='shadedrelief4';
m_proj('mercator','long',[118 124],'lat',[21 26]);
[ELEV,LON,LAT]=m_etopo2([118 124 21 26]);
caxis([-5000 3000]) % %caxis要放在colormap之前，colormap要放在m_shadedrelief之前
colormap([m_colmap('blues',200);m_colmap('gland',120)])% 使用两种色标区分陆地海洋，范围比例要和caxis匹配 
hc=colorbar;
set(get(hc,'title'),'string','Elevation(m)')	%  色标单位
m_shadedrelief(LON(1,:),LAT(:,1),ELEV) 
m_etopo2('shadedrelief','lightangle',45,'gradient',2);
m_gshhs('ic','color','k')
m_grid('box','fancy','tickdir','in','gridlines','no','fontsize',14)
set(gcf,'position',[100 100 800 500])
%地图上标注每个站点的名字
m_scatter(L,B);
for i=1:38
    m_text(L(i),B(i),Sname{i},'vertical','top',)
end



m_proj('mercator','lon',[73 140],'lat',[10 54]);
m_gshhs('hb1','linewidth',2.0,'color','k');% 陆地国界线
m_gshhs('hc1','linewidth',2.0,'color','k');
m_grid('linestyle','none','tickdir','in','linewidth',3,'fontname','Times New Roman','fontweight','bold','fontsize',16);




















%PCA PRINCIPLE COMPONENTS ANALYSIS
X1=X_detrend_intp'; %X1(38*3652)
[M,T] = size(X1);
average= mean(X1')';
for i=1:M
X1(i,:)=X1(i,:)-average(i)*ones(1,T);
end
C = X1*X1.';
[E,D] = eig(C);
D=diag(D);
D = D(end:-1:1)/max(D);
E = E(:,end:-1:1)';
P=E(1:10,:);
X11=P*X1;



%主成分单独贡献率和累积贡献率
for i=1:length(D)
   D1(i)=D(i)/sum(D);%D(38*1)
end%单独贡献率
D2=zeros(1,length(D));
for i=1:length(D)
    D2(i)=sum(D1(1,1:i));
end%累积贡献率




%贡献率画单坐标轴折线图
x=1:1:10;
y=[0.179721928941714,0.0835900187123608,0.0551681793485837,0.0499787741246863,0.0419684781326573,0.0372262942391094,0.0358279860263477,0.0331868427683266,0.0303749966731335,0.0298179078870903]; %a数据y值
plot(x,y,'-*b'); 
axis([0,11,0,0.2]) 
set(gca,'XTick',[0:1:11]) 
set(gca,'YTick',[0:0.1:0.2])
for i=1:10
   text(x(i),y(i),num2str(y(i)));
end




%贡献率画双坐标轴折线图
x=1:1:10;
y1=D1(1:10);
y2=D2(1:10);
figure;
[AX,h1,h2]=plotyy(x,y1,x,y2);
set(AX(1),'ycolor','r')
set(AX(2),'ycolor','b')
set(h1,'linestyle','-','marker','diamond','color','r')
set(h2,'linestyle','-','marker','pentagram','color','b')
for i=1:10
   text(AX(1),x(i)+0.1,y1(i)+0.002,num2str(y1(i),3));
   text(AX(2),x(i)+0.12,y2(i)-0.004,num2str(y2(i),3));
end
legend('Variance ','Cumulative');
title('Principal Component contribution rate');




%简单画台湾地图轮廓以及空间响应-------(对照数据SRs.mat)
load('/Users/maxiaojun/Desktop/TW/TW_mat/SRs.mat')
for i=1:5
    subplot(1,5,i);mapshow(readL); hold on;
    set(gca,'fontweight','bold','color',[0.83,0.82,0.78]);
    title(['SR',num2str(i)],'fontweight','bold');
    axis([119.67,122.5,21.7,25.8]);
    for j=1:38
        if SRs(i,j)>0;
            quiver(L(j),B(j),0,0.45*SRs(i,j)','color','r','maxheadsize',0.7,'linewidth',1.0);hold on
        elseif SRs(i,j)<0;
            quiver(L(j),B(j),0,0.45*SRs(i,j)','color','g','maxheadsize',0.7,'linewidth',1.0);hold on
        end
    end
    box on;grid on;
    quiver(119.832,22.127,0,0.45,'color','r','maxheadsize',0.7,'linewidth',1.0);
    quiver(119.832,22.127,0,-0.45,'color','g','maxheadsize',0.7,'linewidth',1.0);
    text(119.832,22.117,'100%','fontweight','bold');
end
%箭头红色朝上绿色朝下
% for i=1:38
%     if P1(4,i)>0;
%         quiver(L(i),B(i),0,P1(4,i)','color','r');box on;grid on; 
%     elseif P1(4,i)<0;
%         quiver(L(i),B(i),0,P1(4,i)','color','g');box on;grid on; 
%     end
% end


%m_map 画台湾地图以及空间响应(代码要运行两次)
for i=1:5
    subplot(2,3,i)
    figname='shadedrelief3';
    m_proj('mercator','long',[119.4 123],'lat',[21.7 25.8]);
    caxis([-5000 1000]);	% 设置显示高程范围，只显示海洋部分
    colormap((m_colmap('blues',200)));  % 设置色标为m_map自带的m_colmap('blues')
   %  色标单位
    m_etopo2('shadedrelief','lightangle',-45,'gradient',1); % 设置光源方位角为45°,坡角阈值为1
    m_gshhs('ic','patch',[.6 .6 .6]);	% 设置中等分辨率海岸线填充
    m_grid('box','fancy','gridlines','no','fontsize',10,'xlabeldir','end');  % 轴设置
    set(gcf,'position',[100 100 300 800]) %图形大小设置
    title(['PC',num2str(i)],'fontweight','bold','fontsize',12);
    for j=1:38
        if P1(i,j)>0;
            m_quiver(L(j),B(j),0,0.45*P1(i,j)','color','r','maxheadsize',0.7,'linewidth',1.0);hold on
        elseif P1(i,j)<0;
            m_quiver(L(j),B(j),0,0.45*P1(i,j)','color','g','maxheadsize',0.7,'linewidth',1.0);hold on
        end
    end
    m_quiver(119.822,25.1,0,0.45,'color','r','maxheadsize',0.7,'linewidth',0.5);
    m_quiver(119.822,25.1,0,-0.45,'color','g','maxheadsize',0.7,'linewidth',0.5);
    m_text(119.832,25.09,'100%','fontweight','bold');
end
hc=colorbar;set(get(hc,'title'),'string','Elevation(m)')%最后一个图加色条




%画主成分时间序列图(final)_plot
ha = tight_subplot(5, 2,[0.02 0.06], [0.1 0.1], [0.1 0.1]);
for i=1:5
    axes(ha(i));
    if mod(i,2)==0;
        plot(Components(i,:),'color',[0.31,0.31,0.31]);axis([0,3652,-20,20]);
        set(ylabel(['IC',num2str(i/2),'/mm']),'fontweight','bold');
    else
        plot(Components(i,:),'color',[0.31,0.31,0.31]);axis([0,3652,-20,20]);
        set(ylabel(['PC',num2str(round(i/2)),'/mm']),'fontweight','bold');
    end
end
set(ha(1:8),'XTickLabel','','fontweight','bold');
set(ha(9:10),'XTick',[1,365,731,1096,1461,1826,2192,2557,2922,3287,3652]);
set(ha(9:10),'XTickLabel',{'2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016'},'fontweight','bold');

%画主成分时间序列图(final)_scatter
ha = tight_subplot(5, 2,[0.02 0.06], [0.1 0.1], [0.1 0.1]);
for i=1:10
    axes(ha(i));
    if mod(i,2)==0;
        scatter(1:3652,Components(i,:),'sizedata',3,'MarkerEdgeColor','k','MarkerFaceColor','k');axis([0,3652,-20,20]);
        set(ylabel(['IC',num2str(i/2),'/mm']),'fontweight','bold');
    else
        scatter(1:3652,Components(i,:),'sizedata',3,'MarkerEdgeColor','k','MarkerFaceColor','k');axis([0,3652,-20,20]);
        set(ylabel(['PC',num2str(round(i/2)),'/mm']),'fontweight','bold');
    end
end
set(ha(1:8),'XTickLabel','','fontweight','bold');
set(ha(9:10),'XTick',[1,365,731,1096,1461,1826,2192,2557,2922,3287,3652]);
set(ha(9:10),'XTickLabel',{'2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016'},'fontweight','bold');


[ICs,SRs] = Sort_ICs(IC,SR);















%求PCA、ICA时空滤波后整体平均RMS变化情况-------(对照数据RMS.mat)
load('mxj12.mat', 'X_detrend_intp')%加载预处理完后并插值的数据
[B Z]=mFastICA(X_detrend_intp',5);%B=W'*Q*P Z=B*X （P为主成分对应的空间响应、B的逆是独立成分对应的空间响应、Z是独立成分、X是X_detrend_intp经过中心化后的数据）
SRs_IC=pinv(B);%SRs_IC是独立成分对应的空间响应，SRs_IC（38*5）的矩阵、B(5*38)的矩阵、Z(5*3652)的矩阵
CME_IC=SRs_IC(:,1)*Z(1,:);%CME_IC是独立分量计算的共模误差
%PCA中原始数据中心化过程
X=X_detrend_intp';[M,T] = size(X); %获取输入矩阵的行/列数，行数为观测数据的数目，列数为采样点数      
 average= mean(X')';  %均值
 for i=1:M
     X(i,:)=X(i,:)-average(i)*ones(1,T); %X是X_detrend_intp经过中心化后的数据 X(38*3652)、X_detrend_intp(3652*38)
 end
%PCA降维
P = pca(X,5,'eig');%P是前五个主成分对应的空间响应（SRs）
X1=P*X;%X1是前五个主成分(PCs)
SRs_PC=P;%SRs_PC(5*38)
CME_PC=SRs_PC(1,:)'*X1(1,:);%CME_PC是第一主分量计算的共模误差
X_IC=X-CME_IC;%X_IC是ICA滤波后的时间序列
X_PC=X-CME_PC;%X_IC是PCA滤波后的时间序列
for i=1:38
    RMS_X(i)=norm(X(i,:))./sqrt(3652);%求原始序列的RMS值
end
RMS(1,1)=mean(RMS_X,2);%求原始序列的平均RMS值
for i=1:38
    RMS_IC(i)=norm(X_IC(i,:))./sqrt(3652);%求ICA滤波后的时间序列的RMS值
end
RMS(2,1)=mean(RMS_IC,2);%求ICA滤波后的时间序列的平均RMS值
for i=1:38
    RMS_PC(i)=norm(X_PC(i,:))./sqrt(3652);%求PCA滤波后的时间序列的RMS值
end
RMS(3,1)=mean(RMS_PC,2);%求PCA滤波后的时间序列的平均RMS值



%计算Ratio值(老师的方法)
[ICs,SRs] = Sort_ICs(IC,SR);
%画图
bar(Ratio,0.5,'facecolor','black');
axis([0.5 5.5 0 10]);
set(gca, 'xticklabels', {'GPS-IC1', 'GPS-IC2', 'GPS-IC3', 'GPS-IC4', 'GPS-IC5'}, 'fontweight','bold', 'Fontsize', 12); 

%计算Ratio值(我的方法)
X00=X_detrend_intp';
for i=1:39
    X01(i)=norm(X00(i,:)).^2;
end
re5=c1(5,:)'*Z1(5,:);
for i=1:39
    re_IC5(i)=norm(re5(i,:)).^2;
end
for i=1:39
    Ratio_IC5(i)=10*(log10(X01(i))-log10(re_IC5(i)));
end

Ratio_IC(1,:)=Ratio_IC1;Ratio_IC(2,:)=Ratio_IC2;
Ratio_IC(3,:)=Ratio_IC3;Ratio_IC(4,:)=Ratio_IC4;Ratio_IC(5,:)=Ratio_IC5;
for i=1:5
    Ratio(i)=mean(Ratio)
end




%sICA处理流程
[B Z]=mFastICA(X_detrend_intp,5);%X_detrend_intp（3652*38）B(5*3652) Z(5*38)
 SRs_IC=pinv(B);%SRs_IC(3652*5)
 SR=Z';%SR(38*5)
 IC=SRs_IC';%IC(5*3652)
 [ICs,SRs] = Sort_ICs(IC,SR);%ICs(5*3652) SRs(38*5)
 load('mxj15.mat', 'B')
 SRs=SRs';%SRs(5*38)
 for i=1:5
    subplot(1,5,i);mapshow(readL); hold on;
    set(gca,'fontweight','bold','color',[0.83,0.82,0.78]);
    title(['IC',num2str(i)],'fontweight','bold');
    axis([119.67,122.5,21.7,25.8]);
    for j=1:38
        if SRs(i,j)>0;
            quiver(L(j),B(j),0,0.45*SRs(i,j)','color','r','maxheadsize',0.7,'linewidth',1.0);hold on
        elseif SRs(i,j)<0;
            quiver(L(j),B(j),0,0.45*SRs(i,j)','color','g','maxheadsize',0.7,'linewidth',1.0);hold on
        end
    end
    box on;grid on;
    quiver(119.832,22.127,0,0.45,'color','r','maxheadsize',0.7,'linewidth',1.0);
    quiver(119.832,22.127,0,-0.45,'color','g','maxheadsize',0.7,'linewidth',1.0);
    text(119.832,22.117,'100%','fontweight','bold');
end
 
 clear B IC ICs i j SR SRs SRs_IC Z;
 
 
 
 %计算ATML
 [Filename,pathname]=uigetfile({'*','data files'},'Importing','MultiSelect', 'on');
cnt=0;
for i=1:length(Filename)
    fid = fopen(char(Filename(i)));
    for j=1:490
        headinfo=fgetl(fid);
    end

    data = textscan(fid,'%s%f%f%f%s%s%f%f%f'); 
    
    ATMLu=data{7};
    ATMLu(end)=[];
    
    for k=1:length(ATMLu)/38;
        ATML_o(1:38,cnt+k)=ATMLu((k-1)*38+1:k*38);
    end
    a(i)=length(ATMLu)/38;
    cnt=cnt+length(ATMLu)/38;
    fclose(fid);
end
ATML_o=ATML_o';


%画空间响应
for i=1:10
    subplot(2,5,i)
    figname='shadedrelief3';
    m_proj('mercator','long',[119.4 123],'lat',[21.7 25.8]);
    [ELEV,LON,LAT]=m_etopo2([119.4 123 21.7 25.8]);
    caxis([-5000 3000]);	% 设置显示高程范围，只显示海洋部分
    colormap([m_colmap('blues',200);m_colmap('gland',120)])  % 设置色标为m_map自带的m_colmap('blues')
   %  色标单位
    m_shadedrelief(LON(1,:),LAT(:,1),ELEV) 
    m_etopo2('shadedrelief','lightangle',-45,'gradient',10); % 设置光源方位角为45°,坡角阈值为1
    m_gshhs('ic','color','k');	% 设置中等分辨率海岸线填充
    m_grid('box','fancy','gridlines','no','fontsize',10,'xlabeldir','end');  % 轴设置
    set(gcf,'position',[100 100 300 800]) %图形大小设置
    if i<=5;
        m_text(['PC',num2str(i)],'fontweight','bold','fontsize',12);
    else
        title(['IC',num2str(i-5)],'fontweight','bold','fontsize',12);
    end
    for j=1:44
        if Y(i,j)>0;
            m_quiver(L(j),B(j),0,0.8*Y(i,j)','color','r','maxheadsize',0.7,'linewidth',1.0);hold on
        elseif P1(i,j)<0;
            m_quiver(L(j),B(j),0,0.8*Y(i,j)','color','g','maxheadsize',0.7,'linewidth',1.0);hold on
        end
    end
    m_quiver(119.822,25.1,0,0.8,'color','r','maxheadsize',0.7,'linewidth',0.5);
    m_quiver(119.822,25.1,0,-0.8,'color','g','maxheadsize',0.7,'linewidth',0.5);
    m_text(119.832,25.09,'100%','fontweight','bold');
end






for i=1:10
    subplot(2,5,i)
    m_proj('mercator','long',[119.3 122.25],'lat',[21.75 25.75]);
    [CS,CH]=m_etopo2('contourf',[-5000:500:0 250:250:3000],'edgecolor','none');
    m_grid('linestyle','none','tickdir','out','linewidth',3);
    colormap([ m_colmap('blues',80); m_colmap('gland',48)]);
    brighten(.5);colormap('bone');m_gshhs_h('color','k');caxis([-8000,4000]);
     if i<=5;
        m_text(121.68,25.6,['PC',num2str(i)],'fontweight','bold','fontsize',16);
    else
        m_text(121.68,25.6,['IC',num2str(i-5)],'fontweight','bold','fontsize',16);
    end
    for j=1:44
        if Y(i,j)>0;
            m_quiver(L(j),B(j),0,0.55*Y(i,j)','color','r','maxheadsize',0.7,'linewidth',0.5);hold on
        elseif Y(i,j)<0;
            m_quiver(L(j),B(j),0,0.55*Y(i,j)','color','g','maxheadsize',0.7,'linewidth',0.5);hold on
        end
    end
    m_quiver(119.5,22.26,0,0.55,'color','r','maxheadsize',0.7,'linewidth',0.5);
    m_quiver(119.5,22.26,0,-0.55,'color','g','maxheadsize',0.7,'linewidth',0.5);
    m_text(119.5,22.26,'100%','fontweight','bold');
end



for i=1:10
     subplot(2,5,i);
     m_proj('mercator','long',[119.3 122.25],'lat',[21.75 25.75]);
     [CS,CH]=m_etopo2('contourf',[-5000:500:0 250:250:3000],'edgecolor','none');
     m_grid('linestyle','none','tickdir','out','linewidth',3,'fontname','Times New Roman','fontweight','bold','fontsize',16);
     colormap([ m_colmap('blues',80); m_colmap('gland',48)]);
     brighten(.5);colormap('bone');m_gshhs_h('color','k');caxis([-8000,4000]);
     if i<=5;
         m_text(121.68,25.6,['PC',num2str(i)],'fontweight','bold','fontsize',16);
     else
         m_text(121.68,25.6,['IC',num2str(i-5)],'fontweight','bold','fontsize',16);
     end
     for j=1:44
         if Y(i,j)>0;
             m_quiver(L(j),B(j),0,0.65*Y(i,j)','color','r','maxheadsize',0.7,'linewidth',2.5);hold on
         elseif Y(i,j)<0;
             m_quiver(L(j),B(j),0,0.65*Y(i,j)','color','g','maxheadsize',0.7,'linewidth',2.5);hold on
         end
     end
     m_quiver(119.5,22.26,0,0.65,'color','r','maxheadsize',0.7,'linewidth',2.5);
     m_quiver(119.5,22.26,0,-0.65,'color','g','maxheadsize',0.7,'linewidth',2.5);
     m_text(119.5,22.26,'100%','fontweight','bold','fontsize',16);
end


fid = fopen('CWB_BLH.txt');
data = textscan(fid,'%s%f%f%f');
B=data{2};
L=data{3};
H=data{4};
Sname=data{1};fclose(fid);
fid = fopen('IES_BLH.txt');
data = textscan(fid,'%s%f%f%f');
B=data{2};
L=data{3};
H=data{4};
Sname=data{1};fclose(fid);
fid = fopen('CGS_BLH.txt');
data = textscan(fid,'%s%f%f%f');
B=data{2};
L=data{3};
H=data{4};
Sname=data{1};fclose(fid);
fid = fopen('NTU_BLH.txt');
data = textscan(fid,'%s%f%f%f');
B=data{2};
L=data{3};
H=data{4};
Sname=data{1};fclose(fid);
load('/Users/maxiaojun/Desktop/zlj_o.mat/zlj_o08.mat', 'L')
load('/Users/maxiaojun/Desktop/zlj_o.mat/zlj_o08.mat', 'B1')
m_proj('mercator','long',[119.3 122.25],'lat',[21.75 25.75]);[CS,CH]=m_etopo2('contourf',[-5000:500:0 250:250:3000],'edgecolor','none');
m_grid('linestyle','none','tickdir','out','linewidth',3);
colormap([ m_colmap('blues',80); m_colmap('gland',48)]);
brighten(.5);
ax=m_contfbar(1,[.5 .8],CS,CH);
title(ax,{'Level/m',''});
h1=m_scatter(CWB_L,CWB_B,'o','filled','MarkerFaceColor','k','MarkeredgeColor','k');
h2=m_scatter(IES_L,IES_B,'o','filled','MarkerFaceColor','g','MarkeredgeColor','g');
h3=m_scatter(CGS_L,CGS_B,'o','filled','MarkerFaceColor','b','MarkeredgeColor','b');
h4=m_scatter(NTU_L,NTU_B,'o','filled','MarkerFaceColor','c','MarkeredgeColor','c');
h5=m_scatter(L,B1,'o','filled','MarkerFaceColor','r','MarkeredgeColor','r');
hl=legend([h1 h2 h3 h4 h5],'CWB','IES','CGS','NTU','44 test stations');
set(hl,'box','off');
m_northarrow(121.75,22.5,.4,'type',4);
m_ruler([.68 0.98],.1,4,'tickdir','out','ticklen',0.005);
m_ruler([.68 0.98],.105,5,'tickdir','out','ticklen',0.005);
m_ruler([.68 0.98],.11,5,'tickdir','out','ticklen',0.005);
m_ruler([.68 0.98],.115,5,'tickdir','out','ticklen',0.005);
m_ruler([.68 0.98],.115,4,'tickdir','out','ticklen',0.005);
m_ruler([.68 0.98],.11,4,'tickdir','out','ticklen',0.005);
hl=legend([h1 h2 h3 h4 h5],'CWB','IES','CGS','NTU','44 test stations');
set(hl,'box','off');
hold on;m_scatter(119.5923767220,23.6552730110);
hold on;m_scatter(119.5637482890,23.5652040330);
hl=legend([h1 h2 h3 h4 h5],'CWB','IES','CGS','NTU','44 test stations');
set(hl,'box','off');

clear,clc;
load H:\Foshan-20\lungui\1126\matlab-l\row_col_ps;
load H:\Foshan-20\lungui\1126\matlab-l\\S;
load H:\Foshan-20\lungui\1126\matlab-l\da\AC_aver;

x1=row_col_ps(:,2);y1=row_col_ps(:,1);
PS_n=size(row_col_ps,1);
for a=1:38
    eval(['Subsidence',num2str(a),'=','S(:,a)']);

% x1=row_col_ps(:,2);y1=row_col_ps(:,1);
% PS_n=size(row_col_ps,1);
eval(['Subsidence=Subsidence',num2str(a),';'])
Subsidence_max=max(Subsidence);
Subsidence_min=min(Subsidence);
% 导入底图:强度图
row=165;col=1552;
 fid1 = fopen('H:\Xinlihuagong_export\XLHG20181129_LXB_SBAS_processing\work\work_super_master\SM_sentinel1_11_20160117_63518566_IW_SIW1_A_VV_cut_slc_list_7_pwr','r');
 [mli1,count1] = fread(fid1,[row,col],'float32');
 fclose(fid1);
 mli1=abs(mli1');M=log10(mli1);
 M=M(655:1185,506:668);
AC=AC_aver;%(1:1800,1:2600);
M=log10(AC);clear ACaver;
figure(a);imagesc(M);axis image off;colormap(gray);
% figure(1);%imagesc(M);axis image off;colormap(gray);
hold on;  

J=jet(29); 
n_1=0;n_2=0;n_3=0;n_4=0;n_5=0;n_6=0;n_7=0;n_8=0;n_9=0;n_10=0;n_11=0;n_12=0;n_13=0;n_14=0;n_15=0;n_16=0;n_17=0;n_18=0;n_19=0;n_20=0;
n_21=0;n_22=0;n_23=0;n_24=0;n_25=0; n_26=0;n_27=0;n_28=0;n_29=0;n_30=0;n_31=0;n_32=0;n_33=0;n_34=0;n_35=0;n_36=0;n_37=0;n_38=0;n_39=0; n_40=0;
n_41=0;n_42=0;n_43=0;n_44=0;n_45=0;n_46=0;n_47=0;n_48=0;n_49=0;n_50=0;n_51=0;n_52=0;
for i=1:PS_n
%       if  Subsidence(i)<=-108-10
%            plot(x1(i),y1(i),'.','MarkerEdgeColor',[J(52,1) J(52,2) J(52,3)],'MarkerSize',8);   % 颜色渐变的等级
%            n_52=n_52+1;
%        elseif Subsidence(i)>-108-10 && Subsidence(i)<=-105-10
%            plot(x1(i),y1(i),'.','MarkerEdgeColor',[J(51,1) J(51,2) J(51,3)],'MarkerSize',8);
%            n_51=n_51+1;      
%       elseif Subsidence(i)>-105-10&& Subsidence(i)<=-102-10
%            plot(x1(i),y1(i),'.','MarkerEdgeColor',[J(50,1) J(50,2) J(50,3)],'MarkerSize',8);
%            n_50=n_50+1; 
%     elseif Subsidence(i)>-102-10 &&Subsidence(i)<=-99-10
%            plot(x1(i),y1(i),'.','MarkerEdgeColor',[J(49,1) J(49,2) J(49,3)],'MarkerSize',8);
%      











m_proj('mercator','lon',[73 140],'lat',[10 54]);
m_gshhs('hb1','linewidth',2.0,'color','k');% 陆地国界线
m_gshhs('hc1','linewidth',2.0,'color','k');
m_grid('linestyle','none','tickdir','in','linewidth',3,'fontname','Times New Roman','fontweight','bold','fontsize',16);
for i=1:194
        if SRs(i,5)>0;
            m_quiver(L1(i),B1(i),0,4*SRs(i,5)','color','r','maxheadsize',0.7,'linewidth',2.0);hold on
        elseif SRs(i,5)<0;
            m_quiver(L1(i),B1(i),0,4*SRs(i,5)','color','b','maxheadsize',0.7,'linewidth',2.0);hold on
        end
end



% 箭头红色朝上绿色朝下
% for i=1:38
%     if P1(4,i)>0;
%         quiver(L(i),B(i),0,P1(4,i)','color','r');box on;grid on; 
%     elseif P1(4,i)<0;
%         quiver(L(i),B(i),0,P1(4,i)','color','g');box on;grid on; 
%     end
% end
SRs_6=SRs_6';
for i=1:6
     subplot(3,2,i);
%      m_proj('mercator','long',[119.3 122.25],'lat',[21.75 25.75]);
%      [CS,CH]=m_etopo2('contourf',[-5000:500:0 250:250:3000],'edgecolor','none');
%      m_grid('linestyle','none','tickdir','out','linewidth',3,'fontname','Times New Roman','fontweight','bold','fontsize',16);
%      colormap([ m_colmap('blues',80); m_colmap('gland',48)]);
%      brighten(.5);colormap('bone');m_gshhs_h('color','k');caxis([-8000,4000]);
%      if i<=5;
%          m_text(121.68,25.6,['PC',num2str(i)],'fontweight','bold','fontsize',16);
%      else
%          m_text(121.68,25.6,['IC',num2str(i-5)],'fontweight','bold','fontsize',16);
%      end
     for j=1:187
         if SRs_6(i,j)>0;
             quiver(L(j),B1(j),0,6*SRs_6(i,j),'color','r','maxheadsize',0.7,'linewidth',1.5);hold on
         elseif SRs_6(i,j)<0;
             quiver(L(j),B1(j),0,6*SRs_6(i,j),'color','b','maxheadsize',0.7,'linewidth',1.5);hold on
         end
     end
%      m_quiver(119.5,22.26,0,0.65,'color','r','maxheadsize',0.7,'linewidth',2.5);
%      m_quiver(119.5,22.26,0,-0.65,'color','g','maxheadsize',0.7,'linewidth',2.5);
%      m_text(119.5,22.26,'100%','fontweight','bold','fontsize',16);
end




ha = tight_subplot(3,2,[0.01 0.03], [0.1 0.1], [0.1 0.1]);
% for i=1:6
%     axes(ha(i));
% %     if mod(i,2)==0;
% %         plot(Components(i,:),'color',[0.31,0.31,0.31]);axis([0,3652,-20,20]);
% %         set(ylabel(['IC',num2str(i/2),'/mm']),'fontweight','bold');
% %     else
% %         plot(Components(i,:),'color',[0.31,0.31,0.31]);axis([0,3652,-20,20]);
% %         set(ylabel(['PC',num2str(round(i/2)),'/mm']),'fontweight','bold');
% %     end
% end
for i=1:6
    axes(ha(i));
    plot(ICs_6(i,:),'color',[0.31,0.31,0.31]);axis([0,2374,-20,20]);
    set(ylabel(['IC',num2str(i),'/mm']),'fontweight','bold');
end
set(ha(1:6),'XTickLabel','','fontweight','bold');
set(ha(5:6),'XTick',[1,365,731,1096,1461,1827,2192]);
set(ha(5:6),'XTickLabel',{'2011','2012','2013','2014','2015','2016','2017'},'fontweight','bold');


T_offset=cell(1,45);
T_offset(cellfun('isempty',T_offset))=num2cell(0);
for i=1:45
[X_detrend(:,i) fitline(:,i)]=fitmodel(Time,X(:,i),T_offset{i});
end




RMS=RMS_GEO;
for i=1:187
    if RMS_GEO(i)>=0
        RMS(i)=RMS_GEO(i);
    elseif RMS_GEO(i)<0
         RMS(i)=NaN;
    end
end
        



















