function [B,Z]=mFastICA(X,num)
%并行运算同时估计多个独立成分
%X每行为混合信号
%num为要估计的个数
%%
%-----------去均值---------
%m为个数
 
% [M,T] = size(X); %获取输入矩阵的行/列数，行数为观测数据的数目，列数为采样点数      
%  average= mean(X')';  %均值
%  for i=1:M
%      X(i,:)=X(i,:)-average(i)*ones(1,T); 
%  end

%PCA降维
P = pca(X,num,'svd');
X=P*X;


%---------白化/球化------
Cx = cov(X',1);    %计算协方差矩阵Cx
[eigvector,eigvalue] = eig(Cx); %计算Cx的特征值和特征向量
Q=eigvalue^(-1/2)*eigvector';   %白化矩阵
Z=Q*X;   %正交矩阵

%

 
%----------迭代-------
Maxcount=100000;        %最大迭代次数
Critical=0.00001;   %判断是否收敛
m=num;                %需要估计的分量的个数
W=rand(m);
%对称正交化（独立成分分析P138）
%W=(W*W').^(-1/2)*W;
[E1,D1] = eig(W*W');
W=E1*diag((diag(D1)).^(-1/2))*E1'*W;

LastW=zeros(m,m);
count=0;

while abs(W-LastW)&abs(W+LastW)>Critical
%while max(max(abs(W-LastW)))>Critical
    LastW=W;
    for n=1:m
 %       WP=W(:,n);  %初始权矢量（任意）
        LastWP=LastW(:,n);
        for i=1:m     
            WP(i)=mean(Z(i,:).*(tanh(LastWP'*Z)))-mean(1-(tanh((LastWP)'*Z)).^2)*LastWP(i);
        end
        W(:,n)=WP;
    end    
    [E1,D1] = eig(W*W');
    W=E1*diag((diag(D1)).^(-1/2))*E1'*W;
    count=count+1;
     str=(['第' num2str(count) '次迭代']);
      disp(str);
end
    
Z=W'*Z;

B=W'*Q*P;
% B=W'*P;
end


% 
% 
% for n=1:m  
%     WP=W(:,n);  %初始权矢量（任意）
% %     Y=WP'*Z;
% %     G=Y.^3;%G为非线性函数，可取y^3等
% %     GG=3*Y.^2;  %G的导数
%     count=0;
%     LastWP=zeros(m,1);
%     W(:,n)=W(:,n)/norm(W(:,n));
%     
%     while abs(WP-LastWP)&abs(WP+LastWP)>Critical
%         count=count+1;   %迭代次数
%         LastWP=WP;      %上次迭代的值
%        % WP=1/T*Z*((LastWP'*Z).^3)'-3*LastWP;
%         for i=1:m     
%             WP(i)=mean(Z(i,:).*(tanh((LastWP)'*Z)))-(mean(1-(tanh((LastWP))'*Z).^2)).*LastWP(i);
%         end
%         WPP=zeros(m,1);
%         for j=1:n-1
%             WPP=WPP+(WP'*W(:,j))*W(:,j);
%         end
%         WP=WP-WPP;
%         WP=WP/(norm(WP));
%                
%         
%         if count==Maxcount
%             fprintf('未找到相应的信号'); 
%             return; 
%         end
%     end
%     W(:,n)=WP;
%     
%        
%     
% end