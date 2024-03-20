%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  功能说明： UKF在目标跟踪中的应用 
%  参数说明： 1．状态6维，x方向的位置、速度、加速度； 
%                y方向的位置、速度、加速度； 
%             2．观测信息为距离和角度； 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function code_ukf_for_track_6_div_system 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=6;                                                                        % 状态位数 
t=0.5;                                                                      % 采样时间 
Q=[1 0 0 0 0 0; 
   0 1 0 0 0 0; 
   0 0 0.01 0 0 0; 
   0 0 0 0.01 0 0; 
   0 0 0 0 0.0001 0; 
   0 0 0 0 0 0.0001];                                                       %过程噪声协方差阵 
R = [100 0; 
    0 0.001^2];                                                             %量测噪声协方差阵  
% 状态方程 
f=@(x)[x(1)+t*x(3)+0.5*t^2*x(5);x(2)+t*x(4)+0.5*t^2*x(6);... 
    x(3)+t*x(5);x(4)+t*x(6);x(5);x(6)]; 
%x1为X轴位置，x2为Y轴位置，x3、x4分别是X、 
%Y轴的速度，x5、x6为X、Y两方向的加速度 
% 观测方程 
h=@(x)[sqrt(x(1)^2+x(2)^2);atan(x(2)/x(1))]; 
s=[1000;5000;10;50;2; -4]; 
x0=s+sqrtm(Q)*randn(n,1);                                                   % 初始化状态 
P0 =[100 0 0 0 0 0; 
        0 100 0 0 0 0; 
        0 0 1 0 0 0; 
        0 0 0 1 0 0; 
        0 0 0 0 0.1 0; 
        0 0 0 0 0 0.1];                                                     % 初始化协方差 
N=50;                                                                       % 总仿真时间步数，即总时间 
Xukf = zeros(n,N);                                                          % UKF滤波状态初始化 
X = zeros(n,N);                                                             % 真实状态 
Z = zeros(2,N);                                                             % 测量值 
for i=1:N 
    X(:,i)= f(s)+sqrtm(Q)*randn(6,1);                                       % 模拟，产生目标运动真实轨迹 
    s = X(:,i); 
end 
ux=x0;                                                                      % ux为中间变量 
for k=1:N 
    Z(:,k)= h(X(:,k)) + sqrtm(R)*randn(2,1);                                % 测量值
    % 保存观测 
    [Xukf(:,k),P0] = ukf(f,ux,P0,h,Z(:,k),Q,R);                             % 调用ukf滤波算法 
    ux=Xukf(:,k); 
end 
% 跟踪误差分析 
% 这里只分析位置误差，速度、加速度误差分析在此略 
for k=1:N 
    RMS(k)=sqrt( (X(1,k) -Xukf(1,k))^2+(X(2,k) -Xukf(2,k))^2 ); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 画图，轨迹图 
figure 
t=1:N; 
hold on;box on; 
plot( X(1,t),X(2,t), 'k-') 
plot(Z(1,t).*cos(Z(2,t)),Z(1,t).*sin(Z(2,t)),'-b.') 
plot(Xukf(1,t),Xukf(2,t),'-r.') 
legend('实际值','测量值','UKF估计值'); 
xlabel('x方向位置/m') 
ylabel('y方向位置/m') 
% 误差分析图
figure 
box on; 
plot(RMS,'-ko','MarkerFace','r') 
xlabel('t/s') 
ylabel('偏差/m') 
%title('跟踪位置偏差') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% UKF子函数 
function [X,P]=ukf(ffun,X,P,hfun,Z,Q,R) 
% 非线性系统中UKF算法 
L=numel(X);                              % 状态维数 
m=numel(Z);                              % 观测维数 
alpha=1e-2;                              % 默认系数，参看UT变换，下同 
ki=0;                                    % 默认系数 
beta=2;                                  % 默认系数 
lambda=alpha^2*(L+ki) -L;                % 默认系数 
c=L+lambda;                              % 默认系数 
Wm=[lambda/c 0.5/c+zeros(1,2*L)];        % 权值 
Wc=Wm; 
Wc(1)=Wc(1)+(1-alpha^2+beta);            % 权值 
c=sqrt(c); 
 
% 第一步：获得一组Sigma点集 
% Sigma点集，在状态X附近的点集，X是6*13矩阵，每列为1样本 
Xsigmaset=sigmas(X,P,c);  
 
% 第二、三、四步：对Sigma点集进行一步预测，得到均值X1means和方差P1和新sigmaX1 
% 对状态UT变换 
[X1means,X1,P1,X2]=ut(ffun,Xsigmaset,Wm,Wc,L,Q);    
 
% 第五、六步：得到观测预测，Z1为X1集合的预测，Zpre为Z1的均值， 
% Pzz为协方差， 
[Zpre,Z1,Pzz,Z2]=ut(hfun,X1,Wm,Wc,m,R);    % 对观测UT变换 
Pxz=X2*diag(Wc)*Z2';                     % 协方差Pxz 
 
% 第七步：计算Kalman增益 
K=Pxz*inv(Pzz); 
 
% 第八步：状态和方差更新 
X=X1means+K*(Z-Zpre);                   % 状态更新 
P=P1-K*Pxz';                             % 协方差更新 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% UT变换子函数 
% 输入：fun为函数句柄，Xsigma为样本集，Wm和Wc为权值， 
%       n为状态维数(n=6)，COV为方差  
% 输出：Xmeans为均值， 
function [Xmeans,Xsigma_pre,P,Xdiv]=ut(fun,Xsigma,Wm,Wc,n,COV) 
LL=size(Xsigma,2);% 得到Xsigma样本个数 
Xmeans=zeros(n,1);%均值 
Xsigma_pre=zeros(n,LL); 
for k=1:LL 
    Xsigma_pre(:,k)=fun(Xsigma(:,k));  % 一步预测 
    Xmeans=Xmeans+Wm(k)*Xsigma_pre(:,k); 
end 
% Xmeans(:,ones(1,LL))将Xmeans扩展成n*LL矩阵，每一列都相等 
Xdiv=Xsigma_pre-Xmeans(:,ones(1,LL));   % 预测减去均值 
P=Xdiv*diag(Wc)*Xdiv'+COV;    % 协方差 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 产生Sigma点集函数 
function Xset=sigmas(X,P,c) 
A = c*chol(P)';% Cholesky分解 
Y = X(:,ones(1,numel(X))); 
Xset = [X Y+A Y-A]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%