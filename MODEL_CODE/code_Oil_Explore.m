%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 功能说明：石油地震勘测输入白噪声估值器算法仿真程序 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function code_Oil_Explore 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 参数初始化 
clear all 
T=300;  % 总时间 
F=[1,0;0.3,-0.5];  % 状态转移矩阵 
L=[-1,2]'; % 噪声矩阵 
H=[1 1];    % 观测矩阵 
R=0.1; % 观测噪声的方差 
n=2; % 状态的维数 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Bernoulli-Gaussian 白噪声生成器 
Qg=49;  % g(t)的方差 
longa=0.3; % longa的取值，b(t)的取值概率 
Q=longa*Qg;  % w(t)的方差 
randn('seed',13) 
g=sqrt(Qg)*randn(1,T+10);  % 生成g(t)
rand('state',1); 
para=rand(1,T+10);  % 产生0-1之间高斯分布的随机数值 
for t=1:T+10 
    if para(t)<longa 
        b(t)=1; 
    else 
        b(t)=0; 
    end 
    w(t)=b(t)*g(t);  % 产生w(t) 
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 状态空间模拟部分 
% 观测噪声v(t)产生 
v=sqrt(R)*randn(1,T+10); 
% 产生状态和观测信息 
X=zeros(2,T+10); 
Z=zeros(1,T+10); 
Z(1)=H*X(:,1)+v(1); 
for t=2:T+10 
    X(:,t)=F*X(:,t-1)+L*w(t-1);  % 状态方程 
    Z(t)=H*X(:,t)+v(t); % 观测方程 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Kalman滤波部分 
P0=eye(n); 
Xe=zeros(n,T+10); 
PP=[]; 
for t=1:T+8 
    XX=F*X(:,t); % 状态预测 
    % 计算协方差矩阵P 
    P=F*P0*F'+L*Q*L'; 
    PP=[PP,P]; 
    % 计算Kalman增益 
    K(:,t)=P*H'*inv(H*P*H'+R); 
    % 计算新息 
    e(:,t)=Z(t)-H*XX; 
    % 状态更新 
    Xe(:,t)=XX+K(:,t)*e(:,t); 
    % 方差更新 
    P0=(eye(n)-K(:,t)*H)*P; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 白噪声估值器部分 
N=3; % 取N步平滑 
for t=1:T+8 
    Persai(:,:,t)=F*(eye(n)-K(:,t)*H); 
    Qe(:,:,t)=H*PP(:,2*(t-1)+1:2*t)*H'+R; 
end 
for t=1:T+5 
    M(1,t)=Q*L'*H'*inv(Qe(:,:,t+1)); 
    M(2,t)=Q*L'*Persai(:,:,t+1)'*H'*inv(Qe(:,:,t+2)); 
    M(3,t)=Q*L'*Persai(:,:,t+2)'*Persai(:,:,t+1)'*H'*inv(Qe(:,:,t+3)); 
end 
for t=1:T 
    wjian(1,t)=M(1,t+1)*e(t+1); % 一步平滑器 
    wjian(2,t)=wjian(1,t)+M(2,t+2)*e(t+2);  % 二步平滑 
    wjian(3,t)=wjian(2,t)+M(3,t+3)*e(t+3);  % 三步平滑 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  画图部分 
for Num=1:N 
    subplot(3,1,Num); 
    t=1:T; 
    plot(t,wjian(Num,t),'b.'); 
    for t=1:T 
        hh=line( [t,t],[0,w(t)] ); 
        set(hh,'color','k'); 
    end 
    xlabel(['w(t)和',num2str(Num),'步平滑器']) 
    ylabel('w的数值') 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%