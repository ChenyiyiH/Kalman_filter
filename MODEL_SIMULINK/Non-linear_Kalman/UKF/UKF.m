%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 功能说明：基于观测距离，无迹Kalman滤波完成对目标状态估计 
function [sys,x0,str,ts]=UKF(t,x,u,flag) 
global Zdist;    % 观测信息 
global Xukf;    % 粒子滤波估计状态 
%randn('seed',20); 
% 过程噪声Q 
Q=diag([0.01,0.04]);  
% 测量噪声R 
R=1; 
switch flag 
    case 0  % 系统进行初始化，调用mdlInitializeSizes函数 
        [sys,x0,str,ts]=mdlInitializeSizes; 
    case 2  % 更新离散状态变量，调用mdlUpdate函数 
sys=mdlUpdate(t,x,u,Q,R); 
    case 3  % 计算S函数的输出，调用mdlOutputs 
        sys=mdlOutputs(t,x,u); 
    case {1,4} 
        sys=[]; 
    case 9   % 仿真结束，保存状态值 
        save('Xukf','Xukf'); 
        save('Zdist','Zdist'); 
    otherwise   % 其他未知情况处理，用户可以自定义 
        error(['Unhandled flag = ',num2str(flag)]); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1．系统初始化子函数 
function [sys,x0,str,ts]=mdlInitializeSizes(N) 
sizes = simsizes; 
sizes.NumContStates  = 0;   % 无连续量 
sizes.NumDiscStates  = 2;   % 离散状态4维 
sizes.NumOutputs     = 2;   % 输出4维，应为状态量是x-y方向的位置和速度 
sizes.NumInputs      = 1;   % 输入维数，因为噪声模型是2维的 
sizes.DirFeedthrough = 0; 
sizes.NumSampleTimes = 1;   % 至少需要的采样时间 
sys = simsizes(sizes); 
x0  = [0,0]';             % 初始条件 
str = [];               % str总是设置为空 
ts  = [-1 0];  % 表示该模块采样时间继承其前的模块采样时间设置 
global Zdist;    % 观测信息 
Zdist=[]; 
global Xukf;    % 粒子滤波估计状态 
Xukf=[x0]; 
global P; 
P=0.01*eye(2);  % 初始化 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2．进行离散状态变量的更新 
function sys=mdlUpdate(t,x,u,Q,R) 
global Zdist;    % 观测信息 
global Xukf; 
global P; 
Zdist=[Zdist,u]; % 保存观测信息 
%————————————————————————————————————— 
% 下面开始用UKF对状态更新 
Xin=Xukf(:,length(Xukf(1,:))); 
[Xnew,P]=GetUkfResult(Xin,u,P,Q,R)  
% 保存最新的状态并输出 
Xukf=[Xukf,Xnew]; 
sys=Xnew;  % 返回给输出 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3．求取系统的输出信号 
function sys=mdlOutputs(t,x,u) 
sys = x;    % 把算得的模块输出向量赋给sys 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 子函数 
function [Xout,P]=GetUkfResult(Xin,Z,P,Q,R) % UKF算法子程序 
% 权值系数初始化 
L=2; % 状态的维数 
alpha=0.01; 
kalpha=0; 
belta=2; 
ramda=alpha^2*(L+kalpha)-L; 
for j=1:2*L+1 
    Wm(j)=1/(2*(L+ramda)); 
    Wc(j)=1/(2*(L+ramda)); 
end 
Wm(1)=ramda/(L+ramda); 
Wc(1)=ramda/(L+ramda)+1-alpha^2+belta; 
% 滤波初始化 
xestimate=Xin; 
% 第一步：获得一组sigma点集 
P 
cho=(chol(P*(L+ramda)))'; 
for k=1:L 
    xgamaP1(:,k)=xestimate+cho(:,k); 
    xgamaP2(:,k)=xestimate-cho(:,k); 
end 
Xsigma=[xestimate,xgamaP1,xgamaP2]; %Sigma点集 
% 第二步：对Sigma点集进行一步预测 
for k=1:2*L+1 
    Xsigmapre(:,k)=ffun(Xsigma(:,k)); 
end 
%第三步：利用第二步的结果计算均值和协方差 
Xpred=zeros(2,1);    % 均值 
for k=1:2*L+1 
    Xpred=Xpred+Wm(k)*Xsigmapre(:,k); 
end 
Ppred=zeros(2,2);    % 协方差阵预测 
for k=1:2*L+1 
    Ppred=Ppred+Wc(k)*(Xsigmapre(:,k)-Xpred)*(Xsigmapre(:,k)-Xpred)'; 
end 
Ppred=Ppred+Q; 
% 第四步：根据预测值，再一次使用UT变换，得到新的sigma点集 
chor=(chol((L+ramda)*Ppred))'; 
for k=1:L 
    XaugsigmaP1(:,k)=Xpred+chor(:,k); 
    XaugsigmaP2(:,k)=Xpred-chor(:,k); 
end 
Xaugsigma=[Xpred XaugsigmaP1 XaugsigmaP2]; 
% 第五步：观测预测 
for k=1:2*L+1    %  观测预测 
    Zsigmapre(1,k)=hfun(Xaugsigma(:,k)); 
end 
% 第六步：计算观测预测均值和协方差 
Zpred=0;         %  观测预测的均值 
for k=1:2*L+1 
    Zpred=Zpred+Wm(k)*Zsigmapre(1,k); 
end 
Pzz=0; 
for k=1:2*L+1 
    Pzz=Pzz+Wc(k)*(Zsigmapre(1,k)-Zpred)*(Zsigmapre(1,k)-Zpred)'; 
end 
Pzz=Pzz+R;  % 得到协方差Pzz 
Pxz=zeros(2,1); 
for k=1:2*L+1 
    Pxz=Pxz+Wc(k)*(Xaugsigma(:,k)-Xpred)*(Zsigmapre(1,k)-Zpred)'; 
end 
% 第七步：计算Kalman增益 
K=Pxz*inv(Pzz);                   % Kalman增益 
%第八步：状态和方差更新 
Xout=Xpred+K*(Z-Zpred);           % 状态更新 
P=Ppred-K*Pzz*K';                 % 方差更新 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%