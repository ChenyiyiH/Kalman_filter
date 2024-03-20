%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 功能说明：基于观测距离，扩展Kalman滤波完成对目标状态估计 
function [sys,x0,str,ts]=EKF(t,x,u,flag) 
global Zdist;    % 观测信息 
global Xekf;    % 粒子滤波估计状态 
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
case 3    % 计算S函数的输出，调用mdlOutputs 
        sys=mdlOutputs(t,x,u); 
    case {1,4} 
        sys=[]; 
    case 9   % 仿真结束，保存状态值 
        save('Xekf','Xekf'); 
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
global Xekf;    % 粒子滤波估计状态 
Xekf=[x0]; 
global P; 
P=zeros(2,2);  % 初始化 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2．进行离散状态变量的更新 
function sys=mdlUpdate(t,x,u,Q,R) 
global Zdist;    % 观测信息 
global Xekf; 
global P; 
Zdist=[Zdist,u]; % 保存观测信息 
%————————————————————————————————————— 
% 下面开始用EKF对状态更新 
x0=0;y0=0;  
% 第一步：状态预测 
Xold=Xekf(:,length(Xekf(1,:))); 
Xpre=ffun(Xold); 
% 第二步：观测预测 
Zpre=hfun(Xpre); 
% 第三步：求F和H 
F=[1,0;0.1*cos(0.1*Xpre(1,1)),1]; 
H=[(Xpre(1,1)-x0)/Zpre,(Xpre(2,1)-y0)/Zpre]; 
% 第四步：协方差预测 
Ppre=F*P*F'+Q; 
% 第五步：计算Kalman增益 
K=Ppre*H'*inv(H*Ppre*H'+R);
% 第六步：状态更新 
Xnew=Xpre+K*(u-Zpre); 
% 第七步：协方差更新 
P=(eye(2)-K*H)*Ppre; 
% 保存最新的状态并输出 
Xekf=[Xekf,Xnew]; 
sys=Xnew;  % 返回给输出 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3．求取系统的输出信号 
function sys=mdlOutputs(t,x,u) 
sys = x;  % 把算得的模块输出向量赋给sys 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%