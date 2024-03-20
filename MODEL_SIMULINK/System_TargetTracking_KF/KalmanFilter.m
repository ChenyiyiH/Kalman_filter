% 功能说明：S函数仿真系统的状态方程X(k+1)=F*x(k)+G*w(k) 
function [sys,x0,str,ts]=KalmanFilter(t,x,u,flag) 
global Xkf; 
global Z; 
switch flag 
    case 0  % 系统进行初始化，调用mdlInitializeSizes函数 
        [sys,x0,str,ts]=mdlInitializeSizes; 
    case 2  % 更新离散状态变量，调用mdlUpdate函数 
        sys=mdlUpdate(t,x,u); 
case 3    % 计算S函数的输出，调用mdlOutputs 
        sys=mdlOutputs(t,x,u); 
    case {1,4} 
        sys=[]; 
    case 9 
        save('Zobserv','Z'); 
        save('Xkalman','Xkf'); 
    otherwise   % 其他未知情况处理，用户可以自定义 
    error(['Unhandled flag = ',num2str(flag)]); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1．系统初始化子函数 
function [sys,x0,str,ts]=mdlInitializeSizes 
global P; 
global Xkf; 
global Z;   %  
P=0.1*eye(4); 
sizes = simsizes; 
sizes.NumContStates  = 0;   % 无连续量 
sizes.NumDiscStates  = 4;   % 离散状态4维 
sizes.NumOutputs     = 4;   % 输出4维，应为状态量是x-y方向的位置和速度 
sizes.NumInputs      = 2;   % 输入维数，因为噪声模型是2维的 
sizes.DirFeedthrough = 0; 
sizes.NumSampleTimes = 1;   % 至少需要的采样时间 
sys = simsizes(sizes); 
x0  = [10,5,12,5]';             % 初始条件 
str = [];               % str总是设置为空 
ts  = [-1 0];  % 表示该模块采样时间继承其前的模块采样时间设置 
Xkf=[]; 
Z=[]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2．进行离散状态变量的更新 
function sys=mdlUpdate(t,x,u) 
global P; 
global Xkf; 
global Z;   
F=[1,1,0,0;0,1,0,0;0,0,1,1;0,0,0,1]; % 状态转移矩阵 
G=[0.5,0;1,0;0,0.5;0,1];             % 过程噪声驱动矩阵 
H=[1,0,0,0;0,0,1,0];                 % 观测矩阵 
Q=diag([0.0001,0.0009]);                 % 过程噪声方差 
R=diag([0.04,0.04]);                 % 测量噪声方差 
% Kalman滤波 
Xpre=F*x;   % 状态预测 
Ppre=F*P*F'+G*Q*G';% 协方差更新 
K=Ppre*H'*inv(H*Ppre*H'+R); % 计算Kalman增益 
e=u-H*Xpre; % 计算新息 
Xnew=Xpre+K*e; % 状态更新 
P=(eye(4)-K*H)*Ppre; % 协方差更新 
sys=Xnew;   % 将计算的结果返回给主函数 
% 保存观测值和滤波结果 
Z=[Z,u]; 
Xkf=[Xkf,Xnew]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3．求取系统的输出信号 
function sys=mdlOutputs(t,x,u) 
sys = x;    % 把算得的模块输出向量赋给sys 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%