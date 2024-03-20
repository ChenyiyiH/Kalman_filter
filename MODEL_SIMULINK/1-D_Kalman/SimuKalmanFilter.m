function [sys,x0,str,ts] = SimuKalmanFilter(t,x,u,flag) 
% 输入参数： 
%     t、x、u分别对应时间、状态、输入信号 
%     flag为标志位，其，取值不同，S函数执行的任务和返回的数据也是不同的 
% 输出参数： 
%     sys为一个通用的返回参数值，其数值根据flag的不同而不同 
%     x0为状态初始数值 
%     str在目前为止的MATLAB版本中并没有什么作用，一般str=[]即可 
%     ts为一个两列的矩阵，包含采样时间和偏移量两个参数 
switch flag 
    case 0  % 系统进行初始化，调用mdlInitializeSizes函数 
        [sys,x0,str,ts]=mdlInitializeSizes; 
    case 1  % 计算连续状态变量的导数，调用mdlDerivatives函数 
        sys=mdlDerivatives(t,x,u); 
    case 2  % 更新离散状态变量，调用mdlUpdate函数 
        sys=mdlUpdate(t,x,u); 
    case 3  % 计算S函数的输出，调用mdlOutputs 
        sys=mdlOutputs(t,x,u); 
    case 4  % 计算下一仿真时刻， 
        sys=mdlGetTimeOfNextVarHit(t,x,u); 
    case 9  % 仿真结束，调用mdlTerminate函数 
        sys=mdlTerminate(t,x,u); 
    otherwise   % 其他未知情况处理，用户可以自定义 
        error(['Unhandled flag = ',num2str(flag)]); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1．系统初始化子函数 
function [sys,x0,str,ts]=mdlInitializeSizes 
sizes = simsizes; 
sizes.NumContStates  = 0; 
sizes.NumDiscStates  = 1; 
sizes.NumOutputs     = 1; 
sizes.NumInputs      = 1; 
sizes.DirFeedthrough = 1; 
sizes.NumSampleTimes = 1;   % 至少需要的采样时间 
sys = simsizes(sizes); 
x0  = 5000+sqrt(5)*randn;             % 初始条件 
str = [];               % str总是设置为空 
ts  = [-1 0];  % 表示该模块采样时间继承其前的模块采样时间设置 
global P;   % 定义协方差P 
P=5;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2．进行连续状态变量的更新 
function sys=mdlDerivatives(t,x,u) 
sys = []; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3．进行离散状态变量的更新 
function sys=mdlUpdate(t,x,u) 
global P; 
F=1;  % 系统的状态方程和观测方程系数矩阵 
B=0; 
H=1;  % 系统的状态方程和观测方程系数矩阵 
Q=0;  % 过程噪声和测量噪声方差 
R=5;       
xpre=F*x+B*u; % 状态预测 
Ppre=F*P*F'+Q;% 协方差预测 
K=Ppre*H'*inv(H*Ppre*H'+R);% 计算Kalman增益 
e=u-H*xpre;   % u是输入的观测值，在此计算新息 
xnew=xpre+K*e; % 状态更新 
P=(eye(1)-K*H)*Ppre; % 协方差更新 
% 将计算的结果返回给主函数 
sys=xnew; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 4．求取系统的输出信号 
function sys=mdlOutputs(t,x,u) 
sys = x;    % 把算得的模块输出向量赋给sys 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 5．计算下一仿真时刻，由sys返回 
function sys=mdlGetTimeOfNextVarHit(t,x,u) 
sampleTime = 1;    %  此处设置下一仿真时刻为1s以后 
sys = t + sampleTime; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 6．结束仿真子函数 
function sys=mdlTerminate(t,x,u) 
sys = []; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%