%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 功能说明：Kalman滤波用在一维温度数据测量系统中 
function code_1D_TEMP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
N=120;%采样点的个数，时间单位是分钟，可理解为试验进行了60分钟的测量 
CON=25;%室内温度的理论值，在这个理论值得基础上受过程噪声会有波动 
% 对状态和测量初始化 
Xexpect=CON*ones(1,N); %期望的温度是恒定的25℃，但真实温度不可能会这样的 
X=zeros(1,N);   % 房间各时刻真实温度值 
Xkf=zeros(1,N); % Kalman滤波处理的状态，也叫估计值 
Z=zeros(1,N);   % 温度计测量值 
P=zeros(1,N);  
%赋初值 
X(1)=25.1; %假如初始值房间温度为25.1℃ 
P(1)=0.01;  %初始值的协方差 
Z(1)=24.9; 
Xkf(1)=Z(1); %初始测量值为24.9℃，可以作为滤波器的初始估计状态 
% 噪声 
Q=0.01; 
R=0.25; 
W=sqrt(Q)*randn(1,N); % 方差决定噪声的大小 
V=sqrt(R)*randn(1,N); % 方差决定噪声的大小 
% 系统矩阵 
F=1; 
G=1; 
H=1; 
I=eye(1); %本系统状态为一维 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 模拟房间温度和测量过程，并滤波 
for k=2:N 
    % 第一步：随时间推移，房间真实温度波动变化 
    % k时刻房间的真实温度，对于温度计来说，这个真实值是不知道的，但是它的 
    % 存在又是客观事实，读者要深刻领悟这个计算机模拟过程 
    X(k)=F*X(k-1)+G*W(k-1); 
     
    % 第二步：随时间推移，获取实时数据
    % 温度计对k时刻房间温度的测量，Kalman滤波是站在温度计角度进行的， 
    % 它不知道此刻真实状态X(k)，只能利用本次测量值Z(k)和上一次估计值 
    % Xkf(k)来做处理，其目标是最大限度地降低测量噪声R的影响，尽可能 
    % 地逼近X(k)，这也是Kalman滤波目的所在 
    Z(k)=H*X(k)+V(k); 
    % 第三步：Kalman滤波 
    % 有了k时刻的观测Z(k)和k-1时刻的状态，那么就可以进行滤波了, 
    % 推导公式
    % x_hat(k+1|k)=Phi*x_hat(k|k) 状态一步预测
    % x_hat(k+1|k+1)=x_hat(k+1|k)+K(k+1)epsode(k+1) 状态更新
    % epsode(k+1)=Y(k+1)-Hx_hat(k+1|k)
    % K(k+1)=P(k+1)*H^t*[H*P(k+1|k)*H^t+R]^-1  滤波增益矩阵
    % P(k+1|k)=Phi*P(k|k)*Phi^T+Gamma*Q*Gamma^T 一步预测协方差阵
    % P(k+1|k+1)=[In-K(k+1)*H]*P(k+1|k)  协方差阵更新
    % x_hat(0|0)=miu_0,P(0|0)=P0
    X_pre=F*Xkf(k-1);            % 状态预测 
    P_pre=F*P(k-1)*F'+Q;          % 协方差预测 
    Kg=P_pre*inv(H*P_pre*H'+R);  % 计算Kalman增益 
    e=Z(k)-H*X_pre;              % 新息 
    Xkf(k)=X_pre+Kg*e;           % 状态更新 
    P(k)=(I-Kg*H)*P_pre;          % 协方差更新 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 计算误差 
Err_Messure=zeros(1,N);% 测量值与真实值之间的偏差 
Err_Kalman=zeros(1,N);% Kalman估计与真实值之间的偏差 
for k=1:N 
    Err_Messure(k)=abs(Z(k)-X(k)); 
    Err_Kalman(k)=abs(Xkf(k)-X(k)); 
end 
t=1:N; 
% figure('Name','Kalman Filter Simulation','NumberTitle','off'); 
figure  % 画图显示 
% 依次输出理论值，叠加过程噪声（受波动影响）的真实值， 
% 温度计测量值，kalman估计值 
plot(t,Xexpect,'-b',t,X,'-r.',t,Z,'-ko',t,Xkf,'-g*'); 
legend('期望值','真实值','观测值','Kalman滤波值');          
xlabel('采样时间/s'); 
ylabel('温度值/℃'); 
% 误差分析图 
figure  % 画图显示 
plot(t,Err_Messure,'-b.',t,Err_Kalman,'-k*'); 
legend('测量偏差','Kalman滤波偏差');          
xlabel('采样时间/s'); 
ylabel('温度偏差值/℃'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%