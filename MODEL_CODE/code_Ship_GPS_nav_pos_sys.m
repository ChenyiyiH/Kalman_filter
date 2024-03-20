%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 功能说明：Kalman滤波在船舶GPS导航定位系统中的应用 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function code_Ship_GPS_nav_pos_sys
clc;clear; 
T=1;%雷达扫描周期 
N=80/T; %总的采样次数 
X=zeros(4,N); % 目标真实位置、速度 
X(:,1)=[-100,2,200,20];% 目标初始位置、速度 
Z=zeros(2,N); % 传感器对位置的观测 
Z(:,1)=[X(1,1),X(3,1)];  % 观测初始化 
delta_w=1e-2;  %如果增大这个参数，目标真实轨迹就是曲线了 
Q=delta_w*diag([0.5,1,0.5,1]) ; % 过程噪声均值 
R=100*eye(2);  %观测噪声均值 
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];  % 状态转移矩阵 
H=[1,0,0,0;0,0,1,0];   % 观测矩阵 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=2:N
    X(:,t)=F*X(:,t-1)+sqrtm(Q)*randn(4,1);%目标真实轨迹 
    Z(:,t)=H*X(:,t)+sqrtm(R)*randn(2,1); %对目标观测 
end 
% Kalman滤波 
Xkf=zeros(4,N); 
Xkf(:,1)=X(:,1); % Kalman滤波状态初始化 
P0=eye(4); % 协方差阵初始化 
for i=2:N 
    Xn=F*Xkf(:,i-1); %预测 
    P1=F*P0*F'+Q;%预测误差协方差 
    K=P1*H'*inv(H*P1*H'+R);%增益 
    Xkf(:,i)=Xn+K*(Z(:,i)-H*Xn);%状态更新 
    P0=(eye(4)-K*H)*P1;%滤波误差协方差更新 
end 
% 误差分析 
for i=1:N 
    Err_Observation(i)=RMS(X(:,i),Z(:,i)); % 滤波前的误差 
    Err_KalmanFilter(i)=RMS(X(:,i),Xkf(:,i)); % 滤波后的误差 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 画图 
figure 
hold on;box on; 
plot(X(1,:),X(3,:),'-k'); % 真实轨迹 
plot(Z(1,:),Z(2,:),'-b.'); % 观测轨迹 
plot(Xkf(1,:),Xkf(3,:),'-r+'); % Kalman滤波轨迹 
legend('真实轨迹','观测轨迹','滤波轨迹'); 
xlabel('横坐标X/m'); 
ylabel('纵坐标Y/m'); 
figure 
hold on; box on; 
plot(Err_Observation,'-ko','MarkerFace','g') 
plot(Err_KalmanFilter,'-ks','MarkerFace','r') 
legend('滤波前误差','滤波后误差') 
xlabel('观测时间/s'); 
ylabel('误差值'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 计算欧氏距离子函数 
function dist=RMS(X1,X2); 
if length(X2)<=2 
    dist=sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(2))^2 ); 
else 
    dist=sqrt( (X1(1)-X2(1))^2 + (X1(3)-X2(3))^2 ); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%