%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  功能说明：扩展Kalman滤波在纯方位目标跟踪中的应用实例 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function code_EKF_angle 
clc;clear; 
T=1;          %雷达扫描周期 
N=40/T;          %总的采样次数 
X=zeros(4,N);        % 目标真实位置、速度 
X(:,1)=[0,2,1400,-10];       % 目标初始位置、速度 
Z=zeros(1,N);         % 传感器对位置的观测 
delta_w=1e-4;         % 如果增大这个参数，目标真实轨迹就是曲线了
Q=delta_w*diag([1,1]) ;       % 过程噪声均值 
G=[T^2/2,0;T,0;0,T^2/2;0,T];     % 过程噪声驱动矩阵 
R=0.1*pi/180;     %观测噪声方差,读者可以修改值观察其对角度测量的影响 
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];     % 状态转移矩阵 
x0=0;         % 观测站的位置，可以设为其他值 
y0=1000;  
Xstation=[x0;y0]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
w=sqrtm(R)*randn(1,N);            %均值为0，方差为1的高斯噪声 
for t=2:N 
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1);      %目标真实轨迹 
end 
for t=1:N 
    Z(t)=hfun(X(:,t),Xstation)+w(t);          %对目标观测 
    %  对sqrtm(R)*w(t)转化为角度sqrtm(R)*w(t)/pi*180可以看出噪声的大小（单位：度） 
end 
% EKF滤波 
Xekf=zeros(4,N); 
Xekf(:,1)=X(:,1);        % Kalman滤波状态初始化 
P0=eye(4);             % 协方差阵初始化 
for i=2:N 
    Xn=F*Xekf(:,i-1);        %预测 
    P1=F*P0*F'+G*Q*G';     %预测误差协方差 
    dd=hfun(Xn,Xstation);     % 观测预测 
    % 求雅可比矩阵H 
    D=Dist(Xn,Xstation); 
    H=[-(Xn(3,1)-y0)/D,0,(Xn(1,1)-x0)/D,0];     % 即为所求一阶近似 
    K=P1*H'*inv(H*P1*H'+R);               %增益 
    Xekf(:,i)=Xn+K*(Z(:,i)-dd);                  %状态更新 
    P0=(eye(4)-K*H)*P1;                       %滤波误差协方差更新 
end 
% 误差分析 
for i=1:N 
  Err_KalmanFilter(i)=sqrt(Dist(X(:,i),Xekf(:,i)));%滤波后误差 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 画图 
figure 
hold on;box on; 
plot(X(1,:),X(3,:),'-k.');              % 真实轨迹 
plot(Xekf(1,:),Xekf(3,:),'-r+');        % 扩展Kalman滤波轨迹 
legend('真实轨迹','EKF轨迹'); 
xlabel('横坐标X/m'); 
ylabel('纵坐标Y/m'); 
figure 
hold on; box on; 
plot(Err_KalmanFilter,'-ks','MarkerFace','r') 
xlabel('时间/s'); 
ylabel('位置估计偏差/m'); 
figure  
hold on;box on; 
plot(Z/pi*180,'-r.','MarkerFace','r');            % 真实角度值 
plot(Z/pi*180+w/pi*180,'-ko','MarkerFace','g');   % 受噪声污染的观测值 
legend('真实角度','观测角度'); 
xlabel('时间/s'); 
ylabel('角度值/°'); 
figure    % 观测噪声大小 
hold on;box on; 
plot(w,'-ko','MarkerFace','g');     % 受噪声污染的观测值 
xlabel('时间/s'); 
ylabel('观测噪声'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % 子函数 
function cita=hfun(X1,X0)        % 需要注意各个象限角度的变化 
if X1(3,1)-X0(2,1)>=0  %  y1-y0>0 
    if X1(1,1)-X0(1,1)>0 % x1-x0>0 
        cita=atan(abs( (X1(3,1)-X0(2,1))/(X1(1,1)-X0(1,1)) )); 
    elseif X1(1,1)-X0(1,1) == 0 
        cita=pi/2; 
    else 
cita=pi/2+atan(abs( (X1(3,1)-X0(2,1))/(X1(1,1)-X0(1,1)) )); 
    end 
else 
    if X1(1,1)-X0(1,1)>0 % x1-x0>0 
cita=3*pi/2+atan(abs( (X1(3,1)-X0(2,1))/(X1(1,1)-X0(1,1)) )); 
    elseif X1(1,1)-X0(1,1) == 0 
        cita=3*pi/2; 
    else 
cita=pi+atan(abs( (X1(3,1)-X0(2,1))/(X1(1,1)-X0(1,1)) )); 
    end 
end 
function d=Dist(X1,X2) 
if length(X2)<=2 
    d=( (X1(1)-X2(1))^2 + (X1(3)-X2(2))^2 ); 
else 
    d=( (X1(1)-X2(1))^2 + (X1(3)-X2(3))^2 ); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%