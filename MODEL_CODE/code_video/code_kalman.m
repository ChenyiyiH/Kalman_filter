%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 文件名：kalman.m 
% 功能说明： Kalman滤波算法实现视频目标位置跟踪 
%            主要滤除跟踪过程的观测噪声 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function code_kalman 
clear,clc 
% 计算背景图像 
Imzero = zeros(240,320,3); 
for i = 1:5 
    Im{i} = double(imread(['DATA/',int2str(i),'.jpg'])); 
    Imzero = Im{i}+Imzero; 
end 
Imback = Imzero/5; 
[MR,MC,Dim] = size(Imback); 
% Kalman 滤波器初始化 
R=[[0.2845,0.0045]',[0.0045,0.0455]'];   % 观测噪声 
H=[1 0 0 0;0 1 0 0];                    % 观测矩阵 
Q=0.01*eye(4);                       % 过程噪声方差初始化 
P = 100*eye(4);                        % 协方差初始化 
dt=1;                                 % 采样时间间隔 
A=[[1,0,0,0]',[0,1,0,0]',[dt,0,1,0]',[0,dt,0,1]'];    % 状态转移矩阵 
g = 6;               % pixels^2/time step 
Bu = [0,0,0,g]';     % 过程噪声驱动矩阵，也即加速度 
kfinit=0;            % Kalman滤波器初始化标志位 
x=zeros(100,4); 
% 遍历所有图像帧 
for i = 1 : 60 
    % 读取图像帧 
    Im = (imread(['DATA/',int2str(i), '.jpg'])); 
    imshow(Im) 
    imshow(Im) 
    Imwork = double(Im); 
    % 调用目标检测函数，提取视频中的球 
    [cc(i),cr(i),radius,flag] = code_extractball(Imwork,Imback,i); 
    if flag == 0    % 没检测到球，跳到下一帧图像 
        continue 
    end 
    hold on 
        for c = -1*radius: radius/20 : 1*radius 
        r = sqrt(radius^2-c^2); 
        plot(cc(i)+c,cr(i)+r,'g.') 
        plot(cc(i)+c,cr(i)-r,'g.') 
end 
    % Kalman 滤波算法 
    i       % 显示图像帧的进度 
if kfinit == 0     % 如果Kalman滤波标志位为0，说明状态没有初始化 
    xp = [MC/2,MR/2,0,0]' ;   % 给出初始状态 
else 
        xp=A*x(i-1,:)' + Bu;     % 状态预测 
end 
kfinit=1;      % 状态初始化以后，将标志位置成1 
        PP = A*P*A' + Q ;     % 协方差预测 
    K = PP*H'*inv(H*PP*H'+R);   % Kalman增益 
        x(i,:) = (xp + K*([cc(i),cr(i)]' - H*xp))';   % 目标状态更新 
    P = (eye(4)-K*H)*PP;     % 协方差更新 
    hold on 
        for c = -1*radius: radius/20 : 1*radius 
        r = sqrt(radius^2-c^2); 
        % 红色的为Kalman滤波器估计的位置 
        % 以目标区域的中心及半径radius画圆，表示跟踪区域 
        plot(x(i,1)+c,x(i,2)+r,'r.')  ; 
        plot(x(i,1)+c,x(i,2)-r,'.');
    end 
    pause(0.3); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 画图，对比观测值和滤波器估计值的位置比较 
figure 
% x方向 
hold on; box on; 
plot(cc,'-r*') 
plot(x(:,1),'-b.') 
figure 
% y方向 
hold on; box on; 
plot(cr,'-g*') 
plot(x(:,2),'-b.') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 从球的状态中，估计图像的观测噪声方差R 
posn = [cc(55:60)',cr(55:60)']; 
mp = mean(posn); 
diffp = posn - ones(6,1)*mp; 
Rnew = (diffp'*diffp)/5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%