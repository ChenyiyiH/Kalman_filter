%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 跟踪偏差分析 
function DeviationAnalysis 
load Xstate;
load Xkalman;
load Zobserv; 
% Xstate % 在命令窗口查看真实状态值 
% Xkf    % 在命令窗口查看Kalman滤波结果 
% Z      % 在命令窗口查看观测结果 
% 计算误差 
T1=length(Xstate(1,:));  % 虽然仿真时间设100，但是Xstate的列数却为102 
T2=length(Z(1,:));       % 虽然仿真时间设100，但是Z的列数却为101 
T=min(T1,T2); 
Div_Observ_Real=zeros(1,T); 
Div_Kalman_Real=zeros(1,T); 
for i=2:T 
    Div_Observ_Real(i)=sqrt((Z(1,i)-Xstate(1,i))^2+(Z(2,i)-Xstate(3,i))^2 ); 
    Div_Kalman_Real(i)=sqrt((Xkf(1,i)-Xstate(1,i))^2+(Xkf(3,i)-Xstate(3,i))^2 ); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 画轨迹图 
figure   
hold on;box on; 
plot(Xstate(1,:),Xstate(3,:),'-r.'); 
plot(Xkf(1,:),Xkf(3,:),'-k+'); 
plot(Z(1,:),Z(2,:),'-k*'); 
% 偏差图 
figure 
hold on;box on; 
plot(Div_Observ_Real,'-ko','MarkerFace','g'); 
plot(Div_Kalman_Real,'-ks','MarkerFace','b'); 
legend('Observ','Kalman');