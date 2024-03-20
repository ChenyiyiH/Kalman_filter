%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 功能描述：读取并显示AVI 视频程序 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function code_ReadAndShowAVI 
% 读取存放在 work 目录下的 myAVI.avi 视频文件
mov = VideoReader('C:\Users\Jhchen\OneDrive - FAU\Desktop\Kalman_Filter\MODEL_CODE\code_video\myAVI.avi');

% 读取视频的总帧数 
totalFrame = mov.NumFrames;

% 将读取的视频显示 
figure('Name','show the movie');
for k = 1:totalFrame
    currentFrame = read(mov, k);
    image(currentFrame);
    drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 对于单幅图像常用imshow函数。imshow是MATLAB中用于显示图像的函数，具
% 体的调用格式可以为
% imshow（I，n）
% imshow（I，[low high]）