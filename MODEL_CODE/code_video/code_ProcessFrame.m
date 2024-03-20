%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 功能描述：读取AVI，并将AVI 视频的每一帧转为bmp图片存储 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function code_ProcessFrame
% 读取存放在 work 目录下的 myAVI.avi 视频文件
mov = VideoReader('C:\Users\Jhchen\OneDrive - FAU\Desktop\Kalman_Filter\MODEL_CODE\code_video\myAVI.avi');

% 读取关于该视频的一些参数，比如总的帧数 
totalFrame = mov.NumFrames; % 读取视频的总帧数

% 将每一帧转为 BMP 图片序列，在 work 文件夹下创建空文件夹 imageFrame
outputFolder = 'C:\Users\Jhchen\OneDrive - FAU\Desktop\Kalman_Filter\MODEL_CODE\code_video\imageFrame'; % 保存图片的文件夹路径
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% 把每一帧保存为 BMP 图片
for i = 1:totalFrame 
    % 获取视频帧 
    frameData = read(mov, i);
    % 对每一帧序列命名并且编号 
    bmpName = fullfile(outputFolder, ['image', num2str(i), '.bmp']); 
    % 把每帧图像存入硬盘  
    imwrite(frameData, bmpName, 'BMP'); 
    pause(0.02); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%