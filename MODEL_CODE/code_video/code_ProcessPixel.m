%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 功能描述：操作视频帧中的像素，在每一帧中打上标签 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function code_ProcessPixel
% 读取存放在 work 目录下的 myAVI.avi 视频文件
mov = VideoReader('C:\Users\Jhchen\OneDrive - FAU\Desktop\Kalman_Filter\MODEL_CODE\code_video\video.avi');

% 读取关于该视频的一些参数，比如总的帧数 
totalFrame = mov.NumFrames; % 读取视频的总帧数

% 读取 Logo 图片
imageLogo = imread('logo.png');
imageSize = imresize(imageLogo, 0.2); % 缩小 Logo 图片为原来的一半

[height, width, ~] = size(imageSize); % 获取 Logo 图片的尺寸
figure('Name','Processing Pixel'); 

for i = 1:totalFrame
    % 获取视频帧
    frameData = read(mov, i);
    
    % 将 Logo 图片添加到视频帧的左上角
    frameData(1:height, 1:width, :) = imageSize;
    
    % 显示原视频和处理后的视频
    subplot(1,2,1);
    imshow(read(mov, i));
    xlabel('The original video');
    
    subplot(1,2,2);
    imshow(frameData);
    xlabel('The processed video');

    % 需要暂停一下，否则画面太快，会感觉是静止的 
    pause(0.02);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%