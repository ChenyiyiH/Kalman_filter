%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 程序说明：这是一个视频捕获并录制的程序
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function code_VideoCapture
% 生成一个 AVI 对象，并取名字 myAVI
aviobj = VideoWriter('myAVI.avi');
open(aviobj); % 打开 AVI 文件以准备写入

% 创建一个视频输入对象
obj = videoinput('winvideo',1,'YUY2_320x240');
% 预览视频
preview(obj);

% 设置视频录制的帧数（时间），以便让录制停止
T = 100;
k = 0;
while k < T
    % 捕获图像并显示在左边小窗口中
    frame = getsnapshot(obj);
    subplot(1,2,1)
    imshow(frame);
    
    % 转成彩色,显示在右边，这个 frame 就可以按照自己意愿处理了
    frameRGB = ycbcr2rgb(frame);
    subplot(1,2,2);
    imshow(frameRGB);
    drawnow;
    
    % 将每一帧图像保存在创建的 AVI 视频里面
    writeVideo(aviobj, frameRGB);
    
    % 序列自增
    k = k + 1;
end

% 关闭 AVI 对象，同时将截取的图像统一写入 AVI 文件中
close(aviobj);
% 删除对象
delete(obj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%