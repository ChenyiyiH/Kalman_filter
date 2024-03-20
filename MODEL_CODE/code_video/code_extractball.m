%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 文件名：extractball.m 
% 功能说明： 
%          目标提取子函数，提取目标区域的中心和能包含目标的最大半径 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [cc,cr,radius,flag]=code_extractball(Imwork,Imback,index) 
% 初始化目标区域中心的坐标、半径 
cc = 0;      % 目标中心x坐标 
cr = 0;       % 目标中心y坐标 
radius = 0;    % 目标区域半径 
flag = 0;        % 检测到目标标志 
[MR,MC,Dim] = size(Imback); 
% 除去背景，找到最大的不同区域（也即目标区域） 
fore = zeros(MR,MC); 
%  背景相减，得到目标 
fore = (abs(Imwork(:,:,1)-Imback(:,:,1)) > 10) ... 
    | (abs(Imwork(:,:,2) - Imback(:,:,2)) > 10) ... 
    | (abs(Imwork(:,:,3) - Imback(:,:,3)) > 10); 
% 图像腐蚀，除去微小的白噪声点 
% bwmorph该函数的功能是：提取二进制图像的轮廓 
foremm = bwmorph(fore,'erode',2);    %  “2”为次数 
% 选取最大的目标 
labeled = bwlabel(foremm,4);  % 黑背景中甄别有多少白块块，4-联通（上下左右） 
% labeled 是标记矩阵，图像分割后对不同的区域进行不同的标记 
stats = regionprops(labeled,['basic']);  %basic mohem nist 
[N,W] = size(stats); 
if N < 1 
    return    % 一个目标区域也没检测到就返回 
end 
% 在N个区域中，冒泡算法（从大到小）排序 
id = zeros(N); 
for i = 1 : N 
    id(i) = i; 
end 
% 将检测到的目标区域按照从大到小排列（冒泡排序） 
for i = 1 : N-1 
        for j = i+1 : N 
        if stats(i).Area < stats(j).Area 
            % 冒泡算法程序 
            tmp = stats(i); 
            stats(i) = stats(j); 
            stats(j) = tmp; 
            tmp = id(i); 
            id(i) = id(j); 
            id(j) = tmp; 
        end
    end 
end 
% 确保至少有一个较大的区域（具体如下，最大区域面积要大于100） 
if stats(1).Area < 100 
    return  % 如果第一个（也是最大的目标区域）都不符合，则没有检测到目标，返回
end 
selected = (labeled == id(1)); 
% 计算最大区域的中心和半径 
centroid = stats(1).Centroid; 
radius = sqrt(stats(1).Area/pi); 
cc = centroid(1); % 相当于x,中心坐标位置 
cr = centroid(2); % 相当于y 
flag = 1; % 检测到目标，将标志设置为1 
return 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%