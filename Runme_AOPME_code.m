%固定O点,新算法
clc;clear;close all
inFile = fullfile('Case3_Simulation.mp4');
vr = VideoReader(inFile);%打开视频
[~, writeTag, ~] = fileparts(inFile);%定义结果输出文件名称
FrameRate = vr.FrameRate; %读帧率
vid = vr.read(); %读视频
nF = vr.NumFrames; %读帧数

%% 预定义:帧数nF 一次循环的计算帧数frames 像素与相位的比例系数pixels
%% 参与计算的坐标点个数dgs
frames=2;
pixels=7.2;
steps=1;
ym=80;
zhen=nF-2;
yz=1;
t=0.6157;
size=5;
for n=1:nF
    I1=vid(:,:,:,n);
    I1=rgb2gray(I1);
    I(:,:,n)=im2double(I1);
end
tic

% 1. 交互式初始化
figure('Name', '请设置追踪参数');
imshow(I(:,:,1));
title('1.Click Center  2.Click Outer Edge');
[x_pts, y_pts] = ginput(2);
xc = x_pts(1); yc = y_pts(1);
close;

zswzx(1,1)=[x_pts(2)];
zswzy(1,1)=[y_pts(2)];
origin=[xc yc];
R=abs(sqrt((zswzx(1,1)-origin(1,1))^2+(zswzy(1,1)-origin(1,2))^2));
%% 角度+象限判断
dx=zswzx(1,1)-origin(1,1);
dy=origin(1,2)-zswzy(1,1);
gabor_rad(1,1)=atan2(dx,dy);
if gabor_rad(1,1)>=0
    gabor_rad(1,1)=gabor_rad(1,1);
else
    gabor_rad(1,1)=gabor_rad(1,1)+pi;
end
%% 卷积模板
x = 0;
for i = linspace(-ym/2,ym/2,ym)
    x = x + 1;
    y = 0;
    for j = linspace(-ym/2,ym/2,ym)
        y = y + 1;
        z(y,x)=compute1(i,j,gabor_rad(1,1));
    end
end

%% 卷积过程 先卷第一批次
filtered=cell(1,nF);
fai1=cell(1,nF);
for jj=1:frames
    filtered{jj}=conv2(I(:,:,jj),z,'same');
    fai1{jj}=angle(filtered{jj});
end

%%i的取值连续两组不需要取到 因为i=i+1;
for e=1:zhen %zhen=fix求商 45/4=11 即为e=1:11
    f1=fai1{e}(round(zswzy(e,1)),round(zswzx(e,1)));
    f2=fai1{e+1}(round(zswzy(e,1)),round(zswzx(e,1)));
      if f1>=0
          if f2<0 && f1-f2>=pi
            d(e,1)=f2-f1+2*pi;
          else
            d(e,1)=f2-f1;
          end
      elseif f1<0
          if f2>0 && f2-f1>=pi
            d(e,1)=f2-f1-2*pi;
          else
            d(e,1)=f2-f1;
          end
       end
wzx(e+1,1)=zswzx(e,1)-d(e,1)*cos(gabor_rad(e,1))*pixels;
wzy(e+1,1)=zswzy(e,1)-d(e,1)*sin(gabor_rad(e,1))*pixels;
  
    if mod(e,steps)==0  
        gray_image=I(:,:,e+1);
        rwzx=round(wzx(e+1,1));rwzy=round(wzy(e+1,1));
        LY=[gray_image(rwzy-1,rwzx);gray_image(rwzy+1,rwzx);gray_image(rwzy,rwzx-1);gray_image(rwzy,rwzx+1)];
        tt=0;
        if any(LY)>t&&any(LY)<t
            tt==1;
        else
            tt==0;
        end
        if tt==0
            % 搜索所有可能的点
             minDist=inf;
        for it = rwzy-size: rwzy+size
            for jt = rwzx-size : rwzx+size
                if ~(it == rwzy && jt == rwzx) % 避免检查原始点
                    % 计算四邻域
                    LY = [gray_image(it-1, jt);
                          gray_image(it+1, jt);
                          gray_image(it, jt-1);
                          gray_image(it, jt+1)];
                    % 判断是否满足条件
                    if any(LY > t) && any(LY < t)
                        dist = sqrt((it - rwzy)^2 + (jt - rwzx)^2);
                        if dist < minDist
                            minDist = dist;
                            nearestPoint = [it, jt];
                        end
                    end
                end
            end
        end
       %%%%%%
        if minDist==Inf;
            size=15;
 for it = rwzy-size: rwzy+size
            for jt = rwzx-size : rwzx+size
                if ~(it == rwzy && jt == rwzx) % 避免检查原始点
                    % 计算四邻域
                    LY = [gray_image(it-1, jt);
                          gray_image(it+1, jt);
                          gray_image(it, jt-1);
                          gray_image(it, jt+1)];
                    % 判断是否满足条件
                    if any(LY > t) && any(LY < t)
                        dist = sqrt((it - rwzy)^2 + (jt - rwzx)^2);
                        if dist < minDist
                            minDist = dist;
                            nearestPoint = [it, jt];
                        end
                    end
                end
            end
            end
        end
%%%%%%%%%%
          swzx(e+1,1)=nearestPoint(1,2);
          swzy(e+1,1)=nearestPoint(1,1);
        else
          swzx(e+1,1)=wzx(e+1,1);
          swzy(e+1,1)=wzy(e+1,1);
        end
        s(e+1)=sqrt((swzx(e+1,1)-wzx(e+1))^2+(swzy(e+1,1)-wzy(e+1))^2);
       if s(e+1)<=1
         dx=wzx(e+1,1)-origin(1,1);
         dy=origin(1,2)-wzy(e+1,1); 
         ss(e+1)=1;
       else
         dx=swzx(e+1,1)-origin(1,1);
         dy=origin(1,2)-swzy(e+1,1);
         ss(e+1)=0;
       end
        gabor_rad(e+1,1)=atan2(dx,dy);       
%%
       zswzx(e+1,1)=origin(1,1)+R*sin(gabor_rad(e+1,1));
       zswzy(e+1,1)=origin(1,2)-R*cos(gabor_rad(e+1,1));
    else
       dx=wzx(e+1,1)-origin(1,1);
       dy=origin(1,2)-wzy(e+1,1);
       gabor_rad(e+1,1)=atan2(dx,dy);       
       zswzx(e+1,1)=origin(1,1)+R*sin(gabor_rad(e+1,1));
       zswzy(e+1,1)=origin(1,2)-R*cos(gabor_rad(e+1,1));
    end     

    if gabor_rad(e+1,1)>=0 %角度数值不影响迭代 影响卷积 故此处加判断条件
        gabor_rad(e+1,1)=gabor_rad(e+1,1);
    else
        gabor_rad(e+1,1)=gabor_rad(e+1,1)+pi;
    end
%% 新卷积模板 gabor_mean
x = 0;
for i = linspace(-ym/2,ym/2,ym)
    x = x + 1;
    y = 0;
    for j = linspace(-ym/2,ym/2,ym)
        y = y + 1;
        z(y,x)=compute1(i,j,gabor_rad(e+1,1));
    end
end
% 新卷积过程 e
for xjj=e+1:e+2 %起始5-9；9-13…
    filtered{xjj}=conv2(I(:,:,xjj),z,'same');
    fai1{xjj}=angle(filtered{xjj});
end
end
toc
plot(gabor_rad)
% save case3.mat gabor_rad

function gabor_k = compute1(x,y,theta)
lambda=45;
b=0.5;
sigma=(1/pi)*sqrt(log10(2)/2)*((2^b+1)/(2^b-1))*lambda;
psi=0;
gamma=0.5;
x1 = x*cos(theta) + y*sin(theta);
y1 = -x*sin(theta) + y*cos(theta);
gabor_k = exp(-(x1^2/sigma^2+gamma^2*y1^2/sigma^2))*exp(i*(2*pi*x1/lambda+psi));
end
