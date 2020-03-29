% ISM算法计算宽带信号DOA
clc  
clear all  
close all
%% 参数
M = 12;      %阵元数  
N = 200;     % 快拍数
fs = 1000;
f0 = 100;    %中心频率   
f1 = 80;     % 最低频率
f2 = 120;    %最高频率
c = 1500;                                    
lambda = c/f0;                               
d = lambda/2; %阵元间距 
J = 33;     % 子带数
SNR = 15;
num_src = 2; % 声源数
theat1 = 30*pi/180;                              
theat2 = 40*pi/180;                              
n = 0:1/fs:N/fs;                              
theat = [theat1 theat2]';  
dtheta = 0.5; 
thetas = -90:dtheta:90;
%% 生成信号  
nfft = 2048;  
s1 = chirp(n,80,1,120);                 % (t,f0,t1,f1) f0是0时刻的瞬时频率，f1是t1时刻的瞬时频率   
sa = fft(s1,nfft);                      

s2 = chirp(n+0.100,80,1,120);                
sb = fft(s2,nfft);                           
   
%% ISM算法 
P=1:num_src;%两个角度  
startfftindex = (nfft/2)/(fs/2)*f1; % 1024(点)->fs/2=500(Hz)；1024/500*80=163.84点    120Hz->245.76点   
endfftindex = (nfft/2)/(fs/2)*f2;
dfftindex = (endfftindex-startfftindex+1)/J;% 每个子带的近似频点数
a=zeros(M,num_src);  
sump=zeros(1,length(thetas));  
for i=1:J %每个子带
    %% 1. 计算每个阵元接收信号的一个子带的协方差矩阵
    fftindex = ceil(startfftindex+(i-1)*dfftindex+1);%该子带内任意一个频点
    f=fftindex/((nfft/2)/(fs/2)); % 该频点对应频率
    s=[sa(fftindex) sb(fftindex)]';  
    for m=1:M  
        a(m,P)=exp(-j*2*pi*f*d/c*sin(theat(P))*(m-1))'; %阵列流形 
    end  
    R=a*(s*s')*a';  % 每个频点的协方差矩阵
    %% 2. 特征值分解，保留噪声子空间
    [em,zm]=eig(R);  % 特征值分解 em:特征向量   zm：特征值对角矩阵
    [zm1,pos1]=max(zm);  %zm1 所有特征值
    for l=1:num_src % 去除两个声源对应的特征值和特征向量，只保留噪声子空间em 
        [zm2,pos2]=max(zm1);  
        zm1(:,pos2)=[];  
        em(:,pos2)=[];  
    end  
    %% 3. 计算p(k)=A'*em*em'*A
    k=1;  %离散角度索引
    for ii=-90:dtheta:90  
        arfa=sin(ii*pi/180)*d/c;  
        for iii=1:M  
            tao(1,iii)=(iii-1)*arfa;  
        end  
        A=[exp(-j*2*pi*f*tao)]';  
        p(k)=A'*em*em'*A;  
        k=k+1;  
    end  
    %% 4. 所有子带的p(k)加和sump
    sump=sump+abs(p);  
end 
%% 5. 计算空间谱 pm = 1/((1/J)*sump)
pmusic=1/J*sump;  
pm=1./pmusic; 

plot(thetas,20*log(abs(pm)));  
xlabel('入射角/度');  
ylabel('空间谱/dB');  
grid on 
