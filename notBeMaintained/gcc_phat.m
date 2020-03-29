clc
clear
%%%
    %用于求gcc-phat的方法
    %Created by Vincent Lau
    %final version: Dec 5, 2018
%%%
x = [1,2,3,7,9,8,3,7];
y = [4,5,6,5,4,3,8,2];
fs = 1;
crossCorrelate = xcorr(x,y);%自带求互相关函数，用做验证

Nfft = 2*length(x) -1 ;%保证Nfft点数大于length(x)+length(y)-1 即可，最好为2^N
%%
%信号和参考信号的fft，并对参考信号取共轭
b1f=fft(x,Nfft);
b2f=fft(y,Nfft);
b2fc=conj(b2f);
%%
%频谱相乘带加权后，fft反变换
neuma=(b1f).*(b2fc);%频谱相乘

deno=abs(neuma);%GCC的加权系数
%deno = 1;%互相关的加权系数

GPHAT=neuma./deno;%频谱相乘带加权
GPHATi=ifft(GPHAT);  %fft反变换

GCCPHAT = fftshift(real(GPHATi));%GCC-PHAT
tao = -length(x)+1:length(x)-1;
tao = tao/fs;
%归一化
crossCorrelate = crossCorrelate./max(crossCorrelate);
GCCPHAT = GCCPHAT./max(GCCPHAT);
plot(tao,crossCorrelate,'r')
hold on
%figure()
plot(tao,GCCPHAT)
