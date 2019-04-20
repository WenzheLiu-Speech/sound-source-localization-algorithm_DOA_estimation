clc,clear all
close all
%%% 互相关法求时延
N=1000;  %长度
Fs=500;  %采样频率
n=0:N-1;
t=n/Fs;   %时间序列
%%
%造CW信号
%信号1
x = 0.7*sin(2*pi*10*t) ;
x = [x,zeros(1,2*N)];%要保证时延在两倍信号长度之内，不然就调整2*N的2

%信号2
y = 0.7*sin(2*pi*10*t);
tao = 1.05;%延时
Ndelay = fix(tao*Fs);%换算成延时点数
y = [zeros(1,Ndelay),y,zeros(1,length(x)-Ndelay-length(y))];%,zeros(1,500)];

%互相关函数
[c,lags]=xcorr(x,y);%互相关
subplot(211);
xaxis = 1:1:length(x);
plot(xaxis,x,'r');
hold on;
plot(xaxis,y,':');
legend('x信号', 'y信号');
xlabel('时间/s');ylabel('x(t) y(t)');
title('原始信号');grid on;
hold off

subplot(212);
plot(lags/Fs,c,'r');
title('相关函数');
xlabel('时间(s)');ylabel('Rxy');
grid on
[Am,Lm]=max(c);
d = Lm - (length(c)+1)/2;
phy=(2*10*d*180/Fs);
phy = rem(phy,360)%取余360
Delay=d/Fs
