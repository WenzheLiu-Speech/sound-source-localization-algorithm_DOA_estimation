%% SRP Estimate of Direction of Arrival at Microphone Array
% Estimate the direction of arrival of a signal using the SRP-PHAT
% algorithm. 
%%

% x = filter(Num,1,x0);
c = 340.0;
d = 0.0420; %圆阵半径
% more test audio file in ../../TestAudio/ folder
path = '../../TestAudio/XMOS/room_mic5-2/';
[s1,fs] = audioread([path,'音轨-2.wav']);
s2 = audioread([path,'音轨-3.wav']);
s3 = audioread([path,'音轨-4.wav']);
s4 = audioread([path,'音轨-5.wav']);
s5 = audioread([path,'音轨-6.wav']);
s6 = audioread([path,'音轨-7.wav']);

signal = [s1,s2,s3,s4,s5,s6]; % 所有阵元接收语音合成一个接收信号矩阵，nsamples x channels Matrix
M = size(signal,2);%阵元数（通道数）
%%
t = 0;

% minimal searching grid
step = 1;

P = zeros(1,length(0:step:360-step));
tic
h = waitbar(0,'Please wait...');
for i = 0:step:360-step
    % Delay-and-sum beamforming
    [ DS, x1] = DelaySumURA(signal,fs,512,512,256,d,i/180*pi);
    t = t+1;
    %beamformed output energy
    P(t) = DS'*DS;
    waitbar(i / length(step:360-step))
end
toc
close(h) 
[m,index] = max(P);
figure,plot(0:step:360-step,P/max(P))
ang = (index)*step

%[ DS, x1] = DelaySumURA(signal,fs,1024,1024,512,d,(index)*step/180*pi);
% audiowrite('DS.wav',real(DS),fs)
% audiowrite('signal1.wav',signal(:,1),fs)

% [ z ] = postprocessing(x1,DS,fs,(index)*step);
% audiowrite('z9.wav',z,fs)
