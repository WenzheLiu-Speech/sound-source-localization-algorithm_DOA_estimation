clc;clear;close all;

addpath(genpath('./../'));
addpath('./wav files');
%% Input File & Mic config
fileName = 'male_female_mixture.wav';  
micPos = ... 
...%  mic1	 mic2   mic3   mic4   mic5   mic6   mic7  mic8
    [ 0.037 -0.034 -0.056 -0.056 -0.037  0.034  0.056 0.056;  % x
      0.056  0.056  0.037 -0.034 -0.056 -0.056 -0.037 0.034;  % y
    -0.038   0.038 -0.038  0.038 -0.038  0.038 -0.038 0.038]; % z


azBound = [-180 180]; % 方位角搜索范围
elBound = [-90 90];   % 俯仰角搜索范围。若只有水平面：则elBound=0;
gridRes = 1;          % 方位角/俯仰角的分辨率
alphaRes = 5;          % Resolution (? of the 2D reference system defined for each microphone pair

% method = 'SRP-PHAT';
% method = 'SNR-MVDR';
% method = 'SNR-FWMVDR';
method = 'MUSIC';
wlen = 512;
window = hann(wlen);
noverlap = 0.5*wlen;
nfft = 512;
nsrc = 2;          % 声源个数
c = 343;        % 声速
freqRange = [];         % 计算的频率范围 []为所有频率
pooling = 'max';      % 如何聚合各帧的结果：所有帧取最大或求和{'max' 'sum'}
%% 读取音频文件(fix)
[x,fs] = audioread(fileName);
[nSample,nChannel]=size(x);
if nChannel>nSample, error('ERROR:输入信号为nSample x nChannel'); end
[~,nMic,~] = size(micPos);
if nChannel~=nMic, error('ERROR:麦克风数应与信号通道数相等'); end
%% 保存参数(fix)
Param = pre_paramInit(c,window, noverlap, nfft,pooling,azBound,elBound,gridRes,alphaRes,fs,freqRange,micPos);
%% 定位(fix)
if strfind(method,'SRP')
    specGlobal = doa_srp(x,method, Param);
elseif strfind(method,'SNR')
    specGlobal = doa_mvdr(x,method,Param);
elseif strfind(method,'MUSIC')
    specGlobal = doa_music(x,Param,nsrc);
else 
end
% save('n.mat','specGlobal');
% ppfSpec2D = (reshape(specGlobal,length(Param.azimuth),length(Param.elevation)))';
% imagesc(ppfSpec2D)
%% 计算角度
minAngle                   = 10;         % 搜索时两峰之间最小夹角
specDisplay                = 1;          % 是否展示角度谱{1,0}
% pfEstAngles = post_sslResult(specGlobal, nsrc, Param.azimuth, Param.elevation, minAngle);
% 绘制角谱
% [pfEstAngles,figHandle] = post_findPeaks(specGlobal, Param.azimuth, Param.elevation, Param.azimuthGrid, Param.elevationGrid, nsrc, minAngle, specDisplay);
[pfEstAngles,figHandle] = post_findPeaks(specGlobal, Param.azimuth, Param.elevation, Param.azimuthGrid, Param.elevationGrid, nsrc, minAngle, specDisplay);

azEst = pfEstAngles(:,1)';
elEst = pfEstAngles(:,2)';
for i = 1:nsrc
    fprintf('Estimated source %d : \n Azimuth (Theta): %.0f \t Elevation (Phi): %.0f \n\n',i,azEst(i),elEst(i));
end