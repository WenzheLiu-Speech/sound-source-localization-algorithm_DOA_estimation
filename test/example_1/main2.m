function [specGlobal] = main2()
addpath(genpath('./../'));
addpath('./wav files');
%% Input File & Mic config
fileName = 'male_female_mixture.wav';  
micPos = ... 
...%  mic1	 mic2   mic3   mic4   mic5   mic6   mic7  mic8
    [ 0.037 -0.034 -0.056 -0.056 -0.037  0.034  0.056 0.056;  % x
      0.056  0.056  0.037 -0.034 -0.056 -0.056 -0.037 0.034;  % y
    -0.038   0.038 -0.038  0.038 -0.038  0.038 -0.038 0.038]; % z


azBound                    = [40 140]; % 方位角搜索范围?
elBound                    = [-5 50];   % 俯仰角搜索范围。若只有水平面：则elBound=0;
gridRes                    = 5;          % 方位角/俯仰角的分辨率
alphaRes                   = 5;          % Resolution (? of the 2D reference system defined for each microphone pair

% method = 'SRP-PHAT';
method = 'SNR-MVDR';
% method = 'SNR-FWMVDR';
wlen = 512;
window = hann(wlen);
noverlap = 0.5*wlen;
nfft = 512;

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
specGlobal = ssl(x,Param);

end
function [specGlobal] = ssl(x,Param)
lf=8;lt=2;
Rxx = ioa_RxxCompute(x,Param.fs,Param.window,Param.noverlap, Param.nfft,lf,lt);
Rxx = permute(Rxx(:,:,2:end,:),[3 4 1 2]); % nbin x nFrames x nChan x nChan
specGlobal = MVDR_spec(Rxx,Param);
end
function [P] = MVDR_spec(hatRxx,Param)
[nbin,nFrames,nmic,nmic] = size(hatRxx); % nbin x nFrames x nmic x nmic
invhatRxx = zeros(nbin,nFrames,nmic,nmic);
for iframe = 1:nFrames
    for ibin = 1:nbin
        invhatRxx(ibin,iframe,:,:)=inv(shiftdim(hatRxx(ibin,iframe,:,:)));
    end
end
P = zeros(length(Param.azimuth),length(Param.elevation),nFrames);
for iaz = 1:length(Param.azimuth)
    for iel = 1:length(Param.elevation)
        az= Param.azimuth(iaz);el=Param.elevation(iel);
        fprintf('%d %d\n',az,el)
        v = [sind(el)*cosd(az);sind(el)*sind(az);cosd(el)];% 3 x 1
        tau = v'*(Param.micPos-repmat(Param.arrayCentroid,[1,nmic]))./Param.c; % 1 * nmic
        tau = tau.';
        a = exp(-2*1i*pi*tau*Param.f.');  % nmic x nbin
        
        for iframe = 1:nFrames
            for ibin = 1:nbin
                Ptemp = 1./(a(:,ibin)'*(shiftdim(invhatRxx(ibin,iframe,:,:)))*a(:,ibin));
                P(iaz,iel,iframe) = P(iaz,iel,iframe) + Ptemp;
            end
        end
     end
 end

P = shiftdim(sum(P,3));
end
