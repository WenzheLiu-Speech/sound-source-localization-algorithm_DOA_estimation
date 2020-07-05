function specGlobal = doa_music(x,Param,nsrc)
if(size(x,2)<2)
    error('ERROR[MUSIC]:信号通道数必须大于等于2');
end
%% STFT
X = ssl_stft(x.',Param.window,Param.noverlap,Param.nfft,Param.fs);%nbin,nfram,nchan
X = X(2:end,:,:);
X = X(Param.freqBins,:,:);
[nbin,~,nmic] = size(X);
%% MUSIC
% linspace包含端点，保证插值时不会出现NaN
aziGrid = linspace(Param.azimuth(1),Param.azimuth(end),round((Param.azimuth(end)-Param.azimuth(1))/Param.alphaRes)+1);
eleGrid = linspace(Param.elevation(1),Param.elevation(end),round((Param.elevation(end)-Param.elevation(1))/Param.alphaRes)+1);
power = zeros(nbin, length(aziGrid), length(eleGrid));

for ibin = 1:nbin % 对于每个频点 
    Rxx = (transpose(squeeze(X(ibin,:,:)))*conj(squeeze(X(ibin,:,:))));% 自相关矩阵
    [U,~,~] = svd(Rxx);    % SVD分解   Rxx = U * S * U^H
    En = U(:,nsrc+1:end);  % 噪声子空间
    fprintf('%d\n',ibin)
    for iaz = 1 :length(aziGrid)
        for iel = 1 :length(eleGrid)
            v = [cosd(eleGrid(iel))*cosd(aziGrid(iaz));cosd(eleGrid(iel))*sind(aziGrid(iaz));sind(eleGrid(iel))];% 3 x 1
            tau = v'*(Param.micPos-repmat(Param.micPos(:,1),[1,nmic]))./Param.c; % 1 * nmic  参考麦克为1：Param.micPos(:,1)         
            a = exp(1i*2*pi*Param.f(ibin).*transpose(tau));%nmic x 1     SV = exp(-2*1i*pi*tau*Param.f.');  % nmic x nbin
            power(ibin,iaz,iel) = 1./(sum(abs( ctranspose(a) * En * ctranspose(En) * a )));
        end
    end
end

% 对所有频率的空间谱加在一起:
spec = squeeze(sum(power,1)); %nAzi x nEle
[az,el]=meshgrid(aziGrid,eleGrid);
[azi,eli]=meshgrid(Param.azimuth,Param.elevation);

specInterp = interp2(az,el,spec.',azi,eli);
% specInterp = interp2(azOri,elOri,spec.',azInterp,elInterp);
specGlobal = reshape(specInterp.',1,[]);
end

function X=ssl_stft(x,window,noverlap,nfft,fs)

% Inputs:x: nchan x nsampl  window = blackman(wlen);
% Output:X: nbin x nfram x nchan matrix 

[nchan,~]=size(x);
[Xtemp,F,T,~] = spectrogram(x(1,:),window,noverlap,nfft,fs);%S nbinxnframe
nbin = length(F);
nframe = length(T);
X = zeros(nbin,nframe,nchan);
X(:,:,1) = Xtemp;
for ichan = 2:nchan
    X(:,:,ichan) = spectrogram(x(ichan,:),window,noverlap,nfft,fs); 
end

end
