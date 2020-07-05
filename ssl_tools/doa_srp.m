function [specGlobal] = doa_srp(x,method, Param)
%% 
if(~any(strcmp(method, {'SRP-PHAT' 'SRP-NON'})))
    error('ERROR[doa_srp]: method参数错误');   
end
%% STFT
X = ssl_stft(x.',Param.window, Param.noverlap, Param.nfft, Param.fs);
X = X(2:end,:,:);
%% 
if strcmp(method,'SRP-PHAT')
    specGlobal = ssl_srpPhat(X,Param);
else
    specGlobal = ssl_srp_nonlin(X,Param);
end

end

function X=ssl_stft(x,window,noverlap,nfft,fs)

% Inputs:x: nchan x nsampl  window = blackman(wlen);
% Output:X: nbin x nfram x nchan matrix 

[nchan,~]=size(x);
[Xtemp,F,T,~] = spectrogram(x(1,:),window,noverlap,nfft,fs); % S nbin x nframe
nbin = length(F);
nframe = length(T);
X = zeros(nbin,nframe,nchan);
X(:,:,1) = Xtemp;
for ichan = 2:nchan
    X(:,:,ichan) = spectrogram(x(ichan,:),window,noverlap,nfft,fs); 
end

end

function [specGlobal] = ssl_srpPhat(X,Param)
[~,nFrames,~] = size(X);
specInst = zeros(Param.nGrid, nFrames);

for i = 1:Param.nPairs
    spec = srpPhat_spec(X(Param.freqBins,:,Param.pairId(i,:)), Param.f(Param.freqBins), Param.tauGrid{i}); % NV % [freq x fram x local angle for each pair]
    specSampledgrid = (shiftdim(sum(spec,1)))';
    specCurrentPair = interp1q(Param.alphaSampled{i}', specSampledgrid, Param.alpha(i,:)');
    specInst(:,:) = specInst(:,:) + specCurrentPair;
end

switch Param.pooling
    case 'max'
        specGlobal = shiftdim(max(specInst,[],2));
    case 'sum'
        specGlobal = shiftdim(sum(specInst,2));
end
end

function [specGlobal] = ssl_srp_nonlin(X,Param)

alpha_meth = (10*Param.c)./(Param.d*Param.fs);
[~,nFrames,~] = size(X);
specInst = zeros(Param.nGrid, nFrames);

for i = 1:Param.nPairs
    spec = srpNonlin_spec(X(Param.freqBins,:,Param.pairId(i,:)), Param.f(Param.freqBins), alpha_meth(i), Param.tauGrid{i});
    specSampledgrid = (shiftdim(sum(spec,1)))';
    specCurrentPair = interp1q(Param.alphaSampled{i}', specSampledgrid, Param.alpha(i,:)');
    specInst = specInst + specCurrentPair;
end

switch Param.pooling
    case 'max'
        specGlobal = shiftdim(max(specInst,[],2));
    case 'sum'
        specGlobal = shiftdim(sum(specInst,2));
end
end

