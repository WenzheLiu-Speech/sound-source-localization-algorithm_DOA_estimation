function [specGlobal] = doa_mvdr(x,method,Param)
%%
if(~any(strcmp(method, {'SNR-MVDR' 'SNR-FWMVDR' 'SNR-DS' 'SNR-FWDS'})))
    error('ERROR : method参数错误');   
end
%%
lf=8;lt=2;
Rxx = ssl_Rxx(x,Param.fs,Param.window,Param.noverlap, Param.nfft,lf,lt);
Rxx = permute(Rxx(:,:,2:end,:),[3 4 1 2]); % nbin x nFrames x nChan x nChan
%% 
if strcmp(method,'SNR-MVDR')
    specGlobal = ssl_MVDR(Rxx,Param);
elseif strcmp(method,'SNR-FWMVDR')
    specGlobal = ssl_FWMVDR(Rxx,Param);
elseif strcmp(method,'SNR-DS')
    specGlobal = ssl_DS(Rxx,Param);
elseif strcmp(method,'SNR-FWDS')
    specGlobal = ssl_FWDS(Rxx,Param);
else
    error('ERROR :错误的method');
end
end
function [hatRxx]=ssl_Rxx(x,fs,window,noverlap, nfft,lf,lt)
%%
if nargin<3, error('Not enough input arguments.'); end
[nsampl,nchan]=size(x);
if nchan>nsampl, error('The input signal must be in columns.'); end
if nargin<4, lf=2; end
if nargin<5, lt=2; end

%% STFT 
X = ioa_stftCompute(x.',window,noverlap,nfft,fs);
[nbin,nfram,nchan]=size(X);
%% 
winf=hanning(2*lf-1);
wint=hanning(2*lt-1).';
hatRxx = zeros(nchan,nchan,nbin,nfram);

pairId = nchoosek(1:nchan,2);
[nPairs,~] = size(pairId);

for f=1:nbin,
    for t=1:nfram,
        indf=max(1,f-lf+1):min(nbin,f+lf-1);
        indt=max(1,t-lt+1):min(nfram,t+lt-1);
        nind=length(indf)*length(indt);
        wei=ones(nchan,1)*reshape(winf(indf-f+lf)*wint(indt-t+lt),1,nind);
        XX=reshape(X(indf,indt,:),nind,nchan).';
        local_Cx = (XX.*wei)*XX'/sum(wei(1,:));
        for idPair = 1:nPairs
            hatRxx(pairId(idPair,:),pairId(idPair,:),f,t) = local_Cx(pairId(idPair,:),pairId(idPair,:));
        end
    end
end
end
function [specGlobal] = ssl_DS(hatRxx,Param)
[~,nFrames,~,~] = size(hatRxx); % nbin x nFrames x 2 x 2
specInst = zeros(Param.nGrid, nFrames);

for i = 1:Param.nPairs
    spec = ds_spec(hatRxx(Param.freqBins,:,Param.pairId(i,:),Param.pairId(i,:)), Param.f(Param.freqBins), Param.tauGrid{i}); %
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
function [specGlobal] = ssl_FWDS(hatRxx, Param)

[~,nFrames,~,~] = size(hatRxx); % nbin x nFrames x 2 x 2
specInst = zeros(Param.nGrid, nFrames);

for i = 1:Param.nPairs
    spec = fwDs_spec(hatRxx(Param.freqBins,:,Param.pairId(i,:),Param.pairId(i,:)), Param.f(Param.freqBins), Param.d(i), Param.tauGrid{i},Param.c); %
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
function [specGlobal] = ssl_MVDR(hatRxx,Param)

[~,nFrames,~,~] = size(hatRxx); % nbin x nFrames x nmic x nmic
specInst = zeros(Param.nGrid, nFrames);

for i = 1:Param.nPairs
    spec = mvdr_spec(hatRxx(Param.freqBins,:,Param.pairId(i,:),Param.pairId(i,:)), Param.f(Param.freqBins), Param.tauGrid{i}); % 取一对麦克风进行MVDR
    specSampledgrid = (shiftdim(sum(spec,1)))'; % sum on frequencies
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
function [specGlobal] = ssl_FWMVDR(hatRxx,Param)

[~,nFrames,~,~] = size(hatRxx); % nbin x nFrames x 2 x 2
specInst = zeros(Param.nGrid, nFrames);

for i = 1:Param.nPairs
    spec = fwMvdr_spec(hatRxx(Param.freqBins,:,Param.pairId(i,:),Param.pairId(i,:)), Param.f(Param.freqBins), Param.d(i), Param.tauGrid{i}, Param.c); %
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
