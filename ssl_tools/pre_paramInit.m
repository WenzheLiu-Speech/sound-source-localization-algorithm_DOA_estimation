function Param = pre_paramInit(c,window, noverlap, nfft,pooling,azBound,elBound,gridRes,alphaRes,fs,freqRange,micPos)
Param = struct;
%% 
if(isempty(micPos))
    error('ERROR : 请输入micPos');
else
    [dim1,~,~] = size(micPos);
    if(dim1~=3),error('ERROR : micPos必须是三维坐标');end
end

Param.window = window;
Param.noverlap = noverlap;
Param.nfft = nfft;
Param.fs = fs;
Param.f = Param.fs/Param.nfft*(1:Param.nfft/2).';
if(isempty(freqRange))
    Param.freqBins = 1:length(Param.f);
elseif(freqRange(1) < 0 || freqRange(2) > Param.fs/2)
    error('ERROR : 频率范围freqRange应在 0Hz 到 fs/2 之间');
else
    binMin = find(Param.f >= freqRange(1),1,'first');
    binMax = find(Param.f<freqRange(2),1,'last');
    Param.freqBins = binMin:binMax;
end
Param.c = c;
Param.pooling = pooling;
Param.micPos = micPos;
% Param.arrayCentroid = squeeze(mean(Param.micPos,2));

if(isempty(alphaRes))
    Param.alphaRes = 5;
elseif(alphaRes < 0)
    error('ERROR : alphaRes应为正值');
else
    Param.alphaRes = alphaRes;
end
if(isempty(gridRes))
    gridRes = 1;
end

if(isempty(azBound))
    azBound = [-180 180];
elseif((length(azBound) == 1) && azBound >= -90 && azBound <= 90)
    azBound = [azBound,azBound];
elseif(length(azBound) == 2 && azBound(1) >= -180 && azBound(2) <= 180 && azBound(1)<=azBound(2))
    % nothing to do
else
    error('ERROR : azBound输入不合法，应为在-/+ 180范围内的一个标量或一个二维向量');   
end

if(isempty(elBound))
    elBound = [-90 90];
elseif(length(elBound) == 1 && elBound >= -90 && elBound <= 90)
    elBound = [elBound,elBound];
elseif(length(elBound) == 2 && elBound(1) >= -90 && elBound(2) <= 90 && elBound(1)<=elBound(2))
    % nothing to do
else
    error('ERROR : elBound输入不合法，应为在-/+ 90范围内的一个标量或一个二维向量');   
end

if(length(unique(elBound)) == 1 && length(unique(azBound)) == 1)
    error('ERROR : azBound和elBound至多有一个为标量');
end

Param.azimuth = (azBound(1) : gridRes : azBound(2))';
Param.elevation   = (elBound(1) : gridRes : elBound(2));
nAz = length(Param.azimuth);
nEl = length(Param.elevation);
Param.azimuthGrid = repmat(Param.azimuth,nEl,1)';
Param.elevationGrid = (reshape(repmat(Param.elevation,nAz,1),1,nAz*nEl));

%% 将所有候选方位转换为笛卡尔坐标
Param.nGrid = length(Param.azimuthGrid);    % (nAlxnEl) x 1
directionCoordinate = zeros(3,Param.nGrid); % 3 x (nAlxnEl)
[directionCoordinate(1,:), directionCoordinate(2,:), directionCoordinate(3,:)] = sph2cart(Param.azimuthGrid*pi/180, Param.elevationGrid*pi/180, 1);
% 所有的麦克风对都初始化一个所有方位的笛卡尔坐标矩阵 3 x nMicPair x nDirction
micPost = (Param.micPos)';
nMic = size(micPost,1);
Param.pairId = nchoosek(1:nMic,2);
Param.nPairs = size(Param.pairId,1);
coordinate_pair = repmat(directionCoordinate,[1 1 Param.nPairs]);
coordinate_pair = permute(coordinate_pair,[1 3 2]);
%% 所有麦克风对之间的间距
delta12 = micPost(Param.pairId(:,1),:) - micPost(Param.pairId(:,2),:);
Param.d = sqrt(sum(delta12.^2,2));
delta12_pair = repmat(delta12',[1 1 Param.nGrid]);

Param.alpha = real(acosd(shiftdim(sum(coordinate_pair.*delta12_pair),1)./repmat(Param.d,[1 Param.nGrid])));
Param.alphaSampled = cell(1,Param.nPairs);
Param.tauGrid = cell(1,Param.nPairs);
for index = 1:Param.nPairs
    Param.alphaSampled{index} = floor(min(Param.alpha(index,:))/Param.alphaRes) * Param.alphaRes : Param.alphaRes : ceil(max(Param.alpha(index,:))/Param.alphaRes) * Param.alphaRes;
    Param.tauGrid{index} = Param.d(index)*cos(Param.alphaSampled{index}.*pi/180)./Param.c; % 时延
end
end
