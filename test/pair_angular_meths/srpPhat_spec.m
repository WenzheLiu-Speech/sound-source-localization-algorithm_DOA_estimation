function spec = srpPhat_spec(X, f, tauGrid)

X1 = X(:,:,1);
X2 = X(:,:,2);
[nbin,nFrames] = size(X1);
ngrid = length(tauGrid);

P = X1.*conj(X2);
P = P./abs(P);
spec = zeros(nbin,nFrames,ngrid);
for pkInd = 1:ngrid
    EXP = repmat(exp(-2*1i*pi*tauGrid(pkInd)*f),1,nFrames);
    spec(:,:,pkInd) = real(P).*real(EXP) - imag(P).*imag(EXP); % 比直接spec(:,:,pkInd) = real(P.*EXP)计算速度更快
end

end