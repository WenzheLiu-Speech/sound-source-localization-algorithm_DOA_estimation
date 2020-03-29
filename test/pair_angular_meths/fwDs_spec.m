function spec = fwDs_spec(hatRxx, f, d, tauGrid, c)

[nbin,nFrames] = size(hatRxx(:,:,1,1));
ngrid = length(tauGrid);
R11 = hatRxx(:,:,1,1);
R12 = hatRxx(:,:,1,2);
R22 = hatRxx(:,:,2,2);
TR = real(R11 + R22);
SINC = sinc(2*f*d/c);

SNR = zeros(nbin,nFrames,ngrid);
for pkInd=1:ngrid,
    EXP = repmat(exp(-2*1i*pi*tauGrid(pkInd)*f),1,nFrames);
    SNR(:,:,pkInd) = repmat(-(1+SINC)/2,1,nFrames) + repmat((1-SINC)/2,1,nFrames).*(TR + 2*real(R12.*EXP))./(TR - 2*real(R12.*EXP));
end
spec = SNR;

end