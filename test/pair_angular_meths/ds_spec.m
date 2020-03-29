function spec = ds_spec(hatRxx, f, tauGrid)

[nbin,nFrames] = size(hatRxx(:,:,1,1));
ngrid = length(tauGrid);
R11 = hatRxx(:,:,1,1);
R12 = hatRxx(:,:,1,2);
R22 = hatRxx(:,:,2,2);
traceRxx = real(R11 + R22);

SNR = zeros(nbin,nFrames,ngrid);
for pkInd=1:ngrid,
    EXP = repmat(exp(-2*1i*pi*tauGrid(pkInd)*f),1,nFrames);
    SNR(:,:,pkInd) = (traceRxx + 2*real(R12.*EXP))./(traceRxx - 2*real(R12.*EXP));
end
spec = SNR;

end