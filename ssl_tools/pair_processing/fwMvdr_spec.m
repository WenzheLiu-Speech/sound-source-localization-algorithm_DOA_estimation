function spec = fwMvdr_spec(hatRxx, f, d, tauGrid, c)

[nbin,nFrames] = size(hatRxx(:,:,1,1));
ngrid = length(tauGrid);
R11 = hatRxx(:,:,1,1);
R12 = hatRxx(:,:,1,2);
R21 = hatRxx(:,:,2,1);
R22 = hatRxx(:,:,2,2);
traceRxx = real(R11 + R22);
SINC = sinc(2*f*d/c);

SNR = zeros(nbin,nFrames,ngrid);
for pkInd=1:length(tauGrid),
    EXP = repmat(exp(-2*1i*pi*tauGrid(pkInd)*f),1,nFrames);
    power_y = real(R11.*R22 - R12.*R21)./(traceRxx - 2*real(R12.*EXP));
    SNR(:,:,pkInd) = repmat(-(1+SINC)/2,1,nFrames) + repmat((1-SINC)/2,1,nFrames).*power_y./(.5*traceRxx-power_y);
end
spec = SNR;

end
