function spec = mvdr_spec(hatRxx, f, tauGrid)

[nbin,nFrames] = size(hatRxx(:,:,1,1));
ngrid = length(tauGrid);
R11 = hatRxx(:,:,1,1);
R12 = hatRxx(:,:,1,2);
R21 = hatRxx(:,:,2,1);
R22 = hatRxx(:,:,2,2);
traceRxx = real(R11 + R22);   % tr(hatRxx)

SNR = zeros(nbin,nFrames,ngrid);
for pkInd=1:ngrid,
    EXP = repmat(exp(-2*1i*pi*tauGrid(pkInd)*f),1,nFrames); % d = [1 EXP],EXP=exp(ix)
    power_y = real(R11.*R22 - R12.*R21)./(traceRxx - 2*real(R12.*EXP)); 
    % power_y = (d'(hatRxx^-1)d)^-1，是MVDR的方位谱函数，代表提纯信号的功率(见论文two_decades_of_array_signal_processing_research)
    % 定理：对于二阶矩阵：M = [a b; c d];逆矩阵M^-1 = (1/det(M))*[d -b; -c a]
    % det(hatRxx) = real(R11.*R22 - R12.*R21)  hatRxx^-1 = [R22 -R12; -R21 R11]/det(hatRxx)
    % d'(hatRxx^-1)d = [1 exp(-ix)][R22 -R12; -R21 R11][1;exp(ix)]/det(hatRxx)
    %                =[R22-exp(-ix)R12,-R21+exp(-ix)R11][1;exp(ix)]/det(hatRxx)
    %                = (R22-exp(-ix)R12-R21exp(ix)+R11)/det(hatRxx)
    %                        |  协方差矩阵为实对称阵：R12=R21，trace=R11+R22, EXP=exp(ix)
    %                =(trace-2*R12(cos(x)-isin(x)+cos(x)+isin(x)))/det(hatRxx)
    %                = (trace-2*real(EXP)*R12)/det(hatRxx)
    % 因为d'(hatRxx^-1)d是标量，所以power_y = (d'(hatRxx^-1)d)^-1 = det(hatRxx)/(trace-2*real(EXP)*R12)
    SNR(:,:,pkInd) = power_y./(.5*traceRxx-power_y); % 输出信噪比 SNR_out = 提纯信号功率/残留噪声功率；
                                             % 残留噪声功率 = 接收信号功率-提纯信号功率
                                             % 协方差矩阵对角线为各通道接收信号的方差（功率）,采用trace(Rxx)/nmic(所有通道能量对通道求平均)的方法计算一个通道的接收信号功率，而不是只取任意一个通道方差的方法
                                             % 本程序是逐对pair计算MVDR角谱，所以nmic=2，接收信号功率=0.5*trace(Rxx)
end
spec = SNR;

end