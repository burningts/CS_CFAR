%%
function [u1,s] = genFMCW(n,r,signal_SNR)
    %工作信号参数
    f0   = 77e9;                                %工作频率，单位：Hz
    fsamp= 20e6;                               %采样频率，单位：Hz
    ts = 1/fsamp; 
    ti = (0:n-1)*ts; 
    band = 3.6e9;                                %扫频带宽,单位：Hz     
    tp = 60e-6;                              %扫频周期,单位：s  
    u = band/tp; 
    tao = 2.*r./3e8;                             %时延
    freq_gap = u * tao;  %频差

    % 产生信号
    xs = zeros(length(r),n);
    s = zeros(1,n);  %生成的中频IF信号s
    for i = 1:length(r)
        xs(i,:) = exp(1j*2*pi*(u*tao(i)*ti - 0.5*u*tao(i)^2 + f0*tao(i)));
        s = s + xs(i,:);
    end
    s = awgn(0.05*s,signal_SNR);
    fsampu = (0:n/2-1)*(fsamp/n);%取前1024个点
    r = fsampu*3e8/(2*u);
    u = fft(real(s));
    u1 = abs(2*u(1:n/2));

    figure, 
    subplot(1,2,1);
    plot(real(s));
    grid on;
    title('中频IF信号时域波形');
    subplot(1,2,2);
    plot(r,20*log10(u1/max(u1)));       %横轴以MHz为单位，纵轴是dB形式
%     xlim([0,25]);
%     ylim([-100,10]);
    title('中频IF信号频谱');
    xlabel('Range(m)');
    ylabel('Magnitude(dB)');
    grid on;     
end

