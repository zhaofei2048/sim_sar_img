% author: Zhao Fei 
% Date: 2019-2-20
% Description: 匹配滤波器3种实现方式弃置区出现位置的探究 
% This program references 《合成孔径雷达成像算法与实现》, chapter 3.
close all;clear all;

%% 1. 生成被压缩信号序列
% 参数设置
N = 401;    % 每个样本含401点
N_ZD = 60;  % 零频点位于目标中心右侧N_ZD点
GAP = 400;  % 目标间隔400点
K = 90;    % 调频率设为100
T = 1;  % 每个目标持续时间设为1s，假设信号从-T/2持续到T/2
Fs = T * N; % 采样率
t0 = N_ZD / Fs; % 零频出现时间点
fc = -K * t0;   % 计算中心频率
t_target = (-N/2:N/2-1) / Fs;
s_target = exp(1j*(2*pi*fc.*t_target + pi*K.*t_target.^2));
figure;
subplot(2,1,1);
plot(real(s_target));xlabel('时间采样点');ylabel('幅度');
title('单个目标信号时域');
set(gca, 'YLim', [-1.5, 2]);
subplot(2,1,2);
plot(abs(fft(s_target)));xlabel('频率采样点');ylabel('幅度');
title('单个目标信号频域');

% 生成信号序列
% s = [zeros(1,GAP), s_target, zeros(1,GAP), s_target, zeros(1,GAP), s_target];
s = s_target;
S = fft(s);
s_len = size(s, 2);

%% 2. 匹配压缩
% 方式1：将复制脉冲信号时域反褶，补零然后傅里叶变换得到匹配滤波器
h1 = [conj(s_target(end:-1:1)), zeros(1, s_len-N)];
H1 = fft(h1);
h1s = ifft(H1 .* S);


% 方式2：将复制脉冲信号补0，直接傅里叶变换并取复共轭得到匹配滤波器
h2 = [s_target, zeros(1, s_len - N)];
H2 = conj(fft(h2));
h2s = ifft(H2 .* S);


% 方式3：在频域直接生成匹配滤波器，滤波器应和脉冲信号频域的二次相位频谱部分互为共轭，线性相位不必弥补，只是引起时域平移，对压缩效果没有影响
f = (-s_len/2:s_len/2-1) * Fs / s_len;
H3_tmp = exp(1j * pi .* f.^2 / K);
W3 = abs(f+K*t0)<=(abs(K)*T/2); % 这里加这个窗限制的目的在于：（复）信号带宽并不一定等于采样带宽
H3 = W3 .* H3_tmp;    % 滤波器频谱的有效值分布宽度也应等于信号带宽，而不仅仅由f的取值宽度（Fs）决定
H3 = fftshift(H3);
h3s = ifft(H3 .* S);

%% 3. 绘图展示
t = (-s_len/2:s_len/2-1) / Fs;
figure;
subplot(4,1,1);
plot(real(s));ylabel('幅度');
title('信号序列');
set(gca, 'YLim', [-1.5, 2]);

subplot(4,1,2);
plot(abs(h1s));ylabel('幅度');
title('方式1匹配滤波输出结果');

subplot(4,1,3);
plot(abs(h2s));ylabel('幅度');
title('方式2匹配滤波输出结果');

subplot(4,1,4);
plot(abs(h3s));xlabel('时间采样点');ylabel('幅度');
title('方式3匹配滤波输出结果');

figure;
plot(real(H3_tmp));
hold on;
plot(W3);
title('被W3图中红色矩形截断的才是滤波器频谱有效值，之外的为无效值');

