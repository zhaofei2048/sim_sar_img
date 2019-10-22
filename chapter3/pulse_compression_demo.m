% author: Zhao Fei 
% Date: 2019-10-20
% Description: Pulse compression demo on LFM signal
close all;clear all;
Fs = 100;  % sample rate
f0 = 0; % center frequency
k = 10;    % chrip rate
T = 4;  % keep time = 4s

t_scale = linspace(-T / 2, T / 2, T * Fs);
s_lfm = 1.0 * exp(1j * (2 * pi * f0 .* t_scale + pi * k .* t_scale .^ 2));
hn = conj(flip(s_lfm));

y = conv(s_lfm, hn);
t_scale_out = linspace(-T, T, size(y, 2));

figure;
subplot(121); plot(t_scale, real(s_lfm)); xlabel('t/s'); ylabel('实部');
title('调频率10Hz/s，中心频率为0的复线性调频信号实部');
subplot(122); plot(t_scale_out, real(y)); xlabel('t/s'); ylabel('实部');
title('匹配滤波输出信号实部');
suptitle('复线性调频信号脉冲压缩前后对比(Fs=100Hz）');