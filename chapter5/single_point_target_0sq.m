% author: Zhao Fei 
% Date: 2019-10-16
% Description: simulation of single point target(zero squint), show the result in time
% , range-dopler, 2d frequence field. 
% This program references 《合成孔径雷达成像算法与实现》, p125.
close all;clear all;
%% simulation parameters (references p124, table 5.2)
Fr = 7.5e6; % range sample rate
Tr = 25e-6;
Kr = 0.25e12;   % range chirp rate, multiply -1 if negative chrip-rate is needed
Fa = 104;   % azimuth sample rate
f0 = 5.3e9; % radar work rate = 5.3 GHz
c = 3e8;    % light speed
lambda = c / f0;
N_rg = 256; % the smaples of range
N_az = 256;
A0 = 1;
Vr = 150;   % m/s
Vg = Vr;
Vs = Vr;
theta_rc = 0;
delta_fdop = 80;
La = 0.886 * 2 * Vs * cos(theta_rc) / delta_fdop;
beta_bw = 0.886 * lambda / La;    % we suppose
R_etac = 20e3;  % 景中心斜距
R0 = R_etac;  % because theta_sqc = 0
Ka = 2 * Vr.^2 / lambda / R0;   % azimuth chrip rate

%% construct the signal
% the time scale
tau = ((-N_rg / 2) : (N_rg / 2 - 1)) / Fr + 2 * R0 / c;
eta = ((-N_az / 2 : N_az / 2 - 1)) / Fa;
[tauX, etaY] = meshgrid(tau, eta);

% the envelope
R_eta = R0 + Vr.^2 * etaY.^2 / (2*R0);
w_r = (abs(tauX - 2 * R_eta / c) <= Tr / 2);
w_a = sinc(0.886 / beta_bw * atan(Vg / R0 * etaY)).^2;  % we know eta_c = 0

% the phase
theta1 = -1j * 4 * pi * R0 / lambda;
theta2 = -1j * pi * Ka * etaY.^2;
theta3 = 1j * pi * Kr * (tauX - 2 * R_eta / c).^2;

% references 《合成孔径雷达成像算法与实现》 p105， (5.2)
s0 = A0 * w_r .* w_a .* exp(theta1) .* exp(theta2) .* exp(theta3);
s0_neg = A0 * w_r .* w_a .* exp(theta1) .* exp(theta2) .* exp(-theta3);

%% show the result in range-azimuth time
Amp = abs(s0);
Phi = angle(s0);
Phi_neg = angle(s0_neg);
figure; % show amplitude， phase of the signal
subplot(131);imagesc(Amp);xlabel('距离向（采样点）');ylabel('方位向（采样点）');title('(a)幅度');
subplot(132);imagesc(Phi);xlabel('距离向（采样点）');ylabel('方位向（采样点）');title('(b)相位（正扫频）');
subplot(133);imagesc(Phi_neg);xlabel('距离向（采样点）');ylabel('方位向（采样点）');title('(c)相位（负扫频）');
suptitle('零斜视角下单目标点时域信号');


%% show the result in range-dopler field
s_rd = fft(s0);
s_rd_neg = fft(s0_neg);
Amp_rd = abs(s_rd);
Phi_rd = angle(s_rd);
Phi_rd_neg = angle(s_rd_neg);
figure; % show amplitude， phase of the signal
subplot(131);imagesc(Amp_rd);xlabel('距离向（采样点）');ylabel('方位向频率（采样点）');title('(a)幅度');set(gca, 'YDir', 'normal');
subplot(132);imagesc(Phi_rd);xlabel('距离向（采样点）');ylabel('方位向频率（采样点）');title('(b)相位（正扫频）');set(gca, 'YDir', 'normal');
subplot(133);imagesc(Phi_rd_neg);xlabel('距离向（采样点）');ylabel('方位向频率（采样点）');title('(c)相位（负扫频）');set(gca, 'YDir', 'normal');
suptitle('零斜视角下单目标点距离多普勒域信号');
set(gca, 'YDir', 'normal');

%% show the result in 2d frequnce field
s_2d = fft2(s0);
s_2d_neg = fft2(s0_neg);
Amp_2d = abs(s_2d);
Phi_2d = angle(s_2d);
Phi_2d_neg = angle(s_2d_neg);
figure; % show amplitude， phase of the signal
subplot(131);imagesc(Amp_2d);xlabel('距离向频率（采样点）');ylabel('方位向频率（采样点）');title('(a)幅度');set(gca, 'YDir', 'normal');
subplot(132);imagesc(Phi_2d);xlabel('距离向频率（采样点）');ylabel('方位向频率（采样点）');title('(b)相位（正扫频）');set(gca, 'YDir', 'normal');
subplot(133);imagesc(Phi_2d_neg);xlabel('距离向频率（采样点）');ylabel('方位向频率（采样点）');title('(c)相位（负扫频）');set(gca, 'YDir', 'normal');
suptitle('零斜视角下单目标点二维频域信号');
