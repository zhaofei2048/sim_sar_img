% Author: Zhao Fei 
% Date: 2019-2-16
% Description: 多视处理仿真 
% This program references 《合成孔径雷达成像算法与实现》, chapter 6.
close all;clear all;
%% 1. 仿真参数 (部分参考 p142, table 6.1)
R_etac = 20e3;  % 景中心斜距
Vr = 150;   % 等效雷达速度
Tr = 2.5e-6;    % 发射脉冲时宽
Kr = 20e12; % 距离调频率
f0 = 5.3e9; % 雷达工作频率
delta_fdop = 60;    % 多普勒带宽
Fr = 60e6;  % 距离采样率
Fa = 100;   % 方位采样率
Naz = 512;  % 方位向采样点数（距离线条数）
Nrg = 1;  % 距离向采样点数（距离线采样点数）
theta_rc_deg = 0; % 低斜视角0度
eta_c = 0;   % 波束中心偏移时间
f_etac = 0;   % 多普勒中心频率
tc = 0; % tc为所选距离向脉冲零频时刻相对脉冲中心的时延：s(t)=rect(t/T)exp{j*pi*K*(t-tc)^2}
c = 3e8;    % 光速
F_L = 20;   % 子视滤波器带宽
overlap = 20e-2; % 子视20%的重叠
kaiser_beta = 2;    % 子视窗是beta=2的kaiser窗
Num_son = 3;    % 子视数目

% 导出参数
lambda = c / f0;
theta_rc = theta_rc_deg * pi / 180;
Vs = Vr;
Vg = Vr;
La = 0.886 * 2 * Vs * cos(theta_rc) / delta_fdop;
beta_bw = 0.886 * lambda / La;    % we suppose
R0c = R_etac * cos(theta_rc);   % 景中心斜距对应最短斜距（距离徙动校正后图像中心对应斜距）
pr = c/2/Fr;    % 距离向采样间距
Ta = lambda * R0c * delta_fdop / 2 / Vr^2 / cos(theta_rc)^3;    % 目标照射时间
Np = round(Tr * Fr);   % 脉冲序列长度（采样点数）
Npa = round(Ta * Fa);  % 方位向脉冲长度
% 设定3个点目标：A, B, C的位置（最短斜距（相对景中心距离对应最短斜距），方位距离(相对A点)）
% A(0m, -315), B(0m, 0m), C(0m, 315m)，
% 这三个点在同一距离上的不同方位处
NUM_TARGETS = 3;    % 仿真的目标数为3
delta_R0 = [0, 0, 0];  % 目标相对斜距
R_az = [-315, 0, 315]; % 目标相对方位向距离

% 以C点的零多普勒时刻作为方位时间0点，则A,B,C三点的绝对零多普勒时刻分别为：
eta_ca = zeros(1, NUM_TARGETS);
for i = 1:NUM_TARGETS
    eta_ca(i) = R_az(i) / Vr;
end

% 设定原始数据观测时间轴及范围：距离向以景中心斜距所对应距离时间为观测时间中心；
% 方位向以A点的波束中心穿越时刻为观测时间中心（A点的绝对零多普勒时刻记为0）
tau = R_etac * 2 / c;
eta = ((-Naz / 2 : Naz / 2 - 1)) / Fa + eta_c;

%% 2. 构造雷达原始数据（A,B,C三点最短斜距相同，方位距离不同）
% 时间坐标网格
[tauX, etaY] = meshgrid(tau, eta);
% 计算A,B,C三点相对各自零多普勒时刻的方位向时间eta
etaYs = zeros(NUM_TARGETS, Naz, Nrg);
for i = 1:NUM_TARGETS
    etaYs(i,:, :) = etaY - eta_ca(i);
end
% 计算A，B，C三点的最近斜距R0（相当于放置A,B,C三点）
R0 = R0c + delta_R0;
% 计算A，B,C三点瞬时斜距
R_eta = zeros(NUM_TARGETS, Naz, Nrg);
for i = 1:NUM_TARGETS
    R_eta(i, :, :) = (R0(i)^2 + Vr^2 * etaYs(i,:,:).^2 ).^0.5;
end

A0 = 1;
s0 = zeros(Naz, Nrg);
for i = 1:NUM_TARGETS
    % 包络
    w_r = (abs(tauX - 2 * reshape(R_eta(i, :, :), Naz, Nrg) / c) <= Tr / 2);
    % w_a = sinc(0.886 / beta_bw * atan(Vg * (reshape(etaYs(i, :, :), Naz, Nrg) - eta_c) / R0(i))).^2;
    w_a = (abs(reshape(etaYs(i, :, :), Naz, Nrg)-eta_c) < Ta/2);
    % 相位
    theta1 = -4 * pi * f0 * reshape(R_eta(i, :, :), Naz, Nrg) / c;
    theta2 = pi * Kr * (tauX - 2 * reshape(R_eta(i, :, :), Naz, Nrg) / c).^2;
    % 信号多点累加
    s0 = s0 + A0 * w_r .* w_a .* exp(1j*theta1) .* exp(1j*theta2);
end

%% 3. 单视处理：方位向压缩（方式3）
S0 = fft(s0);   % 变换到方位频域
Ka = 2 * Vr^2 / lambda ./ R0c;  % 方位调频率
f_eta = (ifftshift((-Naz/2 : Naz/2-1) * Fa / Naz)).';   % 生成方位频率轴，这里多普勒中心频率为0，所以不怕轨迹卷绕
f_eta = f_eta + round((f_etac - f_eta) / Fa) * Fa;
% 匹配滤波
Haz = exp(-1j*pi*f_eta.^2./Ka);
Sac = S0 .* Haz;
s_ac = ifft(Sac);

%  弃置区
TAL = (Npa-1)/2;    % 左，因为多普勒中心频率为0，因此N_ZD=0
TAR = (Npa-1)/2;    % 右

figure;
subplot(211);plot(real(s0));set(gca,'YLim',[-1.5,1.5]);ylabel('幅度');title('原始数据实部');
subplot(212);plot(real(abs(s_ac)));xlabel('方位向时间（采样点）');ylabel('幅度');title('压缩后的数据（单视）');
suptitle('单视处理');

%% 4. 多视处理
% 子视抽取滤波器位置与方位信号频谱之间的关系
f_eta_normal = (linspace(-Fa/2,Fa/2,Naz).');
S_amp = fftshift(abs(S0));  % 信号fftshift后的幅度谱

% 1. 构造子视窗
Lsvw = round(F_L / Fa * Naz * (1 + overlap));
Wkaiser = kaiser(Lsvw, kaiser_beta);
Wsv = zeros(Naz, Num_son);
Fcen = [-F_L, 0, F_L];
for i = 1:Num_son
    si = floor(Fcen(i) / Fa * Naz + Naz/2 - Lsvw / 2);
    Wsv(si: si+Lsvw-1, i) = Wkaiser;
end
figure;
subplot(211);plot(f_eta_normal, S_amp);ylabel('幅度');title('(a)信号频谱');
subplot(212);plot(f_eta_normal, Wsv(:,1), f_eta_normal,Wsv(:,2), f_eta_normal,Wsv(:,3));xlabel('方位频率（采样点）');ylabel('幅度');title('(b)子视抽取滤波器的位置');
legend('子视1','子视2','子视3');
suptitle('子视抽取滤波器位置与方位信号频谱之间的关系');

% 2. 子视抽取及复数据检测
s1 = ifft(S0 .* Haz .* ifftshift(Wsv(:,1)));
s2= ifft(S0 .* Haz .* ifftshift(Wsv(:,2)));
s3 = ifft(S0 .* Haz .* ifftshift(Wsv(:,3)));
T_L = F_L / Ka; % 子视滤波器时域长度
Lsvt = round(T_L * Fa); % 子视滤波器时域采样点数
abandon_s1 = ones(Naz, 1);  % 弃置区掩膜
abandon_s2 = ones(Naz, 1);
abandon_s3 = ones(Naz, 1);
delta_eta1_N = round((Fcen(1)-f_etac)/abs(Ka) * Fa);   % 子视弃置区延时点数
delta_eta3_N = round((Fcen(3)-f_etac)/abs(Ka) * Fa);   % 子视弃置区延时点数
abandon_s1(end-(Lsvt-1)+1+delta_eta1_N:end+delta_eta1_N) = zeros(Lsvt-1,1);
s2_AL = round((Lsvt-1)/2);
s2_AR = Lsvt-1-s2_AL;
abandon_s2(1:s2_AL) = zeros(s2_AL, 1);
abandon_s2(end-s2_AR+1:end) = zeros(s2_AR, 1);
abandon_s3(1+delta_eta3_N:Lsvt-1+delta_eta3_N) = zeros(Lsvt-1,1);

% 3. 弃置区丢弃
s1_a = abs(s1).* abandon_s1;
s2_a = abs(s2).* abandon_s2;
s3_a = abs(s3).* abandon_s3;

% 4. 子视求和
s_a = sqrt(abs(s1_a).^2 + abs(s2_a).^2 + abs(s3_a).^2);


figure;
subplot(421);plot(abs(s1));ylabel('幅度');title('(a)压缩后的子视1');
subplot(422);plot(s1_a);title('(b)数据舍弃后的子视1');
subplot(423);plot(abs(s2));ylabel('幅度');title('(c)压缩后的子视2');
subplot(424);plot(s2_a);title('(d)数据舍弃后的子视2');
subplot(425);plot(abs(s3));xlabel('方位时间（采样点）');ylabel('幅度');title('(e)压缩后的子视3');
subplot(426);plot(s3_a);title('(f)数据舍弃后的子视3');
subplot(428);plot(s_a);xlabel('方位时间（采样点）');ylabel('幅度');title('(g)子视求和');

% figure
% plot(abandon_s1)
% hold on
% plot(abandon_s2)
% hold on
% plot(abandon_s3)

% 5. 目标E及时域子视滤波器
% 截取目标E
s_E = s0(170:335,1);
figure;
subplot(411);plot(real(s_E));ylabel('幅度');title('(a)目标E的实部');
subplot(412);plot(real(ifft(Haz .* ifftshift(Wsv(:,1)))));ylabel('幅度');title('(b)子视1');
subplot(413);plot((real(ifft(Haz .* ifftshift(Wsv(:,2))))));ylabel('幅度');title('(c)子视2');
subplot(414);plot(real(ifft(Haz .* ifftshift(Wsv(:,3)))));xlabel('方位时间（采样点）');ylabel('幅度');title('(d)子视3');
