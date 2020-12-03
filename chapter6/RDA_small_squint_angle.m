% Author: Zhao Fei 
% Date: 2019-2-16
% Description: 低斜视角下的距离多普勒算法仿真 
% This program references 《合成孔径雷达成像算法与实现》, chapter 6.
close all;clear all;
%% 1. 仿真参数 (参考 p142, table 6.1)
R_etac = 20e3;  % 景中心斜距
Vr = 150;   % 等效雷达速度
Tr = 2.5e-6;    % 发射脉冲时宽
Kr = 20e12; % 距离调频率
f0 = 5.3e9; % 雷达工作频率
delta_fdop = 80;    % 多普勒带宽
Fr = 60e6;  % 距离采样率
Fa = 100;   % 方位采样率
Naz = 256;  % 方位向采样点数（距离线条数）
Nrg = 256;  % 距离向采样点数（距离线采样点数）
theta_rc_deg = 3.5; % 低斜视角3.5度
eta_c = -8.1;   % 波束中心偏移时间
f_etac = 320;   % 多普勒中心频率
tc = 0; % tc为所选距离向脉冲零频时刻相对脉冲中心的时延：s(t)=rect(t/T)exp{j*pi*K*(t-tc)^2}
c = 3e8;    % 光速

% 导出参数
lambda = c / f0;
theta_rc = theta_rc_deg * pi / 180;
Vs = Vr;
Vg = Vr;
La = 0.886 * 2 * Vs * cos(theta_rc) / delta_fdop;
beta_bw = 0.886 * lambda / La;    % we suppose
Np = Tr * Fr;   % 脉冲序列长度（采样点数）
R0c = R_etac * cos(theta_rc);   % 景中心斜距对应最短斜距（距离徙动校正后图像中心对应斜距）
pr = c/2/Fr;    % 距离向采样间距

% 设定3个点目标：A, B, C的位置（最短斜距（相对景中心距离对应最短斜距），方位距离(相对A点)）
% A(-25m, 0), B(-25m, 25m), C(+25m, 50+BC*tan(theta_rc))，
NUM_TARGETS = 3;    % 仿真的目标数为3
delta_R0 = [-25, -25, 25];  % 目标相对斜距
R_rg = [0, 0, 40];  % 这是假定的BC间的地距，它介于（0,50）之间
R_az = [-50, 0 - R_rg(3)*tan(theta_rc), 0]; % 目标相对方位向距离


% 以C点的零多普勒时刻作为方位时间0点，则A,B,C三点的绝对零多普勒时刻分别为：
eta_ca = zeros(1, NUM_TARGETS);
for i = 1:NUM_TARGETS
    eta_ca(i) = R_az(i) / Vr;
end

% 设定原始数据观测时间轴及范围：距离向以景中心斜距所对应距离时间为观测时间中心；
% 方位向以A点的波束中心穿越时刻为观测时间中心（A点的绝对零多普勒时刻记为0）
tau = ((-Nrg / 2) : (Nrg / 2 - 1)) / Fr + R_etac * 2 / c;
eta = ((-Naz / 2 : Naz / 2 - 1)) / Fa + eta_c;

%% 2. 构造雷达原始数据
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
    w_a = sinc(0.886 / beta_bw * atan(Vg * (reshape(etaYs(i, :, :), Naz, Nrg) - eta_c) / R0(i))).^2;
    % 相位
    theta1 = -4 * pi * f0 * reshape(R_eta(i, :, :), Naz, Nrg) / c;
    theta2 = pi * Kr * (tauX - 2 * reshape(R_eta(i, :, :), Naz, Nrg) / c).^2;
    % 信号多点累加
    s0 = s0 + A0 * w_r .* w_a .* exp(1j*theta1) .* exp(1j*theta2);
end

figure; % 绘制滴斜视角情况下的三点雷达原始仿真信号
subplot(221);imagesc(real(s0));ylabel('方位向时间（采样点）');title('(a)实部');
subplot(222);imagesc(imag(s0));title('(b)虚部');
subplot(223);imagesc(abs(s0));xlabel('距离向时间（采样点）');ylabel('方位向时间（采样点）');title('(c)幅度');
subplot(224);imagesc(angle(s0));xlabel('距离向时间（采样点）');title('(d)相位');
suptitle('3.5度斜视角情况下的三点雷达原始仿真信号（时域）');

%% 3. 距离压缩
% 计算匹配滤波器
tau_p = (-Np/2 : Np/2-1) / Fr;  % 脉冲序列的时间坐标
h_rc = (abs(tau_p) <= Tr/2) .* (kaiser(Np, 2.5).') .* exp(1j * pi * Kr * tau_p.^2);
% 方式1：复制脉冲，反褶并复共轭，补零DFT得到频域滤波器
Hrc1 = repmat(fft(conj(h_rc(end:-1:1)), Nrg), Naz, 1);
% 方式2：复制脉冲，补零DFT，然后取复共轭得到滤波器
Hrc2 = repmat(conj(fft(h_rc, Nrg)), Naz, 1);
% 方式3：根据脉冲频谱特性在频域生成滤波器
f_tau = ifftshift((-Nrg/2:Nrg/2-1) * Fr / Nrg); % 生成距离向频率轴
f_tau = f_tau + round((0 - f_tau) / Fr) * Fr;
Hrc3_tmp = exp(1j * pi .* f_tau.^2 / Kr);
W3 = abs(f_tau+Kr*tc) <= (abs(Kr)*Tr/2); % 这里加这个窗限制的目的在于：（复）信号带宽并不一定等于采样带宽
Hrc3 = W3 .* Hrc3_tmp;    % 滤波器频谱的有效值分布宽度也应等于信号带宽，而不仅仅由f的取值宽度（Fs）决定
Hrc3 = repmat(Hrc3, Naz, 1);

% 在频域匹配滤波
S0 = fft(s0.').';
Src = S0 .* Hrc3 .* repmat(ifftshift(kaiser(Nrg, 2.5).'), Naz, 1);   % 选择方式1或2或3进行距离压缩
s_rc = ifft(Src.').';
N_rg = Nrg;
% N_rg = Nrg-Np+1;  % 丢弃弃置区
% s_rc = s_rc(:,1:N_rg);

figure; % 绘制距离压缩后的结果
subplot(121);imagesc(real(s_rc));xlabel('距离向时间（采样点）');ylabel('方位向时间（采样点）');title('(a)实部');
subplot(122);imagesc(abs(s_rc));xlabel('距离向时间（采样点）');title('(b)幅度');
suptitle('3.5度斜视角距离压缩后信号（时域）');

%% 4. 方位向傅里叶变换
Srd = fft(s_rc);
figure; % 绘制距离多普勒域里的距离压缩后的结果
subplot(121);imagesc(real(Srd));xlabel('距离向时间（采样点）');ylabel('方位向频率（采样点）');title('(a)实部');set(gca, 'YDir', 'normal');
subplot(122);imagesc(abs(Srd));xlabel('距离向时间（采样点）');title('(b)幅度');set(gca, 'YDir', 'normal');
suptitle('3.5度斜视角距离压缩后信号（距离多普勒域）');

%% 5. 距离徙动校正
f_eta = (ifftshift((-Naz/2 : Naz/2-1) * Fa / Naz)).';
f_eta = f_eta + round((f_etac - f_eta) / Fa) * Fa;
R0_tau = (-N_rg/2 : N_rg/2-1) * pr + R0c;
[R0_tau_grid, f_eta_grid] = meshgrid(R0_tau, f_eta);
% 计算距离徙动量矩阵
RCM = lambda^2 * R0_tau_grid .* f_eta_grid.^2 / 8 / Vr^2;
RCM = RCM - (R_etac - R0c); % 将距离徙动量转换到原图像坐标系中（因为原坐标系距离轴以R_etac为中心，最后坐标系以R0c为中心）
RCM = RCM / pr;   % 将距离徙动量转换为距离单元偏移量

% 计算插值核系数表
x_tmp = repmat(-4:3, 16, 1);
offset_tmp = (1:16)/16;
x_tmp = x_tmp + repmat(offset_tmp.', 1, 8);
hx = sinc(x_tmp);
x_tmp16 = x_tmp .* 16;
x_tmp16 = round(x_tmp16 + 16 * 8 / 2);
kwin = repmat(kaiser(16*8, 2.5).', 16, 1);
hx = kwin(x_tmp16) .* hx;
hx = hx ./ repmat(sum(hx, 2), 1, 8);

% 插值校正
Srcmc = zeros(Naz, N_rg);  % 存放距离徙动校正后的回波信号
for i = 1:Naz
    for j = 1:N_rg
        offset_int = ceil(RCM(i,j));
        offset_frac = round((offset_int - RCM(i,j)) * 16);
        if offset_frac == 0
            Srcmc(i,j) = Srd(i,ceil(mod(j+offset_int-0.1,N_rg)));   % 利用信号数据S1的周期性假定
        else
            Srcmc(i,j) = Srd(i, ceil(mod((j+offset_int-4:j+offset_int+3)-0.1,N_rg))) * hx(offset_frac,:).';
        end
        
    end
end

figure; % 绘制距离多普勒域里的距离压徙动校正后的结果
subplot(121);imagesc(real(Srcmc));xlabel('距离向时间（采样点）');ylabel('方位向频率（采样点）');title('(a)实部');set(gca, 'YDir', 'normal');
subplot(122);imagesc(abs(Srcmc));xlabel('距离向时间（采样点）');title('(b)幅度');set(gca, 'YDir', 'normal');
suptitle('3.5度斜视角距离徙动校正后信号（距离多普勒域）');

%% 6. 方位压缩
Ka = 2 * Vr^2 / lambda ./ R0_tau_grid;
Haz = exp(-1j*pi*f_eta_grid.^2./Ka);
Srd_ac = Srcmc .* Haz;

%% 7. 得到时域SAR图像
s_ac = ifft(Srd_ac);

figure; % 绘制低斜视角情况下距离压缩且方位压缩后信号
subplot(121);imagesc(real(s_ac));xlabel('距离向时间（采样点）');ylabel('方位向时间（采样点）');title('(a)实部');
subplot(122);imagesc(abs(s_ac));xlabel('距离向时间（采样点）');title('(b)幅度');
suptitle('3.5度斜视角距离压缩且方位压缩后的信号（时域）');

%% 8. 点目标分析
Srcmc_t = ifft(Srcmc);
c_rd_a = Srcmc_t(:,139);

%1
figure; % 绘制理想距离徙动校正与实际距离徙动校正的对比图
subplot(121);
plot(abs(w_a));xlabel('方位向采样点');ylabel('幅度');title('(a)理想的距离徙动校正');
subplot(122);
plot(abs(c_rd_a));xlabel('方位向采样点');ylabel('幅度');title('(b)8点插值距离徙动校正');
suptitle('距离徙动校正不精确时，由调制引入的成对回波');

%2
c_rd_r = Srcmc(52,139-8:139+7);
c_rd_r_interp = interpft(c_rd_r,16*8);
figure; % 绘制目标C的距离、方位相位信息
% 图6.10(a)，(b)，实际上按书上描述的意思应该在Srcmc_t上进行分析，即随方位向时间采样点的分析
subplot(121);
plot((angle((c_rd_r_interp))));xlabel('距离向（采样点）');ylabel('相位（弧度）');title('(a)距离剖面相位');
subplot(122);
plot(unwrap(angle(c_rd_a.')));xlabel('方位向（采样点）');ylabel('相位（弧度）');title('(b)解绕后的方位剖面相位');
suptitle('距离徙动校正后目标C的距离和方位相位');

%3 绘制目标C的频谱，并对其升采样分析
target_c = s_ac(169-7:169+8,139-8:139+7);
target_C = fftshift(fft2(target_c));
[M_C, N_C] = size(target_C);
target_C_padr = [target_C(:,1:14),zeros(size(target_c,1),256-size(target_c,2)),target_C(:,15:end)];    % 距离向补零
target_C_padra = [target_C_padr(1:4,:);zeros(256-size(target_C_padr,1),size(target_C_padr,2));target_C_padr(5:end,:)];   % 方位向补0
target_c_interp = ifft2(circshift(target_C_padra, [-floor(M_C/2), -floor(N_C/2)])); % 注意这里不能直接用ifftshift，否则会发生补零错误
figure; % 绘制以目标C为中心16*16切片的频谱图
colormap('gray');
imagesc((abs(target_C)));xlabel('距离频率（采样点）');ylabel('方位频率（采样点）');title('目标C的频谱图');

figure; % 绘制目标C压缩后时域幅度特性
subplot(121);imagesc(abs(target_c_interp));xlabel('距离向（采样点）');ylabel('方位向（采样点）');title('(a)放大后的目标C');
subplot(122);contour(abs(target_c_interp),64);set(gca, 'YDir', 'reverse');xlabel('距离向（采样点）');ylabel('方位向（采样点）');title('(b)放大后目标C的轮廓图');
suptitle('目标C压缩后时域幅度特性');

%4 绘制升采样后目标C的距离方位包络剖面和相位
figure;
subplot(221);plot(20*log10(abs(target_c_interp(128,:))));xlabel('距离向（采样点）');ylabel('幅度');title('(a)距离向包络剖面图');
subplot(222);plot(20*log10(abs(target_c_interp(:,129))));xlabel('方位向（采样点）');ylabel('幅度');title('(b)方位向包络剖面图');
subplot(223);plot((angle(target_c_interp(128,:))));xlabel('距离向（采样点）');ylabel('相位（弧度）');title('(c)距离向相位剖面图');
subplot(224);plot((angle(target_c_interp(:,129))));xlabel('方位向（采样点）');ylabel('相位（弧度）');title('(d)方位向相位剖面图');
suptitle('目标C的包络相位特性');

