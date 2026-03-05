clc;clear all;
tic;
addpath('./utils/');
params.B = 3;
params.Omega = 10*2*pi/60;
params.R = 87;
params.rho = 1.225;
params.V0 = 9;
params.RootOffset = 0.01;
params.R_hub = 1.56;
params.maxIter = 100;
params.tolerance = 1e-5;
params.relax = 0.4;
params.time = 0:0.01:6;
params.N_harmonics = 4;
params.chi = 0;
params.rPowerMax = 8;
params.plot = true;
params.monitor = true;
params.useParallel = true;
params.useSparse = true;
params.pitch = 0;          % 桨距角(度), 正方向为绕z轴负向旋转
params.gamma = 0;          % 偏航角(度): 风轮旋转轴在水平面内偏离来流方向的角度
                           %   正方向: 从上方俯视，风轮轴顺时针偏离来流(即来流从左侧进入旋转平面)
params.yita  = 0;          % 仰角(度): 风轮旋转轴与水平面之间的夹角(即主轴上仰/下俯角度)
                           %   正方向: 旋转轴朝上倾斜(机头仰起), 来流从下方斜入旋转平面
params.hub_offset = 0;     % 轮毂偏置长度L(m): 旋转中心到叶片根部的额外距离
                           %   整机旋转平面直径变为 2*(R+L)
                           %   叶素旋转半径 r_rot = r + L，切向速度额外增加 Omega*L
params.cone = 0;           % 锥角(度): 叶片相对旋转平面的折转角
                           %   正方向: 顺气流方向折转(叶尖向下风侧偏转)
                           %   锥角仅改变来流速度在叶片局部坐标系中的轴向/切向投影分量
                           %   GDW尾流理论保持平面旋转假设不变
% 基于BEM初始化 + GDW三维非定常气动力计算
% 坐标系: x轴-旋转轴正向, z轴-叶片展向, y轴-叶片弦向

% ===== 必要输入检查 =====
requiredFields = {'Omega','B','R'};
for k = 1:numel(requiredFields)
    if ~isfield(params, requiredFields{k})
        error('缺少必要输入参数: %s', requiredFields{k});
    end
end

% ===== 读取输入与默认参数 =====
Omega = params.Omega;                  % rad/s
B = params.B;                          % 叶片数
R = params.R;                          % 叶片长度
rho = getField(params,'rho',1.225);
V0 = getField(params,'V0',9);          % 来流速度(沿x轴)
RootOffset = getField(params,'RootOffset',0.05); % 叶根偏移
R_hub = getField(params,'R_hub',1.56); % 叶轮半径(用于Prandtl损失计算)

% ===== 轮毂偏置长度 L =====
% 旋转中心到叶片根部的额外距离，旋转平面直径变为 2*(R+L)
% 叶素实际旋转半径: r_rot(i) = r(i) + L，切向速度额外增加 Omega*L
hub_L = getField(params,'hub_offset',0);  % 轮毂偏置长度L(m), 默认0

% ===== 锥角 =====
% 正方向: 顺气流方向折转(叶尖向下风侧)
% 来流 V0_ax 在叶片局部坐标系中:
%   轴向(垂直叶片旋转面)投影: V0_ax * cos(cone) — 保持驱动诱导速度
%   切向(沿叶片展向平面内)投影: V0_ax * sin(cone) — 顺锥角时减小有效切速
% GDW保持平面旋转假设，尾流矩阵不变
cone_deg = getField(params,'cone',0);   % 锥角(度), 默认0
cone_rad = deg2rad(cone_deg);

% ===== 来流速度分解：桨距角 + 偏航角 + 仰角 =====
%
% --- 桨距角 (pitch) ---
% pitch 正方向定义: 绕z轴负向旋转(即前缘向顺桨方向)
pitch_deg = getField(params,'pitch',0);   % 桨距角(度), 默认0度
pitch_rad = deg2rad(pitch_deg);
%
% --- 偏航角 (gamma) ---
% gamma 正方向: 从上方俯视，风轮轴顺时针偏离来流方向
%   => 来流在旋转平面内产生侧向分量 V0_lat = V0 * cos(yita_rad) * sin(gamma_rad)
%   => 轴向分量因偏航减小：参与计算的轴向分量包含 cos(gamma)
gamma_deg = getField(params,'gamma',0);   % 偏航角(度), 默认0度
gamma_rad = deg2rad(gamma_deg);
%
% --- 仰角 (yita) ---
% yita 正方向: 旋转轴朝上倾斜(机头仰起)
%   => 来流在旋转平面内产生俯仰面内分量 V0_elev = V0 * sin(yita_rad)
%   => 轴向分量因仰角减小：参与计算的轴向分量包含 cos(yita)
yita_deg = getField(params,'yita',0);     % 仰角(度), 默认0度
yita_rad = deg2rad(yita_deg);
%
% --- 速度矢量坐标变换 ---
% 来流速度 V0 沿惯性坐标系x轴(与地面平行、与未偏航/未倾斜的旋转轴平行)
% 经偏航(绕竖直轴旋转gamma)和仰角(绕横轴旋转yita)后，
% 在叶片旋转坐标系(x=旋转轴, y/z=旋转平面)的分解为：
%   V0_ax  = V0 * cos(yita) * cos(gamma)   -- 沿旋转轴方向(轴向)
%   V0_lat = V0 * cos(yita) * sin(gamma)   -- 旋转平面内侧向(偏航产生)
%   V0_elev= V0 * sin(yita)                -- 旋转平面内仰俯向(仰角产生)
% 旋转平面内总横向速度(面内分量合矢量):
%   V0_inplane = sqrt(V0_lat^2 + V0_elev^2)
% 再叠加桨距角对轴向/切向分量的投影修正:
%   V0_ax_eff = V0_ax * cos(pitch) + V0_inplane * sin(pitch) （桨距角进一步修正轴向）
%   -- 此处保持原桨距角逻辑，仅在偏斜基础上叠加；
%   -- 工程上桨距角与偏航/仰角为独立控制，本实现中桨距角仅作用于迎角(theta+pitch)，
%      来流分解由偏航+仰角完成，故 V0_ax/V0_ip 仅受偏航+仰角影响，
%      桨距角效果通过 theta_eff = theta + pitch_deg 体现在攻角计算中。
V0_ax  = V0 * cos(yita_rad) * cos(gamma_rad);  % 轴向有效来流分量
V0_lat = V0 * cos(yita_rad) * sin(gamma_rad);  % 面内侧向来流分量(偏航产生)
V0_elev= V0 * sin(yita_rad);                   % 面内仰俯向来流分量(仰角产生)
% 面内来流合矢量大小(用于advance ratio mu)
V0_inplane = sqrt(V0_lat^2 + V0_elev^2);
%
% 桨距角效果：作为附加安装角叠加到几何扭角，使有效扭角增加pitch_deg
% 在攻角计算中使用 theta_eff(i) = theta(i) + pitch_deg
% (正pitch减小攻角，即顺桨方向)

% ===== 偏斜入流参数 (用于BEM偏斜修正 Eq.17-19 和 GDW chi) =====
% 根据 Eq.19 (Burton et al. 2001): chi ≈ (0.6*a_avg + 1) * gamma_eff
% 其中 gamma_eff 为等效偏斜角 = atan2(V0_inplane, V0_ax)
gamma_eff_rad = atan2(V0_inplane, V0_ax + eps);   % 等效来流偏斜角(相对旋转轴)
% chi_eff 将在 BEM 初始化后根据平均诱导因子由 Eq.19 正式计算

% ===== 叶片结构模型与翼型数据 =====
[beamModel, r, c, theta, thickness, dr] = resolve_blade_geometry(params);

% 可调参数
maxIter = getField(params,'maxIter',100);
tolerance = getField(params,'tolerance',1e-5);
relax = getField(params,'relax',0.4);

% 时间步参数
if isfield(params,'time')
    t = params.time(:);
    dt = mean(diff(t));
else
    dt = getField(params,'dt',0.01);
    tEnd = getField(params,'tEnd',1.0);
    t = (0:dt:tEnd).';
end
Nt = numel(t);

% GDW参数
N_h = getField(params,'N_harmonics',4); % 
chi = getField(params,'chi',0); % 尾迹偏斜角
rPowerMax = getField(params,'rPowerMax',8);
S_h = getField(params,'S_h',default_S_h(N_h,rPowerMax));

% ===== 衍生量 =====
Nstations = numel(r);
rhat = r / R;
r_safe = max(r, eps);
sigma = (B * c) ./ (2 * pi * r_safe);

% 叶素实际旋转半径 (含轮毂偏置L)
r_rot = r + hub_L;                      % r_rot(i) = r(i) + L
Omega_r = Omega * r_rot;                % 切向速度基础分量 = Omega * (r + L)

% 锥角引起的来流速度在叶片局部坐标系中的分量修正
% V0_ax_cone: 轴向来流在叶片旋转面法向的投影 (cos因子，减小轴向入流)
% V0_ax_tan:  轴向来流沿叶展方向在旋转面内的投影 (sin因子，顺锥角时沿-切向叠加)
%   正锥角(叶尖顺气流折转) => V0_ax投影到切向为负(减小有效切速)
V0_ax_cone = V0_ax * cos(cone_rad);    % 叶素有效轴向来流 (锥角修正后)
V0_ax_tan  = -V0_ax * sin(cone_rad);   % 锥角引起的切向来流附加量(顺锥为负)

rootCutoff = RootOffset * R;
activeIdx = find(r >= rootCutoff);

% 并行与稀疏设置
useParallel = getField(params,'useParallel',false);
useSparse = getField(params,'useSparse',false);
if useParallel && license('test','Distrib_Computing_Toolbox')
    if isempty(gcp('nocreate'))
        parpool('threads');
    end
else
    useParallel = false;
end

% ===== 1) BEM初始化（稳态） =====
% 传入原始展向坐标 r（叶尖/叶根损失基于 r 和 R，不受 hub_L 影响）
% hub_L 作为额外参数传入，BEM内部切向速度改为 Omega*(r+hub_L)
% 轴向来流使用锥角修正后的 V0_ax_cone
[bem] = bem_initialize(r, c, theta, thickness, B, R, R_hub, RootOffset, V0_ax_cone, Omega, maxIter, tolerance, relax, hub_L);

% BEM初始化完成后，用平均轴向诱导因子更新等效尾迹偏斜角 chi_eff (Eq.19)
a_avg = mean(bem.a(activeIdx));
chi_eff = (0.6 * a_avg + 1) * gamma_eff_rad;  % Eq.19: chi ≈ (0.6a+1)*gamma_eff
fprintf('等效偏斜角 gamma_eff=%.2f deg, chi_eff=%.2f deg\n', ...
    rad2deg(gamma_eff_rad), rad2deg(chi_eff));

% ===== 2) GDW状态初始化 =====
mode = build_modes(N_h, S_h);
phiMat = precompute_shape_functions(mode, rhat);

% 初始化诱导速度系数: 仅用h=0拟合轴向诱导速度
u_ind0 = bem.a .* V0_ax_cone; % 轴向诱导速度(正向与V0_ax_cone相反)
uhat0 = u_ind0 / (Omega * R + eps);
alpha = cell(N_h+1,1);
beta = cell(N_h+1,1);
for h = 0:N_h
    Sh = S_h(h+1);
    alpha{h+1} = zeros(Sh,1);
    if h > 0
        beta{h+1} = zeros(Sh,1);
    else
        beta{h+1} = [];
    end
end
% h=0最小二乘拟合
Phi0 = phiMat{1};
alpha{1} = (Phi0' * Phi0 + 1e-8*eye(size(Phi0,2))) \ (Phi0' * uhat0);

% ===== 3) 预计算矩阵 =====
% GDW矩阵使用等效尾迹偏斜角 chi_eff (由偏航角+仰角综合计算，Eq.19)
% 返回全局矩阵（包含所有跨谐波耦合，Eq.74-76）
% offC(h+1): 谐波h的余弦状态在全局向量的起始偏移(0-based)
% offS(h)  : 谐波h(h>=1)的正弦状态在全局向量的起始偏移(0-based)
[Mdiag_c, Mdiag_s, LtildeC, LtildeS, offC, offS] = precompute_gdw_matrices(mode, chi_eff);
N_c = numel(Mdiag_c);   % 余弦状态总数
N_s = numel(Mdiag_s);   % 正弦状态总数

% 将 cell 数组格式的初始状态打包为全局向量
alphaVec = zeros(N_c, 1);
betaVec  = zeros(N_s, 1);
for h = 0:N_h
    Sh = S_h(h+1);
    if Sh == 0; continue; end
    alphaVec(offC(h+1)+(1:Sh)) = alpha{h+1};
end
for h = 1:N_h
    Sh = S_h(h+1);
    if Sh == 0; continue; end
    betaVec(offS(h)+(1:Sh)) = beta{h+1};
end

% ===== 4) 预分配输出 =====
L_sec = zeros(Nstations, B, Nt);
D_sec = zeros(Nstations, B, Nt);
M_sec = zeros(Nstations, B, Nt);
T_sec = zeros(Nstations, B, Nt);
Ft_sec = zeros(Nstations, B, Nt);
Alpha_sec = zeros(Nstations, B, Nt);
Phi_sec = zeros(Nstations, B, Nt);
Vrel_sec = zeros(Nstations, B, Nt);

Thrust = zeros(Nt,1);
Torque = zeros(Nt,1);
Power = zeros(Nt,1);
% 动态监控设置
monitor = getField(params,'monitor',false);
monitorEvery = getField(params,'monitorEvery',max(1,round(Nt/100)));
if monitor
    monFig = figure('Color','w','Name','动态监控','Position',[60 60 1400 760]);
    tl = tiledlayout(monFig, 2, 3, 'TileSpacing','compact','Padding','compact');

    % 上行: 各叶片x方向力(推力) + 总x方向力
    monColors = lines(B);
    axB = nexttile(tl, 1, [1 2]);   % 跨2列: 各叶片x方向推力
    hold(axB,'on');
    hBlade = gobjects(B,1);
    for q = 1:B
        hBlade(q) = plot(axB, t(1), 0, '-', 'Color', monColors(q,:), 'LineWidth', 1.5);
    end
    grid(axB,'on');
    xlabel(axB,'时间 (s)'); ylabel(axB,'x方向力 (N)');
    title(axB,'各叶片 x 方向力（推力）');
    bladeLegNames = arrayfun(@(q) sprintf('叶片%d',q), 1:B, 'UniformOutput', false);
    legend(axB, hBlade, bladeLegNames, 'Location','best','FontSize',8);

    axTx = nexttile(tl, 3);         % 总x方向力
    hTx = plot(axTx, t(1), 0, 'k-', 'LineWidth', 1.8);
    grid(axTx,'on');
    xlabel(axTx,'时间 (s)'); ylabel(axTx,'x方向力 (N)');
    title(axTx,'总 x 方向力（推力）');

    % 下行: 总推力 / 总扭矩 / 总功率
    ax1 = nexttile(tl); hT = plot(ax1, t(1), 0, 'b-','LineWidth',1.4);
    grid(ax1,'on'); xlabel(ax1,'时间(s)'); ylabel(ax1,'推力(N)'); title(ax1,'总推力');
    ax2 = nexttile(tl); hL = plot(ax2, t(1), 0, 'b-','LineWidth',1.4);
    grid(ax2,'on'); xlabel(ax2,'时间(s)'); ylabel(ax2,'升力(N)'); title(ax2,'总扭矩');
    ax4 = nexttile(tl); hM = plot(ax4, t(1), 0, 'b-','LineWidth',1.4);
    grid(ax4,'on'); xlabel(ax4,'时间(s)'); ylabel(ax4,'力矩(N·m)'); title(ax4,'总功率');
end

% ABM4历史（全局向量形式，与全局矩阵对应）
fHistA = zeros(N_c, 4);   % 余弦状态历史
fHistB = zeros(N_s, 4);   % 正弦状态历史

% 单叶片载荷时间序列预分配 (维度: B × Nt)
% 记录每片叶片在每个时间步的展向积分载荷及方位角，用于方位角-载荷分析
blade_Thrust = zeros(B, Nt);   % 单叶片推力(N)
blade_Torque = zeros(B, Nt);   % 单叶片扭矩(N·m)
blade_psi    = zeros(B, Nt);   % 单叶片方位角(rad)

% ===== 5) 时间推进 =====
dtHat = Omega * dt;
it = 1;
while true
    while it <= Nt
    psi_blade = Omega * t(it) + (0:B-1) * (2*pi/B);

    % ===== 偏斜入流修正系数 (Eq.17: a_skew = a*[1 + K*(r/R)*tan(chi/2)*cos(psi)]) =====
    % psi 定义: 0 为最下风侧位置(Eq.17 注释); 此处 psi_blade 从任意方位角零点计，
    % 对于风力机，通常令最下风侧(6点钟方向)为 psi=0，但本代码 psi 是从初始位置算起的方位角，
    % 物理上偏斜修正的方位角依赖性通过 cos(psi_blade) 体现，方向一致性已满足。
    skew_K = (15*pi/32) * tan(chi_eff/2);  % Eq.17 中的偏斜常数 K
    % skew_factor(i,q) = 1 + K*(r(i)/R)*cos(psi_blade(q))
    % 用于在叶素循环中修正局部轴向诱导因子: a_local = a_BEM * skew_factor

    % 计算诱导速度场(轴向)：从全局向量 alphaVec/betaVec 解包各谐波分量
    uhat = zeros(Nstations, B);
    for h = 0:N_h
        Sh = S_h(h+1);
        if Sh == 0; continue; end
        Phi_h = phiMat{h+1};                          % Nstations × Sh
        alpha_h = alphaVec(offC(h+1)+(1:Sh));         % 当前谐波余弦状态
        if h == 0
            uhat = uhat + Phi_h * alpha_h;
        else
            cosh_psi = cos(h * psi_blade);
            sinh_psi = sin(h * psi_blade);
            beta_h = betaVec(offS(h)+(1:Sh));          % 当前谐波正弦状态
            for q = 1:B
                uhat(:,q) = uhat(:,q) + Phi_h * (alpha_h * cosh_psi(q) + beta_h * sinh_psi(q));
            end
        end
    end
    u_ind = uhat * (Omega * R);

    % 计算气动力并累积推力/阻力/力矩
    L_step = zeros(Nstations, B);
    D_step = zeros(Nstations, B);
    M_step = zeros(Nstations, B);
    Alpha_step = zeros(Nstations, B);
    Phi_step = zeros(Nstations, B);
    Vrel_step = zeros(Nstations, B);
    dT_all = zeros(Nstations, B);
    dFt_all = zeros(Nstations, B);

    if useParallel
        parfor q = 1:B
            Lq = zeros(Nstations,1);
            Dq = zeros(Nstations,1);
            Mq = zeros(Nstations,1);
            Aq = zeros(Nstations,1);
            Pq = zeros(Nstations,1);
            Vq = zeros(Nstations,1);
            dTq = zeros(Nstations,1);
            dFq = zeros(Nstations,1);
            psi_q = psi_blade(q);  % 第q片叶片当前方位角
            for ii = 1:numel(activeIdx)
                i = activeIdx(ii);
                % % --- 偏斜修正局部轴向诱导因子 (Eq.17/29) ---
                % % a_skew = a_BEM * [1 + K*(r/R)*cos(psi)]
                % skew_fac = 1 + skew_K * rhat(i) * cos(psi_q);
                % a_local = bem.a(i) * skew_fac;
                % % 轴向速度: 锥角修正后的轴向来流 - GDW诱导速度 - BEM偏斜诱导
                % Uax = V0_ax_cone - u_ind(i,q) - a_local * V0_ax_cone;

                % 轴向速度: 锥角修正后的轴向来流 - GDW诱导速度
                Uax = V0_ax_cone - u_ind(i,q);

                % 切向速度: Omega*(r+L)*(1+a') + 偏航/仰角面内分量 + 锥角切向附加量
                % V0_ax_tan = -V0_ax*sin(cone): 顺锥角时来流分量沿叶展向，减小有效切速
                Ut = Omega_r(i)*(1+bem.ap(i)) + V0_lat*cos(psi_q) + V0_elev + V0_ax_tan; % 使用BEM的a'

                Vrel = hypot(Uax, Ut);
                Phi_i = atan2(Uax, Ut);
                % 有效扭角 = 几何扭角 + 桨距角 (正pitch减小攻角，即顺桨)
                Alpha_i = Phi_i - deg2rad(theta(i) + pitch_deg);

                [Cl_i, Cd_i, Cm_i] = aeroInterp(thickness(i), Alpha_i);

                dL = 0.5 * rho * Vrel^2 * c(i) * Cl_i * dr(i);
                dD = 0.5 * rho * Vrel^2 * c(i) * Cd_i * dr(i);
                dM = 0.5 * rho * Vrel^2 * c(i)^2 * Cm_i * dr(i);

                dT = dL * cos(Phi_i) - dD * sin(Phi_i);
                dFt = -dL * sin(Phi_i) - dD * cos(Phi_i);

                Lq(i) = dL; 
                Dq(i) = dD; 
                Mq(i) = dM;
                Aq(i) = Alpha_i; 
                Pq(i) = Phi_i; 
                Vq(i) = Vrel;
                dTq(i) = dT; 
                dFq(i) = dFt;
            end
            L_step(:,q) = Lq;
            D_step(:,q) = Dq;
            M_step(:,q) = Mq;
            Alpha_step(:,q) = Aq;
            Phi_step(:,q) = Pq;
            Vrel_step(:,q) = Vq;
            dT_all(:,q) = dTq;
            dFt_all(:,q) = dFq;
        end
    else
        for q = 1:B
            psi_q = psi_blade(q);  % 第q片叶片当前方位角
            for ii = 1:numel(activeIdx)
                i = activeIdx(ii);
                % % --- 偏斜修正局部轴向诱导因子 (Eq.17/29) ---
                % skew_fac = 1 + skew_K * rhat(i) * cos(psi_q);
                % a_local = bem.a(i) * skew_fac;
                % % 轴向速度: 锥角修正后的轴向来流 - GDW诱导速度 - BEM偏斜诱导
                % Uax = V0_ax_cone - u_ind(i,q) - a_local * V0_ax_cone;
                % 轴向速度: 锥角修正后的轴向来流 - GDW诱导速度
                Uax = V0_ax_cone - u_ind(i,q);

                % 切向速度: Omega*(r+L)*(1+a') + 偏航/仰角面内分量 + 锥角切向附加量
                Ut = Omega_r(i)*(1+bem.ap(i)) + V0_lat*cos(psi_q) + V0_elev + V0_ax_tan; % 先忽略a'            
                
                Vrel = hypot(Uax, Ut);
                Phi_i = atan2(Uax, Ut);
                % 有效扭角 = 几何扭角 + 桨距角 (正pitch减小攻角，即顺桨)
                Alpha_i = Phi_i - deg2rad(theta(i) + pitch_deg);

                [Cl_i, Cd_i, Cm_i] = aeroInterp(thickness(i), Alpha_i);

                dL = 0.5 * rho * Vrel^2 * c(i) * Cl_i * dr(i);
                dD = 0.5 * rho * Vrel^2 * c(i) * Cd_i * dr(i);
                dM = 0.5 * rho * Vrel^2 * c(i)^2 * Cm_i * dr(i);

                dT = dL * cos(Phi_i) - dD * sin(Phi_i);
                dFt = -dL * sin(Phi_i) - dD * cos(Phi_i);

                L_step(i,q) = dL;
                D_step(i,q) = dD;
                M_step(i,q) = dM;
                Alpha_step(i,q) = Alpha_i;
                Phi_step(i,q) = Phi_i;
                Vrel_step(i,q) = Vrel;
                dT_all(i,q) = dT;
                dFt_all(i,q) = dFt;
            end
        end
    end

    if useSparse
        dT_all = sparse(dT_all);
        dFt_all = sparse(dFt_all);
    end

    L_sec(:,:,it) = L_step./dr;  % 单位长度升力
    D_sec(:,:,it) = D_step./dr;  % 单位长度阻力
    M_sec(:,:,it) = M_step./dr;  % 单位长度俯仰力矩
    T_sec(:,:,it) = dT_all./dr;  % 单位长度推力
    Ft_sec(:,:,it) = dFt_all./dr;% 单位长度切向力
    Alpha_sec(:,:,it) = Alpha_step;
    Phi_sec(:,:,it) = Phi_step;
    Vrel_sec(:,:,it) = Vrel_step;

    Thrust(it) = sum(dT_all,'all');
    Torque(it) = sum(dFt_all .* r, 'all');
    Power(it) = Torque(it)*Omega;
    % ===== 单叶片展向积分载荷（直接按叶片索引累加，不做平均） =====
    for q = 1:B
        blade_Thrust(q,it) = sum(dT_all(:,q));
        blade_Torque(q,it) = sum(dFt_all(:,q) .* r);
        blade_psi(q,it)    = psi_blade(q);        % 当前方位角(rad)
    end

    if monitor && (mod(it-1, monitorEvery) == 0 || it == Nt)
        % 各叶片x方向力（推力）
        for q = 1:B
            set(hBlade(q), 'XData', t(1:it), 'YData', blade_Thrust(q,1:it));
        end
        % 总x方向力
        set(hTx, 'XData', t(1:it), 'YData', Thrust(1:it));
        % 总推力 / 总扭矩 / 总功率
        set(hT, 'XData', t(1:it), 'YData', Thrust(1:it));
        set(hL, 'XData', t(1:it), 'YData', Torque(1:it));
        set(hM, 'XData', t(1:it), 'YData', Power(1:it));
        drawnow limitrate;
    end

    % ===== 计算GDW压力系数tau =====
    [tauC, tauS] = compute_tau(dT_all, B, Omega, R, rho, psi_blade, mode, phiMat);

    % ===== 计算流动参数 =====
    % 面内速度比 mu: 使用偏航+仰角产生的面内合矢量 V0_inplane (Eq.86-87)
    % (若用户同时指定了 params.U_inplane，则在此基础上叠加)
    mu = V0_inplane / (Omega * R + eps);
    if isfield(params,'U_inplane')
        U_in = params.U_inplane;
        mu = (V0_inplane + norm(U_in)) / (Omega * R + eps);
    end
    % 轴向来流无量纲参数使用等效尾迹偏斜角 chi_eff (Eq.89)
    % 轴向来流取锥角修正后的分量 V0_ax_cone
    lambda_f = (V0_ax_cone * cos(chi_eff)) / (Omega * R + eps);
    alpha10 = alphaVec(offC(1)+1); % h=0, j=1（全局余弦向量第1个元素）
    lambda_m = sqrt(3) * alpha10;
    lambda = lambda_f + lambda_m;
    VT = sqrt(mu^2 + lambda^2);
    V  = (mu^2 + (lambda + lambda_m) * lambda) / (sqrt(mu^2 + lambda^2) + eps);

    % ===== 将 cell 形式的 tau 打包为全局向量 =====
    tauCvec = zeros(N_c, 1);
    tauSvec = zeros(N_s, 1);
    for h = 0:N_h
        Sh = S_h(h+1);
        if Sh == 0; continue; end
        tauCvec(offC(h+1)+(1:Sh)) = tauC{h+1};
    end
    for h = 1:N_h
        Sh = S_h(h+1);
        if Sh == 0; continue; end
        tauSvec(offS(h)+(1:Sh)) = tauS{h+1};
    end

    % ===== 组装全局 V 向量 =====
    % Eq.84: (m,n)=(0,1) 处用 VT，其余用 V
    VvecC = V * ones(N_c, 1);
    VvecC(offC(1)+1) = VT;           % h=0, j=1 对应 (m,n)=(0,1)
    VvecS = V * ones(N_s, 1);        % 正弦方程无 (m,n)=(0,1) 项

    % ===== GDW动态更新（全局矩阵统一求解，含跨谐波耦合） =====
    dAlphaVec = gdw_rhs(alphaVec, tauCvec, Mdiag_c, LtildeC, VvecC);
    [alphaVec, fHistA] = abm4_update(alphaVec, dAlphaVec, fHistA, dtHat, it);

    if N_s > 0
        dBetaVec = gdw_rhs(betaVec, tauSvec, Mdiag_s, LtildeS, VvecS);
        [betaVec, fHistB] = abm4_update(betaVec, dBetaVec, fHistB, dtHat, it);
    end
        it = it + 1;
    end

    % ===== 交互式继续仿真 =====
    cont = input('当前时间步已完成，是否继续仿真? (Y/n) ','s');
    if isempty(cont)
        cont = 'y';
    end
    if lower(cont(1)) ~= 'y'
        break;
    end

    extraDefault = getField(params,'extendDefault',1.0);
    extra = input(sprintf('请输入增加的仿真时长(s) [默认: %.3g]: ', extraDefault));
    if isempty(extra)
        extra = extraDefault;
    end
    if ~isnumeric(extra) || ~isscalar(extra) || ~isfinite(extra) || extra <= 0
        fprintf('输入无效，结束仿真。\n');
        break;
    end

    nExtra = max(1, round(extra / dt));
    tExtra = t(end) + (1:nExtra).' * dt;
    t = [t; tExtra];
    Nt_new = Nt + nExtra;

    L_sec(:,:,Nt+1:Nt_new) = 0;
    D_sec(:,:,Nt+1:Nt_new) = 0;
    M_sec(:,:,Nt+1:Nt_new) = 0;
    Alpha_sec(:,:,Nt+1:Nt_new) = 0;
    Phi_sec(:,:,Nt+1:Nt_new) = 0;
    Vrel_sec(:,:,Nt+1:Nt_new) = 0;

    Thrust(Nt+1:Nt_new,1) = 0;
    Torque(Nt+1:Nt_new,1) = 0;

    blade_Thrust(:,Nt+1:Nt_new) = 0;
    blade_Torque(:,Nt+1:Nt_new) = 0;
    blade_psi(:,Nt+1:Nt_new)    = 0;

    Nt = Nt_new;
end

% ===== 汇总BEM结果 ====

bemTotals = compute_bem_totals(r, c, theta, thickness, dr, B, RootOffset, V0, Omega, rho, bem.a, bem.ap);

% ===== 输出 =====
results.time = t;
results.section.L = L_sec;
results.section.D = D_sec;
results.section.M = M_sec;
results.section.T = T_sec;
results.section.Ft = Ft_sec;
results.section.alpha = Alpha_sec;
results.section.phi = Phi_sec;
results.section.Vrel = Vrel_sec;
% 单叶片展向积分载荷时间序列 (维度: B × Nt)
% 每行对应一片叶片，每列对应一个时间步，可直接按方位角 blade_psi 分析周期载荷
results.blade.Thrust = blade_Thrust;  % 单叶片推力(N)
results.blade.Torque = blade_Torque;  % 单叶片扭矩(N·m)
results.blade.psi    = blade_psi;     % 单叶片方位角(rad)
results.total.Thrust = Thrust;
results.total.Torque = Torque;
results.bem = bem;
results.bemTotals = bemTotals;
results.gdw.alphaVec = alphaVec;   % 余弦状态全局向量（含所有谐波，终态）
results.gdw.betaVec  = betaVec;    % 正弦状态全局向量（含所有谐波，终态）
results.gdw.offC = offC;           % 各谐波在余弦全局向量中的偏移
results.gdw.offS = offS;           % 各谐波在正弦全局向量中的偏移
results.meta = struct('Omega',Omega,'B',B,'R',R);
results.meta.pitch_deg     = pitch_deg;        % 桨距角(度)
results.meta.pitch_rad     = pitch_rad;        % 桨距角(弧度)
results.meta.gamma_deg     = gamma_deg;        % 偏航角(度)
results.meta.gamma_rad     = gamma_rad;        % 偏航角(弧度)
results.meta.yita_deg      = yita_deg;         % 仰角(度)
results.meta.yita_rad      = yita_rad;         % 仰角(弧度)
results.meta.gamma_eff_rad = gamma_eff_rad;    % 等效来流偏斜角(弧度)
results.meta.chi_eff       = chi_eff;          % 等效尾迹偏斜角(弧度, Eq.19)
results.meta.hub_L         = hub_L;            % 轮毂偏置长度L(m)
results.meta.cone_deg      = cone_deg;         % 锥角(度)
results.meta.cone_rad      = cone_rad;         % 锥角(弧度)
results.meta.beamModel = beamModel;
results.meta.geometry = struct('r',r,'chord',c,'twist',theta,'thickness',thickness,'dr',dr);

toc;