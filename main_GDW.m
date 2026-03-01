clc;clear all;
addpath('./utils/');
params.B = 3;
params.Omega = 10*2*pi/60;
params.R = 87;
params.rho = 1.225;
params.V0 = 9;
params.RootOffset = 0.05;
params.R_hub = 1.56;
params.maxIter = 100;
params.tolerance = 1e-5;
params.relax = 0.4;
params.time = 0:0.01:6;
params.N_harmonics = 0;
params.chi = 0;
params.rPowerMax = 2;
params.plot = true;
params.monitor = true;
params.useParallel = true;
params.useSparse = true;
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
Omega_r = Omega * r;
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
[bem] = bem_initialize(r, c, theta, thickness, B, R, R_hub, RootOffset, V0, Omega, maxIter, tolerance, relax);

% ===== 2) GDW状态初始化 =====
mode = build_modes(N_h, S_h);
phiMat = precompute_shape_functions(mode, rhat);

% 初始化诱导速度系数: 仅用h=0拟合轴向诱导速度
u_ind0 = bem.a .* V0; % 轴向诱导速度(正向与V0相反)
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
[Mdiag, LtildeC, LtildeS] = precompute_gdw_matrices(mode, chi);

% ===== 4) 预分配输出 =====
L_sec = zeros(Nstations, B, Nt);
D_sec = zeros(Nstations, B, Nt);
M_sec = zeros(Nstations, B, Nt);
Alpha_sec = zeros(Nstations, B, Nt);
Phi_sec = zeros(Nstations, B, Nt);
Vrel_sec = zeros(Nstations, B, Nt);

Thrust = zeros(Nt,1);
Drag = zeros(Nt,1);
Torque = zeros(Nt,1);
TotalLift_ts = zeros(Nt,1);
TotalDrag_ts = zeros(Nt,1);
TotalMoment_ts = zeros(Nt,1);

% 动态监控设置
monitor = getField(params,'monitor',false);
monitorEvery = getField(params,'monitorEvery',max(1,round(Nt/100)));
if monitor
    monFig = figure('Color','w','Name','动态监控');
    tl = tiledlayout(monFig,2,2,'TileSpacing','compact','Padding','compact');
    ax1 = nexttile(tl); hT = plot(ax1, t(1), 0, 'b-','LineWidth',1.4); grid(ax1,'on'); xlabel(ax1,'时间(s)'); ylabel(ax1,'推力(N)'); title(ax1,'总推力');
    ax2 = nexttile(tl); hL = plot(ax2, t(1), 0, 'b-','LineWidth',1.4); grid(ax2,'on'); xlabel(ax2,'时间(s)'); ylabel(ax2,'升力(N)'); title(ax2,'总升力');
    ax3 = nexttile(tl); hD = plot(ax3, t(1), 0, 'b-','LineWidth',1.4); grid(ax3,'on'); xlabel(ax3,'时间(s)'); ylabel(ax3,'阻力(N)'); title(ax3,'总阻力');
    ax4 = nexttile(tl); hM = plot(ax4, t(1), 0, 'b-','LineWidth',1.4); grid(ax4,'on'); xlabel(ax4,'时间(s)'); ylabel(ax4,'力矩(N·m)'); title(ax4,'总力矩');
end

% ABM4历史
fHistA = cell(N_h+1,1); % 每个h对应 [f_n, f_{n-1}, f_{n-2}, f_{n-3}]
fHistB = cell(N_h+1,1);
for h = 0:N_h
    Sh = S_h(h+1);
    fHistA{h+1} = zeros(Sh,4);
    if h > 0
        fHistB{h+1} = zeros(Sh,4);
    else
        fHistB{h+1} = [];
    end
end

% ===== 5) 时间推进 =====
dtHat = Omega * dt;
it = 1;
while true
    while it <= Nt
    psi_blade = Omega * t(it) + (0:B-1) * (2*pi/B);

    % 计算诱导速度场(轴向)
    uhat = zeros(Nstations, B);
    for h = 0:N_h
        Sh = S_h(h+1);
        if Sh == 0
            continue;
        end
        Phi_h = phiMat{h+1}; % Nstations x Sh
        if h == 0
            uhat = uhat + Phi_h * alpha{h+1};
        else
            cosh = cos(h * psi_blade);
            sinh = sin(h * psi_blade);
            for q = 1:B
                uhat(:,q) = uhat(:,q) + Phi_h * (alpha{h+1} * cosh(q) + beta{h+1} * sinh(q));
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
            for ii = 1:numel(activeIdx)
                i = activeIdx(ii);
                Uax = V0 - u_ind(i,q);

                Ut = Omega_r(i)*(1+bem.ap(i)); % 使用BEM的a'
                % Phi_i = atan2(Uax, Ut);
                % Alpha_i = Phi_i - deg2rad(theta(i));
                % [Cl_i, Cd_i, ~] = aeroInterp(thickness(i), Alpha_i);

                % F_tip = (2/pi) * real(acos(exp(-(B*(R - r(i)))/(2*r_safe(i)*abs(sin(Phi_i))+eps))));
                % F_root = (2/pi) * real(acos(exp(-B/2*(r(i)-R_hub)/(R_hub*abs(sin(Phi_i))+eps))));
                % F = max(F_tip*F_root, 1e-6);

                % sigma_i = sigma(i);
                % Kp = (4 * F * sin(Phi_i) * cos(Phi_i)) / (sigma_i * (Cl_i * sin(Phi_i) - Cd_i * cos(Phi_i) + eps));
                % ap = 1 / (Kp - 1 + eps);
                % Ut = Omega_r(i) * (1 + ap);

                Vrel = hypot(Uax, Ut);
                Phi_i = atan2(Uax, Ut);
                Alpha_i = Phi_i - deg2rad(theta(i));

                [Cl_i, Cd_i, Cm_i] = aeroInterp(thickness(i), Alpha_i);

                dL = 0.5 * rho * Vrel^2 * c(i) * Cl_i * dr(i);
                dD = 0.5 * rho * Vrel^2 * c(i) * Cd_i * dr(i);
                dM = 0.5 * rho * Vrel^2 * c(i)^2 * Cm_i * dr(i);

                dT = dL * cos(Phi_i) + dD * sin(Phi_i);
                dFt = dL * sin(Phi_i) - dD * cos(Phi_i);

                Lq(i) = dL; Dq(i) = dD; Mq(i) = dM;
                Aq(i) = Alpha_i; Pq(i) = Phi_i; Vq(i) = Vrel;
                dTq(i) = dT; dFq(i) = dFt;
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
            for ii = 1:numel(activeIdx)
                i = activeIdx(ii);
                Uax = V0 - u_ind(i,q);

                Ut = Omega_r(i)*(1+bem.ap(i)); % 先忽略a'
                % Phi_i = atan2(Uax, Ut);
                % Alpha_i = Phi_i - deg2rad(theta(i));
                % [Cl_i, Cd_i, ~] = aeroInterp(thickness(i), Alpha_i);
                % 
                % F_tip = (2/pi) * real(acos(exp(-(B*(R - r(i)))/(2*r_safe(i)*abs(sin(Phi_i))+eps))));
                % F_root = (2/pi) * real(acos(exp(-B/2*(r(i)-R_hub)/(R_hub*abs(sin(Phi_i))+eps))));
                % F = max(F_tip*F_root, 1e-6);
                % 
                % sigma_i = sigma(i);
                % Kp = (4 * F * sin(Phi_i) * cos(Phi_i)) / (sigma_i * (Cl_i * sin(Phi_i) - Cd_i * cos(Phi_i) + eps));
                % ap = 1 / (Kp - 1 + eps);
                % Ut = Omega_r(i) * (1 + ap);
                
                Vrel = hypot(Uax, Ut);
                Phi_i = atan2(Uax, Ut);
                Alpha_i = Phi_i - deg2rad(theta(i));

                [Cl_i, Cd_i, Cm_i] = aeroInterp(thickness(i), Alpha_i);

                dL = 0.5 * rho * Vrel^2 * c(i) * Cl_i * dr(i);
                dD = 0.5 * rho * Vrel^2 * c(i) * Cd_i * dr(i);
                dM = 0.5 * rho * Vrel^2 * c(i)^2 * Cm_i * dr(i);

                dT = dL * cos(Phi_i) + dD * sin(Phi_i);
                dFt = dL * sin(Phi_i) - dD * cos(Phi_i);

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

    L_sec(:,:,it) = L_step;
    D_sec(:,:,it) = D_step;
    M_sec(:,:,it) = M_step;
    Alpha_sec(:,:,it) = Alpha_step;
    Phi_sec(:,:,it) = Phi_step;
    Vrel_sec(:,:,it) = Vrel_step;

    Thrust(it) = sum(dT_all,'all');
    Drag(it) = sum(dFt_all,'all');
    Torque(it) = sum(dFt_all .* r, 'all');
    TotalLift_ts(it) = sum(L_step,'all');
    TotalDrag_ts(it) = sum(D_step,'all');
    TotalMoment_ts(it) = sum(M_step,'all');

    if monitor && (mod(it-1, monitorEvery) == 0 || it == Nt)
        set(hT,'XData',t(1:it),'YData',Thrust(1:it));
        set(hL,'XData',t(1:it),'YData',TotalLift_ts(1:it));
        set(hD,'XData',t(1:it),'YData',TotalDrag_ts(1:it));
        set(hM,'XData',t(1:it),'YData',TotalMoment_ts(1:it));
        drawnow limitrate;
    end

    % ===== 计算GDW压力系数tau =====
    [tauC, tauS] = compute_tau(dT_all, B, Omega, R, rho, psi_blade, mode, phiMat);

    % ===== 计算流动参数 =====
    mu = 0; % 默认无横风
    if isfield(params,'U_inplane')
        U_in = params.U_inplane;
        mu = norm(U_in) / (Omega * R + eps);
    end
    lambda_f = (V0 * cos(chi)) / (Omega * R + eps);
    alpha10 = alpha{1}(1); % h=0, j=1
    lambda_m = sqrt(3) * alpha10;
    lambda = lambda_f + lambda_m;
    VT = sqrt(mu^2 + lambda^2);
    V = (mu^2 + (lambda + lambda_m) * lambda) / (sqrt(mu^2 + lambda^2) + eps);

    % ===== GDW动态更新 =====
    for h = 0:N_h
        Sh = S_h(h+1);
        if Sh == 0
            continue;
        end

        % 组装V向量
        Vvec = V * ones(Sh,1);
        if h == 0
            Vvec(1) = VT; % (m,n)=(0,1)
        end

        % 计算导数
        dAlpha = gdw_rhs(alpha{h+1}, tauC{h+1}, Mdiag{h+1}, LtildeC{h+1}, Vvec);

        % ABM4更新
        [alpha{h+1}, fHistA{h+1}] = abm4_update(alpha{h+1}, dAlpha, fHistA{h+1}, dtHat, it);

        if h > 0
            dBeta = gdw_rhs(beta{h+1}, tauS{h+1}, Mdiag{h+1}, LtildeS{h+1}, Vvec);
            [beta{h+1}, fHistB{h+1}] = abm4_update(beta{h+1}, dBeta, fHistB{h+1}, dtHat, it);
        end
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
    Drag(Nt+1:Nt_new,1) = 0;
    Torque(Nt+1:Nt_new,1) = 0;
    TotalLift_ts(Nt+1:Nt_new,1) = 0;
    TotalDrag_ts(Nt+1:Nt_new,1) = 0;
    TotalMoment_ts(Nt+1:Nt_new,1) = 0;

    Nt = Nt_new;
end

% ===== 汇总总升力/阻力/力矩 =====
totalLift = TotalLift_ts;
totalDrag = TotalDrag_ts;
totalMoment = TotalMoment_ts;

bemTotals = compute_bem_totals(r, c, theta, thickness, dr, B, RootOffset, V0, Omega, rho, bem.a, bem.ap);

% ===== 输出 =====
results.time = t;
results.section.L = L_sec;
results.section.D = D_sec;
results.section.M = M_sec;
results.section.alpha = Alpha_sec;
results.section.phi = Phi_sec;
results.section.Vrel = Vrel_sec;
results.total.Thrust = Thrust;
results.total.Tangential = Drag;
results.total.Torque = Torque;
results.total.Lift = totalLift;
results.total.Drag = totalDrag;
results.total.Moment = totalMoment;
results.bem = bem;
results.bemTotals = bemTotals;
results.gdw.alpha = alpha;
results.gdw.beta = beta;
results.meta = struct('Omega',Omega,'B',B,'R',R);
results.meta.beamModel = beamModel;
results.meta.geometry = struct('r',r,'chord',c,'twist',theta,'thickness',thickness,'dr',dr);

if getField(params,'plot',false)
    plotOpts = getField(params,'plotOpts',struct());
    plot_aero_comparison(results, plotOpts);
end

function S_h = default_S_h(N_h, rPowerMax)
S_h = zeros(N_h+1,1);
for h = 0:N_h
    S_h(h+1) = floor((rPowerMax - h)/2) + 1;
    S_h(h+1) = max(S_h(h+1), 0);
end
end

