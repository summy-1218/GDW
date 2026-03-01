function results = BEM(params)
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
params.N_harmonics = 4;
params.chi = 0;
params.rPowerMax = 8;
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
[phiMat, mode] = precompute_shape_functions(mode, rhat);

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

                Ut = Omega_r(i); % 先忽略a'
                Phi_i = atan2(Uax, Ut);
                Alpha_i = Phi_i - deg2rad(theta(i));
                [Cl_i, Cd_i, ~] = aeroInterp(thickness(i), Alpha_i);

                F_tip = (2/pi) * real(acos(exp(-(B*(R - r(i)))/(2*r_safe(i)*abs(sin(Phi_i))+eps))));
                F_root = (2/pi) * real(acos(exp(-B/2*(r(i)-R_hub)/(R_hub*abs(sin(Phi_i))+eps))));
                F = max(F_tip*F_root, 1e-6);

                sigma_i = sigma(i);
                Kp = (4 * F * sin(Phi_i) * cos(Phi_i)) / (sigma_i * (Cl_i * sin(Phi_i) - Cd_i * cos(Phi_i) + eps));
                ap = 1 / (Kp - 1 + eps);
                Ut = Omega_r(i) * (1 + ap);

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

                Ut = Omega_r(i); % 先忽略a'
                Phi_i = atan2(Uax, Ut);
                Alpha_i = Phi_i - deg2rad(theta(i));
                [Cl_i, Cd_i, ~] = aeroInterp(thickness(i), Alpha_i);

                F_tip = (2/pi) * real(acos(exp(-(B*(R - r(i)))/(2*r_safe(i)*abs(sin(Phi_i))+eps))));
                F_root = (2/pi) * real(acos(exp(-B/2*(r(i)-R_hub)/(R_hub*abs(sin(Phi_i))+eps))));
                F = max(F_tip*F_root, 1e-6);

                sigma_i = sigma(i);
                Kp = (4 * F * sin(Phi_i) * cos(Phi_i)) / (sigma_i * (Cl_i * sin(Phi_i) - Cd_i * cos(Phi_i) + eps));
                ap = 1 / (Kp - 1 + eps);
                Ut = Omega_r(i) * (1 + ap);
                
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

end

% 子函数
function [yNew, fHist] = abm4_update(y, fNow, fHist, dtHat, it)
% fHist列: [f_n, f_{n-1}, f_{n-2}, f_{n-3}]
if it < 4
    % 启动阶段: 显式欧拉
    yNew = y + dtHat * fNow;
else
    f_n = fHist(:,1);
    f_n1 = fHist(:,2);

    % Adams-Bashforth 4
    % 预测导数(这里用当前导数近似)
    fPred = fNow;

    % Adams-Moulton 4 修正
    yNew = y + dtHat/24 * (9*fPred + 19*fNow - 5*f_n + f_n1);
end

% 更新历史
fHist = [fNow, fHist(:,1:3)];
end

function [cl,cd,cm] = aeroInterp(thicknessi,alphai)
load cl_cd_cm.mat;
alpha = DU21(:,1);
thickness = [0.1,0.21,0.25,0.3,0.35,0.4,0.9,1.0];

% 规范化原始数据点位


cl_21 = pchip(DU21(:,1),DU21(:,2),alpha);
cd_21 = pchip(DU21(:,1),DU21(:,3),alpha);
cm_21 = pchip(DU21(:,1),DU21(:,4),alpha);


cl_10 = cl_21 ;
cd_10 = cd_21;
cm_10 = cm_21;

cl_25 = pchip(DU25(:,1),DU25(:,2),alpha);
cd_25 = pchip(DU25(:,1),DU25(:,3),alpha);
cm_25 = pchip(DU25(:,1),DU25(:,4),alpha);

cl_30 = pchip(DU30(:,1),DU30(:,2),alpha);
cd_30 = pchip(DU30(:,1),DU30(:,3),alpha);
cm_30 = pchip(DU30(:,1),DU30(:,4),alpha);

cl_35 = pchip(DU35(:,1),DU35(:,2),alpha);
cd_35 = pchip(DU35(:,1),DU35(:,3),alpha);
cm_35 = pchip(DU35(:,1),DU35(:,4),alpha);

cl_40 = pchip(DU40(:,1),DU40(:,2),alpha);
cd_40 = pchip(DU40(:,1),DU40(:,3),alpha);
cm_40 = pchip(DU40(:,1),DU40(:,4),alpha);

cl_90 = pchip(DU90(:,1),DU90(:,2),alpha);
cd_90 = pchip(DU90(:,1),DU90(:,3),alpha);
cm_90 = pchip(DU90(:,1),DU90(:,4),alpha);

cl_100 = cl_90;
cd_100 = cd_90;
cm_100 = cm_90;

[X, Y] = meshgrid(thickness,alpha);  % 转换为网格矩阵
z_cl = [cl_10,cl_21,cl_25,cl_30,cl_35,cl_40,cl_90,cl_100];
z_cd = [cd_10,cd_21,cd_25,cd_30,cd_35,cd_40,cd_90,cd_100];
z_cm = [cm_10,cm_21,cm_25,cm_30,cm_35,cm_40,cm_90,cm_100];

cl = interp2(X,Y,z_cl,thicknessi/100,rad2deg(alphai),'linear');
cd = interp2(X,Y,z_cd,thicknessi/100,rad2deg(alphai),'linear');
cm = interp2(X,Y,z_cm,thicknessi/100,rad2deg(alphai),'linear');

end

function [bem] = bem_initialize(r, c, theta, thickness, B, R, R_hub, RootOffset, V0, Omega, maxIter, tolerance, relax)
N = numel(r); % 径向站数
a = zeros(N,1); % 轴向诱导因子
ap = zeros(N,1); % 周向诱导因子
for i = 1:N
    ri = r(i);
    % 叶根截止点以内不计算诱导因子
    if r(i) < RootOffset*R
        continue;
    end
    if abs(V0) < 1e-8
        a(i) = 0;
        ap(i) = 0;
        continue;
    end
    % 计算局部运动速度,先留好接口
    Ve_op = 0;
    Ve_ip = 0;

    % 局部实度
    sigma = (B * c(i)) / (2 * pi * r(i));

    % 局部尖速比
    lambda = Omega * r(i) / (V0 + eps);

    % 为每个径向站初始化一个猜测值
    a(i) = 0.25*(2+pi*lambda*sigma-sqrt(4-4*pi*lambda*sigma+pi*lambda^2*sigma*(8*deg2rad(theta(i))+pi*sigma)));
    ap(i) = 0;

    for iter = 1:maxIter
        a_old = a(i);
        ap_old = ap(i);

        % 1. 计算入流角 Phi (弧度)
        V_local_axial = V0 * (1 - a(i))+Ve_op;
        V_local_tangential = Omega * ri * (1 + ap(i))+Ve_ip;
        V_local = sqrt(V_local_axial^2 + V_local_tangential^2); % 相对合速度

        if V_local_axial == 0 % 避免除零
            Phi_i = pi/2;
        else
            Phi_i = atan2(V_local_axial, V_local_tangential); % 使用atan2更稳健
        end

        % 2. 计算攻角 Alpha (度)
        Alpha_i = Phi_i - deg2rad(theta(i)); % 几何安装角theta已为度

        % 3. 插值获取当前攻角下的升力系数Cl和阻力系数Cd
        [Cl_i,Cd_i,~] = aeroInterp(thickness(i), Alpha_i);

        % 4. 计算推力系数
        CT = sigma*(1-a_old)^2*(Cl_i*cos(Phi_i)+Cd_i*sin(Phi_i))/(sin(Phi_i))^2;

        % 5. 计算叶尖损失因子和轮毂损失因子
        F_tip = (2/pi) * real(acos(exp(-(B*(R - ri))/(2*ri*abs(sin(Phi_i))))));
        F_root = (2/pi) * real(acos(exp(-B/2*(ri-R_hub)/(R_hub*sin(Phi_i)))));
        F = F_tip*F_root;
        if isnan(F) || F == 0
            F = 1e-6; % 避免除零或无效值
        end

        % 6. 计算新的诱导因子a_new和ap_new
         % 轴向诱导因子更新方程
        if CT > 0.96*F
            % 使用Glauert经验修正
            a_new = (18*F-20-3*sqrt(CT*(50-36*F)+12*F*(3*F-4)))/(36*F-50);
        else
            % 使用经典方法更新
            K = (4 * F_tip * (sin(Phi_i))^2) / (sigma * (Cl_i * cos(Phi_i) + Cd_i * sin(Phi_i)));
            a_new = 1 / (K + 1);
        end
        % 周向诱导因子更新方程 
        Kp = (4 * F_tip * sin(Phi_i) * cos(Phi_i)) / (sigma * (Cl_i * sin(Phi_i) - Cd_i * cos(Phi_i)));
        ap_new = 1 / (Kp - 1);

        % 使用松弛因子更新诱导因子
        a(i) = (1 - relax) * a_old + relax * a_new;
        ap(i) = (1 - relax) * ap_old + relax * ap_new;

        if (abs(a(i) - a_old) < tolerance) && (abs(ap(i) - ap_old) < tolerance)
            converged = true;
            break;
        end
    end
    if converged
        fprintf('径向站 r=%.2f m 迭代收敛，共 %d 次迭代。a=%.4f, a''=%.4f\n', ri, iter, a(i), ap(i));
    else
        fprintf('警告: 径向站 r=%.2f m 未在 %d 次迭代内收敛。a=%.4f, a''=%.4f\n', ri, maxIter, a(i), ap(i));
    end

end
bem.a = a;
bem.ap = ap;
bem.r = r;
end

function mode = build_modes(N_h, S_h)
% mode包含方位角谐波的最高数目Nh，第h次谐波对应的形函数数目Sh
mode = struct();
mode.N_h = N_h;
mode.S_h = S_h;
mode.jList = cell(N_h+1,1);
for h = 0:N_h
    Sh = S_h(h+1);
    j = (h+1):2:(2*Sh+h-1);
    mode.jList{h+1} = j(:);
end
end

function bemTotals = compute_bem_totals(r, c, theta, thickness, dr, B, RootOffset, V0, Omega, rho, a, ap)
N = numel(r);
L = 0; D = 0; M = 0; T = 0; Q = 0;
for i = 1:N
    if r(i) < RootOffset
        continue;
    end
    Vax = V0 * (1 - a(i));
    Vtan = Omega * r(i) * (1 + ap(i));
    Phi = atan2(Vax, Vtan);
    Alpha = Phi - deg2rad(theta(i));
    [Cl, Cd, Cm] = aeroInterp(thickness(i), Alpha);
    Vrel = hypot(Vax, Vtan);
    dL = 0.5 * rho * Vrel^2 * c(i) * Cl * dr(i);
    dD = 0.5 * rho * Vrel^2 * c(i) * Cd * dr(i);
    dM = 0.5 * rho * Vrel^2 * c(i)^2 * Cm * dr(i);
    dT = dL * cos(Phi) + dD * sin(Phi);
    dFt = dL * sin(Phi) - dD * cos(Phi);
    L = L + dL * B;
    D = D + dD * B;
    M = M + dM * B;
    T = T + dT * B;
    Q = Q + dFt * r(i) * B;
end
bemTotals.Lift = L;
bemTotals.Drag = D;
bemTotals.Moment = M;
bemTotals.Thrust = T;
bemTotals.Torque = Q;
end

function [tauC, tauS] = compute_tau(dT_all, B, Omega, R, rho, psi_blade, mode, phiMat)
N_h = mode.N_h;
A = pi * R^2;
scale = 1 / (rho * A * (Omega*R)^2 + eps);

% 归一化叶素推力
Fhat = dT_all * scale; % Nstations x B

% tauC/tauS按谐波存储
tauC = cell(N_h+1,1);
tauS = cell(N_h+1,1);

for h = 0:N_h
    Phi_h = phiMat{h+1};
    Sh = size(Phi_h,2);
    tauC{h+1} = zeros(Sh,1);
    if h > 0
        tauS{h+1} = zeros(Sh,1);
    else
        tauS{h+1} = [];
    end

    for nIdx = 1:Sh
        phi_n = Phi_h(:,nIdx);
        sumC = 0;
        sumS = 0;
        for q = 1:B
            Fq = Fhat(:,q);
            sumC = sumC + sum(Fq .* phi_n) * cos(h * psi_blade(q));
            if h > 0
                sumS = sumS + sum(Fq .* phi_n) * sin(h * psi_blade(q));
            end
        end
        if h == 0
            tauC{h+1}(nIdx) = (1/(2*pi)) * sumC;
        else
            tauC{h+1}(nIdx) = (1/pi) * sumC;
            tauS{h+1}(nIdx) = (1/pi) * sumS;
        end
    end
end
end

function df = double_factorial(n)
if n <= 0
    df = 1;
    return;
end
if mod(n,2)==0
    df = prod(2:2:n);
else
    df = prod(1:2:n);
end
end

function Gamma = Gamma_jn_hm(j,n,h,m)
% Gamma系数: 式(4-77)~(4-79)
if mod(h+m,2)==0
    Gamma = ((-1)^((n+j-2*h)/2)) / sqrt(H_jh(n,m)*H_jh(j,h)) * (2*sqrt((2*n+1)*(2*j+1))) / ((j+n)*(j+n+2)*((j-n)^2-1));
else
    if abs(j-n)==1
        Gamma = pi / sqrt(H_jh(n,m)*H_jh(j,h)) * sign(h-m) / sqrt((2*n+1)*(2*j+1));
    else
        Gamma = 0;
    end
end
end

function val = getField(s, name, defaultVal)
if isfield(s, name)
    val = s.(name);
else
    val = defaultVal;
end
end

function H = H_jh(j,h)
% 详见理论文档式（4-60）
H = double_factorial(j+h-1) * double_factorial(j-h-1) / (double_factorial(j+h) * double_factorial(j-h) + eps);
end

function plot_aero_comparison(results, opts)
if nargin < 2
    opts = struct();
end

t = results.time(:);
L = results.total.Lift(:);
D = results.total.Drag(:);
M = results.total.Moment(:);

bemL = results.bemTotals.Lift;
bemD = results.bemTotals.Drag;
bemM = results.bemTotals.Moment;

smoothWin = getField(opts,'smoothWindow',1);
if smoothWin > 1
    L = smoothdata(L,'movmean',smoothWin);
    D = smoothdata(D,'movmean',smoothWin);
    M = smoothdata(M,'movmean',smoothWin);
end

fig = figure('Color','w','Position',[100 100 1200 700]);

titleStr = getField(opts,'title','气动性能参数对比分析');

subplot(3,1,1);
plot(t, L, 'b-', 'LineWidth', 1.8); hold on;
plot([t(1) t(end)], [bemL bemL], 'r--', 'LineWidth', 1.6);
xlabel('时间(s)'); ylabel('力(N)');
legend('GDW结果','BEM理论','Location','best');
grid on; xlim([t(1) t(end)]);

subplot(3,1,2);
plot(t, D, 'b-', 'LineWidth', 1.8); hold on;
plot([t(1) t(end)], [bemD bemD], 'r--', 'LineWidth', 1.6);
xlabel('时间(s)'); ylabel('力(N)');
legend('GDW结果','BEM理论','Location','best');
grid on; xlim([t(1) t(end)]);

subplot(3,1,3);
plot(t, M, 'b-', 'LineWidth', 1.8); hold on;
plot([t(1) t(end)], [bemM bemM], 'r--', 'LineWidth', 1.6);
xlabel('时间(s)'); ylabel('力矩(N·m)');
legend('GDW结果','BEM理论','Location','best');
grid on; xlim([t(1) t(end)]);

sgtitle(titleStr);

% 导出
if getField(opts,'export',false)
    outDir = getField(opts,'outDir',pwd);
    baseName = getField(opts,'baseName','aero_comparison');
    formats = getField(opts,'formats',{'png','pdf'});
    dpi = getField(opts,'dpi',300);
    for k = 1:numel(formats)
        fmt = lower(formats{k});
        outFile = fullfile(outDir, [baseName '.' fmt]);
        if strcmp(fmt,'png')
            exportgraphics(fig, outFile, 'Resolution', dpi);
        elseif strcmp(fmt,'pdf')
            exportgraphics(fig, outFile, 'ContentType', 'vector');
        else
            exportgraphics(fig, outFile);
        end
    end
end
end

function [Mdiag, LtildeC, LtildeS] = precompute_gdw_matrices(mode, chi)
N_h = mode.N_h;
Mdiag = cell(N_h+1,1);
LtildeC = cell(N_h+1,1);
LtildeS = cell(N_h+1,1);
X = tan(abs(chi)/2);

for h = 0:N_h
    jList = mode.jList{h+1};
    Sh = numel(jList);

    % M矩阵(对角)
    Mdiag{h+1} = zeros(Sh,1);
    for a = 1:Sh
        n = jList(a);
        Mdiag{h+1}(a) = (2/pi) * H_jh(n,h);
    end

    % Ltilde矩阵
    Lc = zeros(Sh,Sh);
    Ls = zeros(Sh,Sh);
    for a = 1:Sh
        j = jList(a);
        for b = 1:Sh
            n = jList(b);
            if h == 0
                Gamma = Gamma_jn_hm(j,n,0,0);
                Lc(a,b) = (X^0) * Gamma; % m=0
            else
                Gamma = Gamma_jn_hm(j,n,h,h);
                Pi = min(h,h);
                Lc(a,b) = (X^0 + (-1)^Pi * X^(2*h)) * Gamma;
                Ls(a,b) = (X^0 - (-1)^Pi * X^(2*h)) * Gamma;
            end
        end
    end
    LtildeC{h+1} = Lc + 1e-8*eye(Sh);
    if h > 0
        LtildeS{h+1} = Ls + 1e-8*eye(Sh);
    else
        LtildeS{h+1} = [];
    end
end
end

function [phiMat, mode] = precompute_shape_functions(mode, rhat)
N_h = mode.N_h;
phiMat = cell(N_h+1,1);
for h = 0:N_h
    jList = mode.jList{h+1};
    Sh = numel(jList);
    Phi = zeros(numel(rhat), Sh);
    for idx = 1:Sh
        j = jList(idx);
        Phi(:,idx) = shape_function(h, j, rhat);
    end
    phiMat{h+1} = Phi;
end
end

function [beamModel, r, c, theta, thickness, dr] = resolve_blade_geometry(params)
beamModel = struct();
if isfield(params,'beamModel')
    beamModel = params.beamModel;
elseif isfield(params,'beamModelFile')
    S = load(params.beamModelFile);
    if isfield(S,'beamModel')
        beamModel = S.beamModel;
    else
        beamModel = S;
    end
elseif exist('beamModel.mat','file')
    S = load('beamModel.mat');
    if isfield(S,'beamModel')
        beamModel = S.beamModel;
    else
        beamModel = S;
    end
end

% 读取几何参数：优先params，其次beamModel
if isfield(params,'r')
    r = params.r(:);
elseif isfield(beamModel,'distance')
    r = beamModel.distance(:);
else
    error('缺少叶片剖面径向位置: 请提供 params.r 或 beamModel.distance');
end

if isfield(params,'chord')
    c = params.chord(:);
elseif isfield(beamModel,'chord')
    c = beamModel.chord(:);
else
    error('缺少弦长分布: 请提供 params.chord 或 beamModel.chord');
end

if isfield(params,'twist')
    theta = params.twist(:);
elseif isfield(beamModel,'twist')
    theta = beamModel.twist(:);
else
    error('缺少扭角分布: 请提供 params.twist 或 beamModel.twist');
end

if isfield(params,'thickness')
    thickness = params.thickness(:);
elseif isfield(beamModel,'Thickness')
    thickness = beamModel.Thickness(:);
elseif isfield(beamModel,'thickness')
    thickness = beamModel.thickness(:);
else
    error('缺少翼型厚度参数: 请提供 params.thickness 或 beamModel.Thickness');
end

% 统一厚度单位(若为比例则转为百分数)
if max(thickness) <= 1.5
    thickness = thickness * 100;
end

% 叶素宽度由相邻距离决定
if isfield(params,'dr')
    dr = params.dr(:);
elseif isfield(beamModel,'distance')
    dr = [diff(r); diff(r(end-1:end))];
else
    dr = [diff(r); diff(r(end-1:end))];
end

% 匹配性检查
n = numel(r);
if numel(c) ~= n || numel(theta) ~= n || numel(thickness) ~= n || numel(dr) ~= n
    error('剖面数据长度不一致: r(%d), chord(%d), twist(%d), thickness(%d), dr(%d)', n, numel(c), numel(theta), numel(thickness), numel(dr));
end

% 若params与beamModel同时存在则校验一致性
if isfield(params,'r') && isfield(beamModel,'distance')
    if max(abs(params.r(:) - beamModel.distance(:))) > 1e-6
        error('params.r 与 beamModel.distance 不一致');
    end
end
if isfield(params,'chord') && isfield(beamModel,'chord')
    if max(abs(params.chord(:) - beamModel.chord(:))) > 1e-6
        error('params.chord 与 beamModel.chord 不一致');
    end
end
if isfield(params,'twist') && isfield(beamModel,'twist')
    if max(abs(params.twist(:) - beamModel.twist(:))) > 1e-6
        error('params.twist 与 beamModel.twist 不一致');
    end
end
if isfield(params,'thickness') && (isfield(beamModel,'Thickness') || isfield(beamModel,'thickness'))
    bmThick = thickness;
    if isfield(beamModel,'Thickness')
        bmThick = beamModel.Thickness(:);
    elseif isfield(beamModel,'thickness')
        bmThick = beamModel.thickness(:);
    end
    if max(bmThick) <= 1.5
        bmThick = bmThick * 100;
    end
    if max(abs(params.thickness(:) - bmThick(:))) > 1e-6
        error('params.thickness 与 beamModel 厚度数据不一致');
    end
end
end

function phi = shape_function(h, j, rhat)
% 形函数: 式(4-59)
% phi_j^h(r) = sqrt((2j+1)H_j^h) * sum_{q=h,h+2,...}^{j-1} r^q * (-1)^((q-h)/2) * (j+q)!! / ((q-h)!! (q+h)!! (j-q-1)!!)
H = H_jh(j, h);
pref = sqrt((2*j+1) * H);
phi = zeros(size(rhat));
for q = h:2:(j-1)
    num = double_factorial(j+q);
    den = double_factorial(q-h) * double_factorial(q+h) * double_factorial(j-q-1);
    coef = (-1)^((q-h)/2) * num / den;
    phi = phi + coef * (rhat.^q);
end
phi = pref * phi;
end