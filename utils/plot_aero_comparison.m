function plot_aero_comparison(results, opts)
% plot_aero_comparison  气动结果可视化
%
% 图1: 三叶片总载荷时间序列（与BEM稳态基准对比）  —— 保持不变
% 图2: 单叶片动态载荷（2×3布局）
%       上行: 各叶片升力vs方位角 | 各叶片阻力vs方位角 | 各叶片扭矩vs方位角
%       下行: 各叶片升力vs时间   | 各叶片阻力vs时间   | 各叶片扭矩vs时间
%
% 数据来源：
%   results.blade.L/D/Torque : B×Nt，每行对应一片叶片，直接展向积分，不做叶片间平均
%   results.blade.psi         : B×Nt，对应方位角(rad)

if nargin < 2
    opts = struct();
end

t     = results.time(:);
Nt    = numel(t);
B     = results.meta.B;
Omega = results.meta.Omega;

% ===== 总载荷 =====
L_tot = results.total.Lift(:);
D_tot = results.total.Drag(:);
M_tot = results.total.Moment(:);

bemL = results.bemTotals.Lift;
bemD = results.bemTotals.Drag;
bemM = results.bemTotals.Moment;

smoothWin = getField(opts,'smoothWindow',1);
if smoothWin > 1
    L_tot = smoothdata(L_tot,'movmean',smoothWin);
    D_tot = smoothdata(D_tot,'movmean',smoothWin);
    M_tot = smoothdata(M_tot,'movmean',smoothWin);
end

titleStr = getField(opts,'title','气动性能参数对比分析');

% =========================================================
% 图1: 三叶片总载荷时间序列（保持不变）
% =========================================================
fig1 = figure('Color','w','Position',[80 80 1200 700],'Name','总载荷时间序列');

subplot(3,1,1);
plot(t, L_tot, 'b-', 'LineWidth', 1.8); hold on;
plot([t(1) t(end)], [bemL bemL], 'r--', 'LineWidth', 1.6);
xlabel('时间 (s)'); ylabel('升力 (N)');
legend('GDW总升力','BEM基准','Location','best');
title('三叶片总升力'); grid on; xlim([t(1) t(end)]);

subplot(3,1,2);
plot(t, D_tot, 'b-', 'LineWidth', 1.8); hold on;
plot([t(1) t(end)], [bemD bemD], 'r--', 'LineWidth', 1.6);
xlabel('时间 (s)'); ylabel('阻力 (N)');
legend('GDW总阻力','BEM基准','Location','best');
title('三叶片总阻力'); grid on; xlim([t(1) t(end)]);

subplot(3,1,3);
plot(t, M_tot, 'b-', 'LineWidth', 1.8); hold on;
plot([t(1) t(end)], [bemM bemM], 'r--', 'LineWidth', 1.6);
xlabel('时间 (s)'); ylabel('力矩 (N·m)');
legend('GDW总力矩','BEM基准','Location','best');
title('三叶片总力矩'); grid on; xlim([t(1) t(end)]);

sgtitle(fig1, titleStr);

% =========================================================
% 图2: 单叶片动态载荷（2×3，升力/阻力/扭矩 × 方位角/时间）
% =========================================================
if ~isfield(results,'blade')
    warning('results.blade 不存在，跳过单叶片载荷图。');
    return;
end

bL  = results.blade.L;       % B×Nt  升力(N)
bD  = results.blade.D;       % B×Nt  阻力(N)
bQ  = results.blade.Torque;  % B×Nt  扭矩(N·m)
psi = results.blade.psi;     % B×Nt  方位角(rad)

% 方位角归一化到 [0°, 360°)
psi_deg = mod(rad2deg(psi), 360);   % B×Nt

% 取末尾一个完整旋转周期用于方位角图（体现完整周期内的载荷变化）
T_rot    = 2*pi / Omega;
dt_sim   = mean(diff(t));
nPerRev  = max(2, round(T_rot / dt_sim));
idx_psi  = max(1, Nt - nPerRev + 1) : Nt;

% 颜色与图例
colors     = lines(B);
bladeNames = arrayfun(@(q) sprintf('叶片%d',q), 1:B, 'UniformOutput', false);

% 统一线宽与字体
LW_blade = 1.6;
FS_label = 10;
FS_title = 10;

fig2 = figure('Color','w','Position',[120 100 1380 820],'Name','单叶片动态载荷');
tl2  = tiledlayout(fig2, 2, 3, 'TileSpacing','compact','Padding','compact');

% ---------- 辅助函数：绘制"vs方位角"子图 ----------
function plot_vs_psi(tl, pos, data_B_Nt, psi_B_Nt, idx, cols, names, ylbl, ttl)
    ax = nexttile(tl, pos);
    hold(ax,'on');
    for qq = 1:size(data_B_Nt,1)
        pq = psi_B_Nt(qq, idx);
        dq = data_B_Nt(qq, idx);
        [ps, si] = sort(pq);
        plot(ax, ps, dq(si), '-', 'Color', cols(qq,:), 'LineWidth', LW_blade);
    end
    xlabel(ax, '方位角 (°)', 'FontSize', FS_label);
    ylabel(ax, ylbl, 'FontSize', FS_label);
    title(ax, ttl, 'FontSize', FS_title);
    legend(ax, names, 'Location','best', 'FontSize', 8);
    grid(ax,'on'); xlim(ax,[0 360]); xticks(ax, 0:60:360);
end

% ---------- 辅助函数：绘制"vs时间"子图 ----------
function plot_vs_time(tl, pos, data_B_Nt, t_vec, cols, names, ylbl, ttl)
    ax = nexttile(tl, pos);
    hold(ax,'on');
    for qq = 1:size(data_B_Nt,1)
        plot(ax, t_vec, data_B_Nt(qq,:), '-', 'Color', cols(qq,:), 'LineWidth', LW_blade);
    end
    xlabel(ax, '时间 (s)', 'FontSize', FS_label);
    ylabel(ax, ylbl, 'FontSize', FS_label);
    title(ax, ttl, 'FontSize', FS_title);
    legend(ax, names, 'Location','best', 'FontSize', 8);
    grid(ax,'on'); xlim(ax,[t_vec(1) t_vec(end)]);
end

% ===== 上行：各叶片 vs 方位角 =====
plot_vs_psi(tl2, 1, bL, psi_deg, idx_psi, colors, bladeNames, ...
    '升力 (N)', '各叶片升力 vs 方位角');

plot_vs_psi(tl2, 2, bD, psi_deg, idx_psi, colors, bladeNames, ...
    '阻力 (N)', '各叶片阻力 vs 方位角');

plot_vs_psi(tl2, 3, bQ, psi_deg, idx_psi, colors, bladeNames, ...
    '扭矩 (N·m)', '各叶片扭矩 vs 方位角');

% ===== 下行：各叶片 vs 时间 =====
plot_vs_time(tl2, 4, bL, t, colors, bladeNames, ...
    '升力 (N)', '各叶片升力 vs 时间');

plot_vs_time(tl2, 5, bD, t, colors, bladeNames, ...
    '阻力 (N)', '各叶片阻力 vs 时间');

plot_vs_time(tl2, 6, bQ, t, colors, bladeNames, ...
    '扭矩 (N·m)', '各叶片扭矩 vs 时间');

sgtitle(tl2, '单叶片动态载荷特性（方位角数据取末尾一个旋转周期）');

% =========================================================
% 导出
% =========================================================
if getField(opts,'export',false)
    outDir   = getField(opts,'outDir', pwd);
    baseName = getField(opts,'baseName','aero_comparison');
    formats  = getField(opts,'formats',{'png','pdf'});
    dpi      = getField(opts,'dpi',300);
    figs     = {fig1, fig2};
    suffixes = {'_total','_blade'};
    for fi = 1:2
        for k = 1:numel(formats)
            fmt     = lower(formats{k});
            outFile = fullfile(outDir, [baseName suffixes{fi} '.' fmt]);
            if strcmp(fmt,'png')
                exportgraphics(figs{fi}, outFile, 'Resolution', dpi);
            elseif strcmp(fmt,'pdf')
                exportgraphics(figs{fi}, outFile, 'ContentType','vector');
            else
                exportgraphics(figs{fi}, outFile);
            end
        end
    end
end
end
