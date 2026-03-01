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