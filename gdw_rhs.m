function dAlpha = gdw_rhs(alpha, tauC, Mdiag, Ltilde, Vvec)
% GDW状态方程右端: M * dalpha/dt = 0.5*tau - Ltilde^{-1} * V * alpha
% 手册 Eq.95/96: [M]{alpha}' + [Ltilde]^{-1}[V]{alpha} = 0.5*{tau}
rhs = 0.5 * tauC;

% 防奇异: 检查Ltilde条件数，过大时加正则化
if rcond(Ltilde) < 1e-12
    reg = 1e-6 * max(abs(diag(Ltilde))) * eye(size(Ltilde));
    Lterm = (Ltilde + reg) \ (Vvec .* alpha);
else
    Lterm = Ltilde \ (Vvec .* alpha);
end

dAlpha = (rhs - Lterm) ./ (Mdiag + eps);

% 发散检测: 若导数过大则截断
maxRate = 1e6;
dAlpha(~isfinite(dAlpha)) = 0;
dAlpha = max(min(dAlpha, maxRate), -maxRate);
end