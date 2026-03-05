function [Mdiag_c, Mdiag_s, LtildeC, LtildeS, offC, offS] = precompute_gdw_matrices(mode, chi)
% 预计算 GDW 全局矩阵 M、LtildeC 和 LtildeS
%
% 符号约定（严格对应手册 Eq.60-81）：
%   r   : 诱导速度方位角谐波阶（行索引，h=0..N_h）
%   m   : 压力系数方位角谐波阶（列索引，h=0..N_h），与 r 独立
%   j   : 诱导速度径向形函数下标（属于谐波 r 的第 a 个元素）
%   n   : 压力系数径向下标（属于谐波 m 的第 b 个元素）
%
% 输出：
%   Mdiag_c  : 余弦状态方程的 M 矩阵对角向量（对应 r=0,1,...,N_h 拼接）
%   Mdiag_s  : 正弦状态方程的 M 矩阵对角向量（对应 r=1,...,N_h 拼接）
%   LtildeC  : 全局余弦 Ltilde 矩阵（N_c × N_c），包含所有 r×m 跨谐波耦合
%   LtildeS  : 全局正弦 Ltilde 矩阵（N_s × N_s），包含所有 r×m 跨谐波耦合
%   offC(h+1): 谐波 h 的余弦状态在全局向量中的起始偏移（0-based）
%   offS(h)  : 谐波 h（h=1..N_h）的正弦状态在全局正弦向量中的起始偏移（0-based）
%
% 物理意义：
%   手册 Eq.74-76 中，Ltilde 的行 (r,j) 与列 (m,n) 完全独立，
%   同一行 r 的元素对所有列 m 均可能非零（跨谐波耦合）。
%   旧代码固定 r=m=h，只保留对角块，丢失了全部跨谐波耦合。
%
% Eq.74  (r=0 行): sigma_{jn}^{0m,c} = X^m * Gamma_{jn}^{0m}
% Eq.75  (r>0 行，余弦): T_{jn}^{rm,c} = (X^|m-r| + (-1)^l * X^|m+r|) * Gamma_{jn}^{rm}
% Eq.76  (r>0 行，正弦): pi_{jn}^{rm,s} = (X^|m-r| - (-1)^l * X^|m+r|) * Gamma_{jn}^{rm}
% Eq.77  : l = min(r, m)
% Eq.78  : X = tan(|chi|/2)，X^0 = 1（约定）

N_h = mode.N_h;
X   = tan(abs(chi) / 2);            % Eq.78: X = tan(|chi|/2)

% ===== 计算各谐波的径向形函数列表和全局偏移 =====
% 余弦状态: h = 0, 1, ..., N_h
offC = zeros(1, N_h+2);             % offC(h+1) = 谐波 h 在余弦全局向量的起始位置(0-based)
for h = 0:N_h
    offC(h+2) = offC(h+1) + numel(mode.jList{h+1});
end
N_c = offC(N_h+2);                  % 余弦状态总数

% 正弦状态: h = 1, 2, ..., N_h
offS = zeros(1, N_h+1);             % offS(h) = 谐波 h（h≥1）在正弦全局向量的起始位置(0-based)
for h = 1:N_h
    offS(h+1) = offS(h) + numel(mode.jList{h+1});
end
N_s = offS(N_h+1);                  % 正弦状态总数

% ===== M 矩阵（对角）Eq.71: zeta_j^h = (2/pi)*H_j^h =====
Mdiag_c = zeros(N_c, 1);
for h = 0:N_h
    jList = mode.jList{h+1};
    for a = 1:numel(jList)
        j_a = jList(a);
        Mdiag_c(offC(h+1) + a) = (2/pi) * H_jh(j_a, h);
    end
end

Mdiag_s = zeros(N_s, 1);
for h = 1:N_h
    jList = mode.jList{h+1};
    for a = 1:numel(jList)
        j_a = jList(a);
        Mdiag_s(offS(h) + a) = (2/pi) * H_jh(j_a, h);
    end
end

% ===== LtildeC 全局矩阵构建 =====
% 行索引: (r, j) — 诱导速度谐波 r 的径向分量 j
% 列索引: (m, n) — 压力谐波 m 的径向分量 n
% 每一行 r 对所有列 m（0..N_h）均可能有耦合贡献
Lc = zeros(N_c, N_c);

for r = 0:N_h
    jList_r = mode.jList{r+1};
    Sr = numel(jList_r);
    if Sr == 0; continue; end

    for m = 0:N_h
        jList_m = mode.jList{m+1};
        Sm = numel(jList_m);
        if Sm == 0; continue; end

        l = min(r, m);                      % Eq.77
        pow_mr = Xpow(X, abs(m - r));       % X^|m-r|，X^0=1
        pow_pr = Xpow(X, abs(m + r));       % X^|m+r|，X^0=1
        sign_l = (-1)^l;

        for a = 1:Sr
            j = jList_r(a);
            row = offC(r+1) + a;
            for b = 1:Sm
                n = jList_m(b);
                col = offC(m+1) + b;

                G = Gamma_jn_hm(j, n, r, m);

                if r == 0
                    % Eq.74: sigma_{jn}^{0m,c} = X^m * Gamma
                    coeff = Xpow(X, m);
                else
                    % Eq.75: T_{jn}^{rm,c} = (X^|m-r| + (-1)^l * X^|m+r|) * Gamma
                    coeff = (pow_mr + sign_l * pow_pr);
                end
                Lc(row, col) = coeff * G;
            end
        end
    end
end

% ===== LtildeS 全局矩阵构建 =====
% 行索引: (r, j) — r=1..N_h
% 列索引: (m, n) — m=1..N_h（正弦方程中 m=0 无正弦项）
Ls = zeros(N_s, N_s);

for r = 1:N_h
    jList_r = mode.jList{r+1};
    Sr = numel(jList_r);
    if Sr == 0; continue; end

    for m = 1:N_h
        jList_m = mode.jList{m+1};
        Sm = numel(jList_m);
        if Sm == 0; continue; end

        l = min(r, m);                      % Eq.77
        pow_mr = Xpow(X, abs(m - r));       % X^|m-r|
        pow_pr = Xpow(X, abs(m + r));       % X^|m+r|
        sign_l = (-1)^l;

        for a = 1:Sr
            j = jList_r(a);
            row = offS(r) + a;
            for b = 1:Sm
                n = jList_m(b);
                col = offS(m) + b;

                G = Gamma_jn_hm(j, n, r, m);
                % Eq.76: pi_{jn}^{rm,s} = (X^|m-r| - (-1)^l * X^|m+r|) * Gamma
                coeff = (pow_mr - sign_l * pow_pr);
                Ls(row, col) = coeff * G;
            end
        end
    end
end

% 加正则化微扰防止精确奇异（eps 量级扰动不影响物理结果）
LtildeC = Lc + 1e-8 * eye(N_c);
if N_s > 0
    LtildeS = Ls + 1e-8 * eye(N_s);
else
    LtildeS = zeros(0);
end

end

% ===== 辅助函数: X^p，处理 X^0 = 1 的约定 =====
function val = Xpow(X, p)
% 手册注释：X^m 和 X^|m-r| 在指数为 0 时取 1（即使 X=0）
if p == 0
    val = 1;
else
    val = X^p;
end
end
