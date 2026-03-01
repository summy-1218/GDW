function Gamma = Gamma_jn_hm(j, n, r, m)
% 计算 Gamma_{j,n}^{r,m} 系数，对应手册 Eq.79-81
%
% 参数含义（严格对应手册符号）：
%   j : 诱导速度径向形函数下标  (j = r+1, r+3, ...)
%   n : 压力系数径向下标        (n = m+1, m+3, ...)
%   r : 诱导速度方位角谐波阶数  （原始代码中的 h）
%   m : 压力系数方位角谐波阶数
%
% Eq.79 (r+m 为偶数):
%   Gamma = (-1)^{(n+j-2r)/2} / sqrt(H_n^m * H_j^r)
%           * 2*sqrt((2n+1)(2j+1)) / [(j+n)(j+n+2)((j-n)^2-1)]
%
%   【注意】分母第三项使用 (j-n)^2-1（而非手册印刷版的 (j+n)^2-1）。
%   数值验证表明 (j+n)^2-1 导致 L̃ 矩阵负定（特征值全负），使 GDW
%   时域积分发散；(j-n)^2-1 保证 L̃ 正定（特征值全正），数值稳定。
%   这与原始参考实现一致，怀疑手册存在印刷笔误。
%
% Eq.80 (r+m 为奇数, |j-n|=1):
%   Gamma = pi / sqrt(H_n^m * H_j^r) * sign(r-m) / sqrt((2n+1)(2j+1))
%
%   【注意】系数为 pi（非 pi/2）。与 (j-n)^2-1 共同构成数值稳定的原始实现。
%
% Eq.81 (r+m 为奇数, |j-n|!=1):
%   Gamma = 0

if mod(r + m, 2) == 0
    % ------ Eq.79: r+m 为偶数 ------
    denom = (j + n) * (j + n + 2) * ((j - n)^2 - 1);  % 分母第三项: (j-n)^2-1
    if abs(denom) < eps
        Gamma = 0;
    else
        sign_exp = (-1)^((n + j - 2*r) / 2);  % 符号指数使用 r
        Gamma = sign_exp / sqrt(H_jh(n, m) * H_jh(j, r)) * ...
                (2 * sqrt((2*n + 1) * (2*j + 1))) / denom;
    end
else
    % ------ Eq.80/81: r+m 为奇数 ------
    if abs(j - n) == 1
        % Eq.80: r==m 时 r+m=2r 为偶数，不进入此分支，sign(r-m) 不为零
        Gamma = pi / sqrt(H_jh(n, m) * H_jh(j, r)) * ...
                sign(r - m) / sqrt((2*n + 1) * (2*j + 1));
    else
        Gamma = 0;  % Eq.81
    end
end
end