function [Mdiag, LtildeC, LtildeS] = precompute_gdw_matrices(mode, chi)
% 预计算 GDW 矩阵 M（表观质量）、Ltilde^c 和 Ltilde^s（尾流偏斜增益矩阵）
%
% 符号约定（严格对应手册 Eq.60-81）：
%   h   : 方位角谐波阶数（对应手册中的 r，此处用 h 区分循环变量与函数参数）
%   j   : 诱导速度径向形函数下标（外层循环）
%   n   : 压力系数径向下标（内层循环）
%   r=h : Gamma_jn_hm 第3参数（诱导速度谐波阶，等于当前谐波 h）
%   m=h : Gamma_jn_hm 第4参数（压力谐波阶，等于当前谐波 h）
%
% M矩阵（Eq.70-71）：对角元素 zeta_j^h = (2/pi)*H_j^h
% Ltilde^c（Eq.74-75）：h=0 用 Eq.74，h>0 用 Eq.75（余弦耦合）
% Ltilde^s（Eq.76）   ：h>0 用 Eq.76（正弦耦合），h=0 时无正弦项

N_h = mode.N_h;
Mdiag  = cell(N_h+1, 1);
LtildeC = cell(N_h+1, 1);
LtildeS = cell(N_h+1, 1);
X = tan(abs(chi) / 2);          % Eq.78: X = tan(|chi|/2)

for h = 0:N_h
    jList = mode.jList{h+1};    % 该谐波对应的径向下标列表
    Sh = numel(jList);

    %% --- M 矩阵（对角）Eq.71: zeta_j^h = (2/pi)*H_j^h ---
    Mdiag{h+1} = zeros(Sh, 1);
    for a = 1:Sh
        j_a = jList(a);                         % 诱导速度径向下标
        Mdiag{h+1}(a) = (2/pi) * H_jh(j_a, h); % H_j^h，j在前，h（谐波阶）在后
    end

    %% --- Ltilde 矩阵 ---
    Lc = zeros(Sh, Sh);
    Ls = zeros(Sh, Sh);
    for a = 1:Sh
        j = jList(a);           % 诱导速度径向下标（行）
        for b = 1:Sh
            n = jList(b);       % 压力系数径向下标（列）

            % Gamma_{jn}^{r,m}，此处 r=h（诱导速度谐波），m=h（压力谐波）
            % 同谐波自耦合：r+m=2h 必为偶数，始终走 Eq.79 分支
            G = Gamma_jn_hm(j, n, h, h);

            if h == 0
                % Eq.74: [~sigma_{jn}^{0,0}]^c = X^0 * Gamma = Gamma
                Lc(a, b) = G;               % X^0 = 1
            else
                % Eq.75: [~T_{jn}^{h,h}]^c = (X^|h-h| + (-1)^h * X^(h+h)) * Gamma
                %       = (1 + (-1)^h * X^(2h)) * Gamma
                % Eq.76: [~pi_{jn}^{h,h}]^s = (X^|h-h| - (-1)^h * X^(h+h)) * Gamma
                %       = (1 - (-1)^h * X^(2h)) * Gamma
                % l = min(r,m) = min(h,h) = h  （Eq.77）
                pow_X = X^(2*h);            % X^|m+r| = X^(2h)
                sign_l = (-1)^h;            % (-1)^l，l=h
                Lc(a, b) = (1 + sign_l * pow_X) * G;   % Eq.75
                Ls(a, b) = (1 - sign_l * pow_X) * G;   % Eq.76
            end
        end
    end

    % 加正则化微扰防止精确奇异（eps 量级扰动不影响物理结果）
    LtildeC{h+1} = Lc + 1e-8 * eye(Sh);
    if h > 0
        LtildeS{h+1} = Ls + 1e-8 * eye(Sh);
    else
        LtildeS{h+1} = [];                  % h=0 无正弦项
    end
end
end