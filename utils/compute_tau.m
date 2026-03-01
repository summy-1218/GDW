function [tauC, tauS] = compute_tau(dT_all, B, Omega, R, rho, psi_blade, mode, phiMat)
N_h = mode.N_h;
Nst = size(dT_all, 1);
if size(dT_all, 2) ~= B
    error('compute_tau: dT_all列数(%d)与叶片数B(%d)不一致', size(dT_all,2), B);
end
if numel(psi_blade) ~= B
    error('compute_tau: psi_blade长度(%d)与叶片数B(%d)不一致', numel(psi_blade), B);
end
if size(phiMat{1}, 1) ~= Nst
    error('compute_tau: phiMat径向站数(%d)与dT_all行数(%d)不一致', size(phiMat{1},1), Nst);
end
% 手册Eq.92: 归一化系数基准 = 1 / (pi * rho * Omega^2 * R^4)
% A = pi*R^2 => rho*A*(Omega*R)^2 = pi*rho*Omega^2*R^4 (正确)
A = pi * R^2;
scale = 1 / (rho * A * (Omega * R)^2 + eps);  % = 1/(pi*rho*Omega^2*R^4)

% 归一化叶素推力 (手册 Eq.91): Fhat = L_i^q / (rho * A * (Omega*R)^2)
% dT_all 单位: N (有量纲轴向推力), scale = 1/(pi*rho*Omega^2*R^4)
% Fhat 为无量纲, 对应手册"normalized element thrust force"
Fhat = dT_all * scale; % Nstations x B, 无量纲

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
        % scale = 1/(pi*rho*Omega^2*R^4), 与手册Eq.92/94分母对应
        % m=0 (Eq.94): tau = sumC / (2*pi*rho*Omega^2*R^4) = scale/2 * sumC
        % m!=0(Eq.92): tau = sumC / (pi*rho*Omega^2*R^4)   = scale   * sumC
        if h == 0
            tauC{h+1}(nIdx) = 0.5 * sumC;
        else
            tauC{h+1}(nIdx) = sumC;
            tauS{h+1}(nIdx) = sumS;
        end
    end
end
end