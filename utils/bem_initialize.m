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

    converged = false;
    for iter = 1:maxIter
        a_old = a(i);
        ap_old = ap(i);

        % 1. 计算入流角 Phi (弧度)
        V_local_axial = V0 * (1 - a(i))+Ve_op;
        V_local_tangential = Omega * ri * (1 + ap(i)) + Ve_ip;
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
        %F_root = (2/pi) * real(acos(exp(-B/2*(ri-R_hub)/(R_hub*sin(Phi_i)))));
        % 手册 Eq.13: f = B/2 * (r - R_hub) / (r * sin(phi))
        F_root = (2/pi) * real(acos(exp(-B/2 * (ri - R_hub) / (ri * abs(sin(Phi_i)) + eps))));
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
            % 手册 Eq.27: 使用组合损失因子 F = F_tip * F_root
            K = (4 * F * (sin(Phi_i))^2) / (sigma * (Cl_i * cos(Phi_i) + Cd_i * sin(Phi_i)));
            a_new = 1 / (K + 1);
        end
        % 手册 Eq.28: 切向诱导因子同样使用 F 
        Kp = (4 * F * sin(Phi_i) * cos(Phi_i)) / (sigma * (Cl_i * sin(Phi_i) - Cd_i * cos(Phi_i)) + eps);
        ap_new = 1 / (Kp - 1 + eps); % 加eps避免除零;

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