clear; clc; close all;
addpath('./utils');
%% 输入参数
% 叶片几何参数
B = 3;              % 叶片数量
R_hub = 1.56;
load beamModel.mat;
r = beamModel.distance; % 从叶根到叶尖的径向站坐标 (m)。
dr = [diff(r); diff(r(end-1:end))]; % 每个叶素的宽度 (m)
R = r(end);           % 叶片半径 (m)
% 各径向站的叶片几何特性 (这些数组的长度应与r相同)
c = beamModel.chord; % 弦长分布示例 (m)，从叶根到叶尖线性减小
theta = beamModel.twist;  % 几何安装角分布示例 (度)，从叶根到叶尖扭角变化
RootOffset = 0.05*R; % 叶根截止位置 (m)，通常为0.2R
% 运行条件
V0 = 10.0;          % 来流风速 (m/s)
Omega = 1.2566; % 风轮角速度 (rad/s)，假设尖速比TSR=2
rho = 1.225;        % 空气密度 (kg/m^3)

% 翼型气动数据表 (关键！需要你提供实际翼型的数据)
% 示例: 假设使用NACA0012翼型，但数据是简化的。请替换为你的翼型数据。
load cl_cd_cm.mat;
% Alpha_Data = cl_cd_cm(:,1)'; % 攻角示例 (度)
% Cl_Data = cl_cd_cm(:,2)'; % 升力系数示例
% Cd_Data = cl_cd_cm(:,3)'; % 阻力系数示例

%% 初始化诱导因子等数组
Nstations = length(r);
a = zeros(Nstations, 1);       % 轴向诱导因子
ap = zeros(Nstations, 1);      % 周向诱导因子

% 迭代参数
maxIter = 100;      % 最大迭代次数
tolerance = 1e-5;   % 收敛容差
relax = 0.4;        % 松弛因子，防止迭代振荡，通常在0.1到0.5之间

% 为存储结果预分配数组
Phi = zeros(Nstations, 1);
Alpha = zeros(Nstations, 1);
Cl = zeros(Nstations, 1);
Cd = zeros(Nstations, 1);
dT = zeros(Nstations, 1); % 单个叶素产生的推力
dQ = zeros(Nstations, 1); % 单个叶素产生的扭矩(对风轮中心)
a = zeros(Nstations, 1);
ap = zeros(Nstations, 1);

fprintf('开始BEMT迭代计算...\n');
for i = 1:Nstations
    ri = r(i);
    if ri < RootOffset
        continue; % 跳过叶根截止点以内的部分
    end
    
    % 计算局部运动速度,先留好接口
    Ve_op = 0;
    Ve_ip = 0;

    % 计算局部实度
    sigma = (B * c(i)) / (2 * pi * ri);
    
    % 计算当地尖速比
    lambda = Omega*ri/V0;

    converged = false;

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
        [Cl_i,Cd_i,~] = aeroInterp(beamModel.Thickness(i), Alpha_i);

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

        % 
        % if Kp <= 0
        %     ap_new = 0;
        % else
        %     
        % end        

        % 使用松弛因子更新诱导因子
        a(i) = (1 - relax) * a_old + relax * a_new;
        ap(i) = (1 - relax) * ap_old + relax * ap_new;


        % 检查收敛性
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
    

    % 叶素理论: 计算微元推力dT和扭矩dQ
    dL = 0.5 * rho * V_local^2 * c(i) * Cl_i * dr(i); % 微元升力
    dD = 0.5 * rho * V_local^2 * c(i) * Cd_i * dr(i); % 微元阻力

    dT_i = dL * cos(Phi_i) + dD * sin(Phi_i); % 推力方向微元力
    dQ_i = (dL * sin(Phi_i) - dD * cos(Phi_i)) * ri; % 扭矩微元力 (力 * 半径)

    % 存储最终结果
    Phi(i) = rad2deg(Phi_i);
    Alpha(i) = Alpha_i;
    Cl(i) = Cl_i;
    Cd(i) = Cd_i;
    dT(i) = dT_i * B; % 单个叶素微元推力 * 叶片数
    dQ(i) = dQ_i * B; % 单个叶素微元扭矩 * 叶片数
end
fprintf('BEMT迭代计算完成。\n');

%% 积分得到总推力和扭矩
Total_Thrust = sum(dT); % 总推力 (N)
Total_Torque = sum(dQ); % 总扭矩 (Nm)，对风轮中心
Power = Total_Torque * Omega; % 气动功率 (W)

fprintf('\n========== 最终结果 ==========\n');
fprintf('总推力: %.2f N\n', Total_Thrust);
fprintf('总扭矩: %.2f Nm\n', Total_Torque);
fprintf('气动功率: %.2f MW\n', Power/1000000);
fprintf('功率系数 Cp: %.4f\n', Power / (0.5*rho*pi*R^2*V0^3));

%% 绘制结果
figure('Position', [100, 100, 1200, 800])

subplot(2, 3, 1)
plot(r/R, a, 'b-', 'LineWidth', 2); hold on;
plot(r/R, ap, 'r-', 'LineWidth', 2);
xlabel('无量纲半径 r/R');
ylabel('诱导因子');
legend('轴向诱导因子 a', '周向诱导因子 a''', 'Location', 'best');
grid on;
title('诱导因子沿展向分布');

subplot(2, 3, 2)
plot(r/R, rad2deg(Alpha), 'g-', 'LineWidth', 2);
xlabel('无量纲半径 r/R');
ylabel('攻角 (度)');
grid on;
title('攻角沿展向分布');

subplot(2, 3, 3)
plot(r/R, Cl, 'c-', 'LineWidth', 2); hold on;
plot(r/R, Cd, 'm-', 'LineWidth', 2);
xlabel('无量纲半径 r/R');
ylabel('气动系数');
legend('升力系数 C_l', '阻力系数 C_d', 'Location', 'best');
grid on;
title('气动系数沿展向分布');

subplot(2, 3, 4)
plot(r/R, dT./dr/B, 'k-', 'LineWidth', 2); % 显示单个叶素的贡献
xlabel('无量纲半径 r/R');
ylabel('单位展长推力 (N/m)');
grid on;
title('推力分布 (单叶片)');

subplot(2, 3, 5)
plot(r/R, dQ./dr/B, 'k-', 'LineWidth', 2); % 显示单个叶素的贡献
xlabel('无量纲半径 r/R');
ylabel('单位展长扭矩 (Nm/m)');
grid on;
title('扭矩分布 (单叶片)');

subplot(2, 3, 6)
CT_local = dT ./ (0.5 * rho * V0^2 * 2*pi*r.*dr); % 局部推力系数
plot(r/R, CT_local, 'Color', [0.7, 0.5, 0.1], 'LineWidth', 2);
xlabel('无量纲半径 r/R');
ylabel('局部推力系数 C_T_,_l_o_c_a_l');
grid on;
title('局部推力系数分布');

sgtitle('风机叶片气动性能BEMT计算结果');

load BEM.mat;
load VLM.mat
figure
plot(r/R, Cl, 'b*', 'LineWidth', 2); hold on;
plot(r/R, Cd, 'm*', 'LineWidth', 2);
plot(VLM(:,1)/VLM(end,1), VLM(:,end-3), 'b+', 'LineWidth', 2);
plot(VLM(:,1)/VLM(end,1), VLM(:,end-2), 'm+', 'LineWidth', 2);
plot(BEM(:,1)/BEM(end,1), BEM(:,end-2), 'b-', 'LineWidth', 2);
plot(BEM(:,1)/BEM(end,1), BEM(:,end-1), 'm-', 'LineWidth', 2);
xlabel('无量纲半径 r/R');
ylabel('气动系数');
legend('BEM:升力系数 C_l', 'BEM:阻力系数 C_d', 'VLM:升力系数 C_l', 'VLM:阻力系数 C_d',...
    'BEM_bladed:升力系数 C_l', 'BEM_bladed:阻力系数 C_d','Location', 'best');
grid on;
title('气动系数沿展向分布');