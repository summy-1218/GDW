function S_h = default_S_h(N_h, rPowerMax)
S_h = zeros(N_h+1,1);
for h = 0:N_h
    S_h(h+1) = floor((rPowerMax - h)/2) + 1;
    S_h(h+1) = max(S_h(h+1), 1);  % 手册Table1: 每个谐波至少保留1个形函数
end
end