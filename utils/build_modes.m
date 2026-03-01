function mode = build_modes(N_h, S_h)
% mode包含方位角谐波的最高数目Nh，第h次谐波对应的形函数数目Sh
mode = struct();
mode.N_h = N_h;
mode.S_h = S_h;
mode.jList = cell(N_h+1,1);
for h = 0:N_h
    Sh = S_h(h+1);
    j = (h+1):2:(2*Sh+h-1);
    mode.jList{h+1} = j(:);
end
end