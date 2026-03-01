function phiMat = precompute_shape_functions(mode, rhat)
N_h = mode.N_h;
phiMat = cell(N_h+1,1);
for h = 0:N_h
    jList = mode.jList{h+1};
    Sh = numel(jList);
    Phi = zeros(numel(rhat), Sh);
    for idx = 1:Sh
        j = jList(idx);
        Phi(:,idx) = shape_function(h, j, rhat);
    end
    phiMat{h+1} = Phi;
end
end
