function [yNew, fHist] = abm4_update(y, fNow, fHist, dtHat, it)
% Adams-Bashforth 4阶显式格式 (自洽，无需回调)
% fHist列: [f_{n}, f_{n-1}, f_{n-2}, f_{n-3}]  (最新在左)
%
% 启动阶段(it<4): 显式欧拉，待历史积累后切换到AB4
% it>=4: 4阶Adams-Bashforth
%   yNew = y + dt/24*(55*f_n - 59*f_{n-1} + 37*f_{n-2} - 9*f_{n-3})
%
% 参考: AeroDyn Theory Manual §2.2.8, Press et al. 1982

if it < 4
    % 启动阶段: 显式欧拉，逐步积累历史
    yNew = y + dtHat * fNow;
else
    f_n  = fHist(:,1);  % f_{n-1}
    f_n1 = fHist(:,2);  % f_{n-2}
    f_n2 = fHist(:,3);  % f_{n-3}

    % Adams-Bashforth 4 (纯显式，稳定且自洽)
    yNew = y + dtHat/24 * (55*fNow - 59*f_n + 37*f_n1 - 9*f_n2);
end

% 更新历史 (保留最近4步，最新在第1列)
fHist = [fNow, fHist(:,1:3)];
end
