function H = H_jh(j,h)
% 详见理论文档式（4-60）
H = double_factorial(j+h-1) * double_factorial(j-h-1) / (double_factorial(j+h) * double_factorial(j-h) + eps);
end