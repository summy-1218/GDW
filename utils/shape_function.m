function phi = shape_function(h, j, rhat)
% 形函数: 式(4-59)
% phi_j^h(r_hat) = sqrt((2j+1)H_j^h) * sum_{q=h,h+2,...}^{j-1} r^q * (-1)^((q-h)/2) * (j+q)!! / ((q-h)!! (q+h)!! (j-q-1)!!)
H = H_jh(j, h);
pref = sqrt((2*j+1) * H);
phi = zeros(size(rhat));
for q = h:2:(j-1)
    num = double_factorial(j+q);
    den = double_factorial(q-h) * double_factorial(q+h) * double_factorial(j-q-1);
    coef = (-1)^((q-h)/2) * num / den;
    phi = phi + coef * (rhat.^q);
end
phi = pref * phi;
end