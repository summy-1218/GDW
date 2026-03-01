function bemTotals = compute_bem_totals(r, c, theta, thickness, dr, B, RootOffset, V0, Omega, rho, a, ap)
N = numel(r);
L = 0; D = 0; M = 0; T = 0; Q = 0;
for i = 1:N
    if r(i) < RootOffset
        continue;
    end
    Vax = V0 * (1 - a(i));
    Vtan = Omega * r(i) * (1 + ap(i));
    Phi = atan2(Vax, Vtan);
    Alpha = Phi - deg2rad(theta(i));
    [Cl, Cd, Cm] = aeroInterp(thickness(i), Alpha);
    Vrel = hypot(Vax, Vtan);
    dL = 0.5 * rho * Vrel^2 * c(i) * Cl * dr(i);
    dD = 0.5 * rho * Vrel^2 * c(i) * Cd * dr(i);
    dM = 0.5 * rho * Vrel^2 * c(i)^2 * Cm * dr(i);
    dT = dL * cos(Phi) + dD * sin(Phi);
    dFt = dL * sin(Phi) - dD * cos(Phi);
    L = L + dL * B;
    D = D + dD * B;
    M = M + dM * B;
    T = T + dT * B;
    Q = Q + dFt * r(i) * B;
end
bemTotals.Lift = L;
bemTotals.Drag = D;
bemTotals.Moment = M;
bemTotals.Thrust = T;
bemTotals.Torque = Q;
end