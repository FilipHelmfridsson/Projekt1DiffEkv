% function ds = fvel(t, s, b, aL, aR, wL, wR)
%     B = (wR + wL) ./ 2;
%     D = (wR - wL) ./ b;
%     C = (aR - aL) ./ (2*b);
%     A = (aR + aL) ./ 2;

%     dxdt = (A * t + B) * cos(C * (t^2) + D * t + theta0);
%     dydt = (A * t + B) * sin(C * (t^2) + D * t + theta0);
%     dthetadt = 2 * C * t + D;

%     ds = [dxdt,dydt,dthetadt];
% end

function ds = fvel(t, s, b, aL, aR, wL, wR)
    theta = s(3);
    vR = aR * t + wR;
    vL = aL * t + wL;
    v = (vR + vL) / 2;
    
    dxdt = v * cos(theta);
    dydt = v * sin(theta);
    dthetadt = (vR - vL) / b;
    
    ds = [dxdt, dydt, dthetadt];
end
