function [c, ceq] = nlcon(x, R, rx, ry, r)
    vrel = x(1)^2 + x(2)^2;
    c = vrel*R^2 - vrel*r^2 + (rx*x(1) + ry*x(2))^2;
    ceq = [];
end        