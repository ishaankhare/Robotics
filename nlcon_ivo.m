function [c, ceq] = nlcon_ivo(x, R, rx, ry, r, vobs)
    vrel = (vobs(1) - x(1))^2 + (vobs(2) - x(2))^2; 
    c = vrel*R^2 - vrel*r^2 + (rx*(vobs(1) - x(1)) + ry*(vobs(2) - x(2)))^2;
    ceq = [];
end

