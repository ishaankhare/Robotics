function [c, ceq] = nlcon_Vo_2(x,R, dx, dy, v,rx,ry,r, r1, rx1, ry1, vobs, vobs1, state,  RX, RY, RX1, RY1,vrx, vry, vrx1, vry1)
    vrel = (x(1) - vobs(1))^2 + (x(2) - vobs(2))^2;
    vrel1 = (x(1) - vobs1(1))^2 + (x(2) - vobs1(2))^2;
    
    if state == 0       
        c(1) = vrel*R^2 - vrel*r^2 + (rx*(x(1) - vobs(1)) + ry*(x(2) - vobs(2)))^2;
        c(2) = RX*vrx + RY*vry;
        ceq = [];       
    end    
    if state == 1 
        c(1) = vrel1*R^2 - vrel1*r1^2 + (rx1*(x(1) - vobs1(1)) + ry1*(x(2) - vobs1(2)))^2;
        c(2) = RX1*vrx1 + RY1*vry1;
        ceq = [];    
    end
    if state == 2  
        c(1) = vrel*R^2 - vrel*r^2 + (rx*(x(1) - vobs(1)) + ry*(x(2) - vobs(2)))^2;
        c(2) = vrel1*R^2 - vrel1*r1^2 + (rx1*(x(1) - vobs1(1)) + ry1*(x(2) - vobs1(2)))^2;
        c(3) = RX*vrx + RY*vry;
        c(4) = RX1*vrx1 + RY1*vry1;
        ceq = [];    
    end
    if state == 3    
        c = [];
        ceq = [];
    end 
    
    
    
end 






