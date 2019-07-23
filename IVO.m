pos = zeros(8, 200);

Rd = [5, 0];
Rg = [5, 20];
Ro = [5, 10];

vc = 10;
v = [vc, vc];

vobs = [0, 0];

Rsense = 3;
R = 2;

t = 0.01;

[d, dx, dy] = distance(Rg, Rd);
[r, rx, ry] = distance(Ro, Rd);

i = 1;

while d>0.5
    pos(1, i) = Rd(1);
    pos(2, i) = Rd(2);
       
    i = i + 1;
     
    if Rsense >= r
        [ux, uy] = optimize(v, Rg, Rd, t, rx, ry, r, R, vobs);   
        v = [v(1) + ux, v(2) + uy];
        Rd(1) = Rd(1) + v(1)*t;
        Rd(2) = Rd(2) + v(2)*t;
        pos(7, i) = 1;
        pos(3, i) = v(1);
        pos(4, i) = v(2);
    else
        [v1, v2] = optimize1(v, Rg, Rd, t);   
        v = [v1, v2];
%         v(1) = vc*(rx/r);
%         v(2) = vc*(ry/r);
        Rd(1) = Rd(1) + v(1)*t;
        Rd(2) = Rd(2) + v(2)*t;
        pos(8, i) = 1;
        pos(3, i) = v(1);
        pos(4, i) = v(2);
    end
    
    
    [d, dx, dy] = distance(Rg, Rd);
    [r, rx, ry] = distance(Ro, Rd);
    pos(5, i) = r;
    pos(6, i) = d;

end    
    
disp(pos);
figure;
plot(pos(1, :), pos(2, :), 'o');
function [d, dx, dy] = distance(X, Y)
    dx = (X(1) - Y(1));
    dy = (X(2) - Y(2));
    d = sqrt((dx)^2 + (dy)^2);
end

function [ux, uy] = optimize(v, Rg, Rd, t, rx, ry, r, R, vobs)
    objective = @(x) ((Rg(1) - (Rd(1) + x(1)*t))^2 + (Rg(2) - (Rd(2) + x(2)*t))^2);
    
    x0 = [v(1), v(2)];
    
    disp(['Initial Objective: ' num2str(objective(x0))])
    
    lb = -10*[1, 1];
    ub = 10*[1, 1];  
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,@(x)nlcon_ivo(x,R,rx,ry,r, vobs));
      
    ux = x(1);
    uy = x(2);
    
end

function [vx, vy] = optimize1(v, Rg, Rd, t)
    objective = @(x) ((Rg(1) - (Rd(1) + x(1)*t))^2 + (Rg(2) - (Rd(2) + x(2)*t))^2);
    
    x0 = [v(1), v(2)];
    
    disp(['Initial Objective: ' num2str(objective(x0))])
    
    lb = -10*[1, 1];
    ub = 10*[1, 1];  
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub);
      
    vx = x(1);
    vy = x(2);
    
end
