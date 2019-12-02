pos = zeros(12, 200);
ob = zeros(1, 200);

Rd = [5, 0];
Rg = [5, 25];
Ro = [5, 10];
Ro1 = [20 ,15];

vc = 10;
v = [0, vc];

Rsense = 5;
R = 2;

t = 0.01;

[d, dx, dy] = distance(Rg, Rd);
[r, rx, ry] = distance(Ro, Rd);
[r1, rx1, ry1] = distance(Ro1, Rd);

i = 1;
j = 1;

state = 3;

pos(1, i) = Rd(1);
pos(2, i) = Rd(2);
pos(3, i) = v(1);
pos(4, i) = v(2);
pos(5, i) = d;
pos(6, i) = state;
pos(7, i) = Ro(1);
pos(8, i) = Ro(2);
pos(9, i) = Ro1(1);
pos(10, i) = Ro1(2);

vobs = [0, -5];
vobs1 = [-5, 0];


while d>0.5
    Ro(1) = Ro(1) + vobs(1)*t;
    Ro(2) = Ro(2) + vobs(2)*t;
        
    Ro1(1) = Ro1(1) + vobs1(1)*t;
    Ro1(2) = Ro1(2) + vobs1(2)*t;    
    
    vrx = pos(3, i) - vobs(1);
    vry = pos(4, i) - vobs(2);
    
    RX = pos(1, i) - pos(7, i);
    RY = pos(2, i) - pos(8, i);
    
    vrx1 = pos(3, i) - vobs1(1);
    vry1 = pos(4, i) - vobs1(2);
    
    RX1 = pos(1, i) - pos(9, i);
    RY1 = pos(2, i) - pos(10, i);
    
    i = i + 1;
         
    if Rsense > r && Rsense < r1  
        state = 0;
    elseif Rsense < r && Rsense > r1
        state = 1;
    elseif Rsense > r && Rsense > r1  
        state = 2;
    else 
        state = 3;
    end
        
    [d, dx, dy] = distance(Rg, Rd);
    [r, rx, ry] = distance(Ro, Rd);
    [r1, rx1, ry1] = distance(Ro1, Rd);
    
    [v1, v2, ob] = optimize(v, Rg, Rd, dx, dy, t, rx, ry, r, r1, rx1, ry1, R, vobs, vobs1, state, j, ob, RX, RY, RX1, RY1,vrx, vry, vrx1, vry1);
    v = [v1, v2];
    j = j+1;
    
    Rd(1) = Rd(1) + v(1)*t;
    Rd(2) = Rd(2) + v(2)*t;
    
    pos(1, i) = Rd(1);
    pos(2, i) = Rd(2);
    pos(3, i) = v(1);
    pos(4, i) = v(2);
    pos(5, i) = d;
    pos(6, i) = state;
    pos(9, i) = Ro(1);
    pos(10, i) = Ro(2);  
    pos(11, i) = Ro1(1);
    pos(12, i) = Ro1(2);
    
end    

%%
figure;
pause(1);
for i=1:length(pos)
xlim([-10,25]);
ylim([0,25]);


hold all;

h = scatter(pos(1, i), pos(2, i));
h1 = scatter(pos(9, i), pos(10, i));
h2 = scatter(pos(11, i), pos(12, i));


s = 1.7;
s1 = 1.35;
s2 = 1.35;

currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2)

currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);
markerWidth1 = s1/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h1, 'SizeData', markerWidth1^2)

currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);
markerWidth2 = s2/diff(xlim)*axpos(3); % Calculate Marker width in points
set(h2, 'SizeData', markerWidth2^2)

drawnow
pause(0.01);
clf;

end

function [d, dx, dy] = distance(X, Y)
    dx = (X(1) - Y(1));
    dy = (X(2) - Y(2));
    d = sqrt((dx)^2 + (dy)^2);
end

function [vx, vy, ob] = optimize(v, Rg, Rd,dx, dy, t, rx, ry, r, r1, rx1, ry1, R, vobs, vobs1, state, j, ob, RX, RY, RX1, RY1,vrx, vry, vrx1, vry1)
    objective = @(x) ((Rg(1) - (Rd(1) + x(1)*t))^2 + (Rg(2) - (Rd(2) + x(2)*t))^2);
       
    x0 = [v(1), v(2)];
    
    disp(['Initial Objective: ' num2str(objective(x0))])
    
    ob(1, j) = objective(x0);
        
    lb = [-3, -3];
    ub = [5, 5];  
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,@(x)nlcon_Vo_2(x,R,dx, dy, v,rx,ry,r, r1, rx1, ry1, vobs, vobs1, state,  RX, RY, RX1, RY1,vrx, vry, vrx1, vry1));
      
    vx = x(1);
    vy = x(2);
    
end
