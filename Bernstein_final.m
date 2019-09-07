syms tc     % for differenciation
syms t

Ro = [5,0];      %initial position
Rg = [15, 15];    %final postion
Rd_inter1 = [10, 25];     %intermediate position

ko = 0;
kf = 0;
to = 0;     %intial time
tf = 25;    %final time
t_inter1 = 12.5;   %time at which drone should reach R_inter
dt = 0.1;  %time interval
n = 5;
i = 0;      % 0 <= i <= 4 

pos = zeros(2, 300);
t_current = 0;
q = 1;

wx0 = Ro(1);
wx5 = Rg(1);

a11 = bernstein_t(to, tf, t_inter1, n, 1, t); 
a12 = bernstein_t(to, tf, t_inter1, n, 2, t);
a13 = bernstein_t(to, tf, t_inter1, n, 3, t);
a14 = bernstein_t(to, tf, t_inter1, n, 4, t);
a21 = bernstein_t(to, tf, tf, n, 1, t);
a22 = bernstein_t(to, tf, tf, n, 2, t);
a23 = bernstein_t(to, tf, tf, n, 3, t);
a24 = bernstein_t(to, tf, tf, n, 4, t);
A = [a11, a12, a13, a14; a21, a22, a23, a24];
b11 = Rd_inter1(1) - wx0*bernstein_t(to, tf, t_inter1, n, 0, t) - wx5*bernstein_t(to, tf, t_inter1, n, 5, t);
b21 = Rg(1) - wx0*bernstein_t(to, tf, tf, n, 0, t) - wx5*bernstein_t(to, tf, tf, n, 5, t);

B = [b11;b21];
W = pinv(A)*B;
    %disp(W);
Xt = X_t(to, tf, tf, n, t, wx0, wx5, W);
disp(Xt);
    
wk0 = ko;
wk5 = kf;

ka11 = Fi(to, tf, tf, n, 1, t, wx0, wx5, W);
ka12 = Fi(to, tf, tf, n, 2, t, wx0, wx5, W);
ka13 = Fi(to, tf, tf, n, 3, t, wx0, wx5, W);
ka14 = Fi(to, tf, tf, n, 4, t, wx0, wx5, W);
ka21 = Fi(to, tf, t_inter1, n, 1, t, wx0, wx5, W);
ka22 = Fi(to, tf, t_inter1, n, 2, t, wx0, wx5, W);
ka23 = Fi(to, tf, t_inter1, n, 3, t, wx0, wx5, W);
ka24 = Fi(to, tf, t_inter1, n, 4, t, wx0, wx5, W);
ka31 = bernstein_t(to, tf, tf, n, 1, t);
ka32 = bernstein_t(to, tf, tf, n, 2, t);
ka33 = bernstein_t(to, tf, tf, n, 3, t);
ka34 = bernstein_t(to, tf, tf, n, 4, t);
ka41 = bernstein_t(to, tf, to, n, 1, t);
ka42 = bernstein_t(to, tf, to, n, 2, t);
ka43 = bernstein_t(to, tf, to, n, 3, t);
ka44 = bernstein_t(to, tf, to, n, 4, t);

kA = [ka11, ka12, ka13, ka14; ka21, ka22, ka23, ka24; ka31, ka32, ka33, ka34; ka41, ka42, ka43, ka44];

kb11 = Rg(2) - wk0*Fi(to, tf, tf, n, 0, t, wx0, wx5, W) - wk5*Fi(to, tf, tf, n, 5, t, wx0, wx5, W);
kb21 = Rd_inter1(2) - wk0*Fi(to, tf, t_inter1, n, 0, t, wx0, wx5, W) - wk5*Fi(to, tf, t_inter1, n, 5, t, wx0, wx5, W);
kb31 = kf - wk0*bernstein_t(to, tf, tf, n, 0, t) - wk5*bernstein_t(to, tf, tf, n, 5, t);
kb41 = ko  - wk0*bernstein_t(to, tf, to, n, 0, t) - wk5*bernstein_t(to, tf, to, n, 5, t);

kB = [kb11; kb21; kb31; kb41];

kW = pinv(kA)*kB;


while t_current<= tf
    disp(t_current);
    tx = t_current;
    Xt = X_t(to, tf, tx, n, t, wx0, wx5, W);
    pos(1, q) = Xt;
      
    ty = t_current;
    F0 = Fi(to, tf, ty, n, 0, t, wx0, wx5, W);
    F1 = Fi(to, tf, ty, n, 1, t, wx0, wx5, W);
    F2 = Fi(to, tf, ty, n, 2, t, wx0, wx5, W);
    F3 = Fi(to, tf, ty, n, 3, t, wx0, wx5, W);
    F4 = Fi(to, tf, ty, n, 4, t, wx0, wx5, W);
    F5 = Fi(to, tf, ty, n, 5, t, wx0, wx5, W);
    Yt = wk0*F0 + kW(1)*F1 + kW(2)*F2 + kW(3)*F3 + kW(4)*F4 + wk5*F5;
    pos(2, q) = Yt;
    
    t_current = t_current + dt;
    q = q + 1;
    
end


disp(pos);
figure(5);
plot(pos(1, :), pos(2, :), 'o');

    % bernstein polynomial

function F = Fi(to, tf, t_current, n, i, t, wx0, wx5, W)
    [~, ~, ~, dBu0] = bernstein_t(to, tf, t_current, n, 0, t);
    [~, ~, ~, dBu1] = bernstein_t(to, tf, t_current, n, 1, t);
    [~, ~, ~, dBu2] = bernstein_t(to, tf, t_current, n, 2, t);
    [~, ~, ~, dBu3] = bernstein_t(to, tf, t_current, n, 3, t);
    [~, ~, ~, dBu4] = bernstein_t(to, tf, t_current, n, 4, t);
    [~, ~, ~, dBu5] = bernstein_t(to, tf, t_current, n, 5, t);
    B = wx0*dBu0 + W(1)*dBu1 + W(2)*dBu2 + W(3)*dBu3 + W(4)*dBu4 + wx5*dBu5;
    [~, ~, Bui, ~] = bernstein_t(to, tf, t_current, n, i, t);
    fi = Bui*(B);
    F = int(fi,t, to, t_current);
end

function [X, Xv, Xt, Xtv] = X_t(to, tf, t_current, n, t, wx0, wx5, W)
    [X0, dX0, Xt0, dXt0] = bernstein_t(to, tf, t_current, n, 0, t);
    [X1, dX1, Xt1, dXt1] = bernstein_t(to, tf, t_current, n, 1, t);
    [X2, dX2, Xt2, dXt2] = bernstein_t(to, tf, t_current, n, 2, t);
    [X3, dX3, Xt3, dXt3] = bernstein_t(to, tf, t_current, n, 3, t);
    [X4, dX4, Xt4, dXt4] = bernstein_t(to, tf, t_current, n, 4, t);
    [X5, dX5, Xt5, dXt5] = bernstein_t(to, tf, t_current, n, 5, t);
    Xt = wx0*Xt0 + W(1)*Xt1 + W(2)*Xt2 + W(3)*Xt3 + W(4)*Xt4 + wx5*Xt5;
    Xtv = wx0*dXt0 + W(1)*dXt1 + W(2)*dXt2 + W(3)*dXt3 + W(4)*dXt4 + wx5*dXt5;
    X = wx0*X0 + W(1)*X1 + W(2)*X2 + W(3)*X3 + W(4)*X4 + wx5*X5;
    Xv = wx0*dX0 + W(1)*dX1 + W(2)*dX2 + W(3)*dX3 + W(4)*dX4 + wx5*dX5;
end

function [Bi, dBi, Bu, dBu] = bernstein_t(to, tf, t_current, n, i, t)
    u = (t - to)/(tf - to);

    Bu = nchoosek(n,i)*((1-u)^(n-i))*(u^i);
    dBu = diff(Bu);

    Bi = vpa(subs(Bu, t, t_current));
    dBi = vpa(subs(dBu, t, t_current));
end

