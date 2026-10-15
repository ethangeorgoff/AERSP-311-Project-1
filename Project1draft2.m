%Load dataset
clear; clc; close all;
data = load ("Project1_data.mat");
g = data.grid;
f = data.flow;

%Calculate the freestream conditions
c = g.chord;
muInf = 3.62e-7;
aInf = sqrt(1.4*f.R*f.Tinf);
Uinf = f.Mach*aInf;
rhoInf = f.Re*muInf/(Uinf*c);

%Select Outermost Boundary
x = g.x(:,end);
y = g.y(:,end);
u = f.u(:,end);
v = f.v(:,end);

%Calculate the circulation
dx = diff(x);
dy = diff(y);

uAvg = (u(1:end-1) + u(2:end))/2;
vAvg = (v(1:end-1) + v(2:end))/2;

dGamma = -(uAvg.*dx + vAvg.*dy);
Gamma = sum(dGamma);

%Calculate the lift coefficient and display results
cl = 2*Gamma/(Uinf*c);
clCFD = 0.9387;
percentDifference = 100*abs(cl-clCFD)/clCFD;

fprintf('Computed lift coefficient: %.4f\n', cl);
fprintf('Percent difference from CFD: %.2f%%\n', percentDifference);

%Plot the airfoil and velocity arrows
plot(g.x(:,1),g.y(:,1),'k-');
hold on;
plot(x,y,'k--');

idx = 1:8:(length(x)-1);
quiver(x(idx),y(idx),u(idx),v(idx),'b');

%Plot each segments circulation contribution
xMid = (x(1:end-1) + x(2:end))/2;

figure;
plot(xMid, dGamma, '.');

