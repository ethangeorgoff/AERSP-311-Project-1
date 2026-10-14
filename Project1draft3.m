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
fprintf('Freestream speed: %.4f ft/s\n', Uinf);
fprintf('Freestream density: %.8f slug/ft^3\n', rhoInf);

results = table(Gamma, cl, clCFD, percentDifference, 'VariableNames', {'Gamma_ft2_s','cl_KJ','cl_CFD','Difference_percent'});
disp(results);


%Plot the airfoil and velocity arrows
figure('Color','w');
plot(g.x(:,1),g.y(:,1),'k-');
hold on;
plot(x,y,'k--');

idx = 1:8:(length(x)-1);
quiver(x(idx),y(idx),u(idx),v(idx),'b');
axis equal;
grid on;
xlabel('x (ft)');
ylabel('y (ft)');
title('Airfoil Boundary and Velocity Field');
legend('Airfoil','Outermost Boundary','Velocity');

%Plot each segments circulation contribution
xMid = (x(1:end-1) + x(2:end))/2;

%Split the contour at its leftmost point
[~,iLE] = min(x);
upper = 1:(iLE-1);
lower = iLE:length(dGamma);

figure('Color','w');

plot(xMid(upper), dGamma(upper),'-','LineWidth',1.5);
hold on;
plot(xMid(lower), dGamma(lower),'-','LineWidth',1.5);

grid on;
xlabel('Segment midpoint x (ft)');
ylabel('\Delta\Gamma_i (ft^2/s)');
title('Circulation contribution of each segment');

legend('Upper Boundary Branch','Lower Boundary Branch', 'Location','best');

