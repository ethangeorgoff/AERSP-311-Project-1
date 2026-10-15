%Load dataset
clear; clc; close all;
data = load ("Project1_data.mat");
g = data.grid;                           %Contains the geometry
f = data.flow;                           %Contains the flow information(vel, Mach#, and temp)

%% Needed Equations/Givens
%c = g.chord;                     % ft
%Re = f.Re;
%M = f.Mach;
%Tinf = f.Tinf;                    % degrees Rankine
%R = f.R;                         % ft*lbf/(slug*degR)
%gamma = 1.4;
%muRef = 3.62e-7;                 % lbf*s/ft^2 = slug/(ft*s)
%Tref = 518.7;                    % degrees Rankine
%Suth = 198.72;                   % degrees Rankine
%clCFD = 0.9387;

%muInf = muRef*(Tinf/Tref)^(3/2)*(Tref + Suth)/(Tinf + Suth);
%aInf = sqrt(gamma*R*Tinf);        % ft/s
%Uinf = M*aInf;                   % ft/s (speed, not the x component)
%rhoInf = Re*muInf/(Uinf*c);       % slug/ft^3
%pInf = rhoInf*R*Tinf;            % lbf/ft^2, ideal-gas check
%qInf = 0.5*rhoInf*Uinf^2;         % lbf/ft^2

%Calculate the freestream conditions
c = g.chord;
muInf = 3.62e-7;
aInf = sqrt(1.4*f.R*f.Tinf);             %speed of sound
Uinf = f.Mach*aInf;                      %freestream speed
rhoInf = f.Re*muInf/(Uinf*c);            %freestream density

%Select Outermost Boundary
x = g.x(:,end);                %collects every position of x on outer boundary
y = g.y(:,end);                %collects every position of y on outer boundary
u = f.u(:,end);                %collects every position of u on outer boundary
v = f.v(:,end);                %collects every position of v on outer boundary

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
plot(g.x(:,1),g.y(:,1),'k-', 'LineWidth',2);
hold on;
plot(x,y,'k--');


ax = gca;
ax.Color = 'w';
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';


idx = 1:8:(length(x)-1);
quiver(x(idx),y(idx),u(idx),v(idx),'b', 'Linewidth', 1.25, 'MaxHeadSize',0.07);
axis equal;
grid on;
xlabel('x (ft)');
ylabel('y (ft)');
title('Airfoil Boundary and Velocity Field', 'Color','k');

% Legend creation and formatting
key = legend('Airfoil','Outermost Boundary','Velocity');

key.Color = 'w';
key.EdgeColor = 'k';
key.TextColor = 'k';

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


% Change colors of axes and background
ax = gca;
ax.Color = 'w';
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';

grid on;
xlabel('Segment midpoint x (ft)');
ylabel('\Delta\Gamma_i (ft^2/s)');
title('Circulation contribution of each segment', 'Color','k');


% Legend 2 creation and formatting
key2 = legend('Upper Boundary Branch','Lower Boundary Branch', 'Location','best');

key2.Color = 'w';
key2.EdgeColor = 'k';
key2.TextColor = 'k';


