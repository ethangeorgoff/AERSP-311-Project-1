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
quiver(x(idx),y(idx),u(idx),v(idx),'b', 'Linewidth', 1.01, 'MaxHeadSize',0.04);
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


% -------------- EXTRA CREDIT ------------------


% determine array size
[num_i, num_j] = size(g.x);


% compute gradient for all necessary values
[x_j, x_i] = gradient(g.x);
[y_j, y_i] = gradient(g.y);
[u_j, u_i] = gradient(f.u);
[v_j, v_i] = gradient(f.v);


% setup matrices for derivatives
dudx = zeros(num_i, num_j);
dudy = zeros(num_i, num_j);
dvdx = zeros(num_i, num_j);
dvdy = zeros(num_i, num_j);



% loop through to solve for real derivatives
for i = 1:num_i
    for j = 1:num_j
        % create matrix to later solve
        % [dx/di, dy/di]
        % [dx/dj, dy/dj]
        A = [ x_i(i,j), y_i(i,j); 
            x_j(i,j), y_j(i,j) ];

        % create vectors to solve the linear system
        % [du/di; du/dj] and [dv/di; dv/dj]
        vec_u = [ u_i(i,j); u_j(i,j) ];
        vec_v = [ v_i(i,j); v_j(i,j) ];

        % solve for derivatives of u using \
        sol_u = A \ vec_u;
        dudx(i,j) = sol_u(1);
        dudy(i,j) = sol_u(2);

        % solve for derivatives of v using \
        sol_v = A \ vec_v;
        dvdx(i,j) = sol_v(1);
        dvdy(i,j) = sol_v(2);
    end
end


% vorticity assignment
vorticity = dvdx - dudy;


% cell geometry and setup
cell_i = num_i - 1;
cell_j = num_j - 1;
cell_area = zeros(cell_i, cell_j);
cell_vorticity = zeros(cell_i, cell_j);

% loop to find cell areas
for i = 1:cell_i
    for j = 1:cell_j

        % find corners of quadrilateral cells
        x_corner = [g.x(i,j), g.x(i+1,j), g.x(i+1,j+1), g.x(i,j+1)];
        y_corner = [g.y(i,j), g.y(i+1,j), g.y(i+1,j+1), g.y(i,j+1)];

        % calculate cell area = 0.5 * ((x1y2 + x2y3 + x3y4 + x4y1) -
        % (y1x2 + y2x3 + y3x4 + y4x1))
        cell_area(i,j) = 0.5 * abs(x_corner(1)*y_corner(2) + x_corner(2)*y_corner(3)... % extend line
            + x_corner(3)*y_corner(4) + x_corner(4)*y_corner(1) ... % extend line
            - (y_corner(1)*x_corner(2) + y_corner(2)*x_corner(3) ... % extend line
            + y_corner(3)*x_corner(4) + y_corner(4)*x_corner(1)));
    


        % Average vorticity over the four cell corners
        cell_vorticity(i,j) = mean([vorticity(i,j), vorticity(i+1,j), ... % extend line
            vorticity(i+1,j+1), vorticity(i,j+1)]);

    end
end


% calculate values
gamma_vorticity = -sum(cell_vorticity .* cell_area, 'all');
cl_vorticity = 2 * gamma_vorticity / (Uinf * c);
percent_diff_2 = 100 * abs(cl_vorticity - clCFD)/ clCFD;


% table def
header = {'Path Integral', 'Vorticity Sum', 'CFD Solution'};
gamma_table = [Gamma; gamma_vorticity; NaN];
cl_table = [cl; cl_vorticity; clCFD];
percent_table = [percentDifference; percent_diff_2; 0];

extra_credit_table = table( ...
    ["Path Integral"; "Vorticity Sum"; "CFD Solution"], ...
    [Gamma; gamma_vorticity; NaN], ...
    [cl; cl_vorticity; clCFD], ...
    [percentDifference; percent_diff_2; 0], ...
    'VariableNames', {'Method', 'Gamma', 'cl', 'PercentDifference'});

disp('-------- Extra Credit ---------');
disp(extra_credit_table);


% plotting final results
figure('Color','w');
levels = linspace(-2000, 2000, 50);
contourf(g.x, g.y, vorticity, levels, 'LineStyle', 'none');
hold on;
fill(g.x(:,1), g.y(:,1), 'k'); % Solid black airfoil patch

cb = colorbar;
clim([-2000, 2000]);
ylabel(cb, 'Vorticity omega_z (1/s)', 'Color', 'k');
colormap(jet);
axis equal; grid on;
xlabel('x (ft)'); ylabel('y (ft)');
title('Vorticity Contour Field (omega_z)', 'Color', 'k');

axEC = gca;
axEC.Color = 'w';
axEC.XColor = 'k';
axEC.YColor = 'k';