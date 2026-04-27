clear
clearvars
close all

%% Spaceship Layout Properties

innerRadius = 10;   % km
outerRadius = 20;   % km
numRadial = 1;      % number of radial divisions
numTheta  = 20;     % number of angular divisions

%% Material Properties

sigma = 5.67 .* 10 .^ -8; % Stephan-Boltzman constant
epsilon = 0.5;            % Emissivity of material 
k = 21.9;                 % thermal conductivity (W/(m*K))

%% Heating Properties

Twall = 300;                          % Outer nodes have a known temperature
Q  = 10000;                           % W/m^3
qIn = 10;                             % heat flux in from nuclear plant (W/m^2)
qOut = sigma .* epsilon .* Twall ^ 4; % heat flux out from radiation (W/m^2)

%% Plotting

% Angular divisions 
theta = linspace(0, 2*pi, numTheta+1);

% Radial divisions
r2 = linspace(innerRadius^2, outerRadius^2, numRadial+1);
r = sqrt(r2);

% Create a mesh based on the following elements
[Theta, R] = meshgrid(theta, r);
X = R .* cos(Theta);
Y = R .* sin(Theta);

% Plot
figure;
hold on;
axis equal;

% Plot radial lines
for i = 1:numTheta+1
    plot(X(:,i), Y(:,i), 'k');
end

% Plot circular arcs
for j = 1:numRadial+1
    plot(X(j,:), Y(j,:), 'b');
end

xlabel('X (km)');
ylabel('Y (km)');
title('Equal-Area Mesh');

hold off;