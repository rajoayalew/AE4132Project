clear
clearvars
close all

%% Spaceship Layout Properties

innerRadius = 10;   % km
outerRadius = 20;   % km
numRadial = 10;     % number of radial divisions
numTheta  = 20;      % number of angular divisions

%% Material Properties

sigma = 5.67 .* 10 .^ -8; % Stephan-Boltzman constant
epsilon = 0.5;            % Emissivity of material 
k = 21.9;                 % thermal conductivity (W/(m*K))

%% Heat Properties

Tinner = 10000;                        % Inner nodes 
Touter = 300;                          % Outer nodes have a known temperature
Q  = 10000;                            % W/m^3
qIn = 10;                              % heat flux in from nuclear plant (W/m^2)
qOut = sigma .* epsilon .* Touter ^ 4; % heat flux out from radiation (W/m^2)

%% Mesh Creation

% Angular divisions 
theta = linspace(0, 2*pi, numTheta+1);

% Radial divisions
r = sqrt(linspace(innerRadius^2, outerRadius^2, numRadial+1));

% Create a mesh based on the following elements
[Theta, R] = meshgrid(theta, r);
X = R .* cos(Theta);
Y = R .* sin(Theta);

% Element Table
nodes = [Y(:), X(:)];   % (N x 2) table of node coordinates in Cartesian
nodeID = reshape(1:numel(X), size(X)); % For a given radial station i and theta location j, nodeID(i, j) returns the index in the nodes table 

numNodes = size(nodes,1);
numElements = numRadial * numTheta;

elements = zeros(numElements, 4); 


%% Plotting 

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
grid on, grid minor
title('Equal-Area Mesh');

hold off;