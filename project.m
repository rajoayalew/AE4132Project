clear
clearvars
close all

%% Spaceship Layout Properties

innerRadius = 10;   % km
outerRadius = 20;   % km
numRadial = 20;     % number of radial divisions
numTheta  = 20;     % number of angular divisions
t = 1;              % thickness (km)

%% Material Properties

sigma = 5.67 .* 10 .^ -8; % Stephan-Boltzman constant
epsilon = 0.5;            % Emissivity of material 
k = 21.9;                 % thermal conductivity (W/(m*K))

%% Heat Properties

Tinner = 10000;                        % Temperature of nodes facing nuclear reactor (r = 1 aka inner ring) 
Touter = 300;                          % Temperature of nodes facing outer space (r = numRadial aka outer ring)
Q  = 10000;                            % W/m^3
qIn = 10;                              % heat flux in from nuclear plant (W/m^2)
qOut = sigma .* epsilon .* Touter ^ 4; % heat flux out from radiation (W/m^2)

%% Graph Properties 

offset = 0.2;                          % Offset for the node number in the plots (visual effect only)

%% Mesh Generation

% Angular divisions 
theta = linspace(0, 2*pi, numTheta+1);
theta = theta(1:end-1);

% Radial divisions
r = sqrt(linspace(innerRadius^2, outerRadius^2, numRadial+1));

% Create a mesh based on the following elements
[Theta, R] = meshgrid(theta, r);
X = R .* cos(Theta);
Y = R .* sin(Theta);

% nodes is a (N x 2) matrix of node coordinates in Cartesian coords
% For a given radial station i and theta location j, nodeID(i, j) returns 
% the node number that can be used to index the nodes table 
nodes = [Y(:), X(:)];                 
nodeID = reshape(1:numel(X), size(X)); 

% Element Table
numNodes = size(nodes,1);           % Number of nodes in the mesh
numElements = numRadial * numTheta; % Number of elements in the mesh
elements = zeros(numElements, 5);   

%{
Each row of the element matrix corresponds to an element in the mesh. Each
column corresponds to a node in the quadliteral element. Figure 1.

n2 (i+1,j) ---- n3 (i,j+1)
     |              |
     |              |
n1 (i,j) -----  n4 (i+1,j+1)

The fifth element indicates if the element has a boundary heat flux applied to it
%}
currElemNumber = 1;
for i = 1:numRadial
    for j = 1:numTheta

        jp1 = j + 1;
        if (jp1 > numTheta)
            jp1 = 1;
        end

        n1 = nodeID(i, j);      % node 1 is bottom left
        n2 = nodeID(i+1, j);    % node 2 is top left
        n3 = nodeID(i+1, jp1);  % node 3 is top right
        n4 = nodeID(i, jp1);    % node 4 is bottom right 

        %{
        1 = inner ring (imposed heat flux in)
        2 = outer ring (imposed heat flux out)
        3 = in between inner and outer (no imposed heat flux)
        %}

        if (i == 1)
            % All elements at radial station 1 have an imposed heat flux into 
            % the element from the reactor.
            elements(currElemNumber,:) = [n1, n2, n3, n4, 1];
        elseif (i == numRadial)
            % All elements at radial station numRadial have an heat flux leaving 
            % the element and heading into space via radiation. 
            elements(currElemNumber,:) = [n1, n2, n3, n4, 2];
        else
            % All other elements have no imposed heat flux
            elements(currElemNumber,:) = [n1, n2, n3, n4, 3];
        end

        currElemNumber = currElemNumber + 1;
    end
end

%% Isoparameteric Mapping

K = zeros(numNodes, numNodes); % Global stiffness matrix
R = zeros(numNodes, 1);        % Global load vector

% Shape Functions

% N_1(zeta, eta) = 1/4 * (1 - zeta) * (1 - eta) = 1/4 * (1 - zeta - eta + zeta*eta)
% N_2(zeta, eta) = 1/4 * (1 - zeta) * (1 + eta) = 1/4 * (1 - zeta + eta - zeta*eta)
% N_3(zeta, eta) = 1/4 * (1 + zeta) * (1 + eta) = 1/4 * (1 + zeta + eta + zeta*eta)
% N_4(zeta, eta) = 1/4 * (1 + zeta) * (1 - eta) = 1/4 * (1 + zeta - eta - zeta*eta)

% Lambda is our conductivity matrix. The first element and the last element typically
% indicate the conductivity in the x and y directions, but in this case they are the same
% as the material is the same throughout the structure
lambda = [k, 0
          0, k];

% This gives the gauss points for a 4-node quad element. Each row corresponds to a particular node
% number. For example, node 1 while being located in reality at zeta = -1, eta = -1, will be located
% at zeta = -1/sqrt(3) and eta = -1/sqrt(3) for the purposes of guass-quadrature intergration. Node
% 2 is located at zeta = -1, eta = 1 and therefore the gauss point will be at zeta = -1/sqrt(3) and
% eta = 1/sqrt(3). The local node number, as shown in the Figure 1 above, is the index into this
% table. Since we are using a 2-point Gauss quadrature rule the weight is equal to 1.
%{
Natural coordinate space (zeta, eta): 
Figure 2
   η
   ↑
   |   Node 2 (-1,+1)        Node 3 (+1,+1)
   |        ●----------------●
   |        |                |
   |        |    •     •     |   ← Gauss points (±0.577)
   |        |                |
   |        |    •     •     |
   |        |                |
   |        ●----------------●
   |   Node 1 (-1,-1)        Node 4 (+1,-1)
   |
   +------------------------------→ ζ
%}
gaussTable2D = [-1/sqrt(3), -1/sqrt(3)
                -1/sqrt(3),  1/sqrt(3)
                 1/sqrt(3),  1/sqrt(3)
                 1/sqrt(3), -1/sqrt(3)];

gaussTable1D = [-1/sqrt(3), 1/sqrt(3)];

for currElemNumber = 1:numElements

    % Get the global coordinates of the nodes for a given element

    % currElemeCoords is a 4x2 matrix where the global coordinates
    % for node 1 is loacted in the first row, node 2's global coordinates
    % are located in the second row, and so on and so forth
    currElemCoords = nodes(elements(currElemNumber, 1:4), :);

    % Element stiffness matrix
    Ke = zeros(4,4);

    % Local load vector
    Re = zeros(4, 1);

    for gaussPoint = 1:4
        % Gauss coords in local coord system
        zeta = gaussTable2D(gaussPoint, 1);    
        eta = gaussTable2D(gaussPoint, 2);

        % NMat is the shape function result for a given gauss point 
        NMat = shape_function_calc(zeta, eta);

        % Get the partial derivatives at the gauss points for each of the shape functions 
        % dN_dzeta is a 4x1 matrix containing \frac{\partial N_i(zeta, eta)}{\partial zeta}
        % dN_deta is a 4x1 matrix containing \frac{\partial N_i(zeta, eta)}{\partial eta}
        [dN_dzeta, dN_deta] = shape_function_derivatives(zeta, eta);

        % Get the Jacobian 
        [J, detJ] = jacobian(currElemCoords, dN_dzeta, dN_deta);

        % Get the partial derivatives of the shape functions with respect to x and y
        %{
        [dN_dzeta'; dN_deta'] is equal to the following:
        
        \begin
        %}
        dN_dxdy = J \ [dN_dzeta'; dN_deta'];

        % Assemble B matrix
        dN_dx = dN_dxdy(1,:);
        dN_dy = dN_dxdy(2,:);
        B = [dN_dx; dN_dy];

        % Calculate stiffness matrix
        % Here we are doing gauss quadrature integration but you will notice there is no
        % weights here and that is because for a 2-point gauss quad integration
        % the weight is equal to 1
        Ke = Ke + (B' * lambda * B * detJ * t);

        % Calculate local load vector due to heat generation by the element
        % We use the shape functions to determine how much of the distributed heat
        % source should be assigned to each node 
        Re = Re + (Q * NMat * detJ * t);
    end

    % Check if we need to impose heat fluxes on this element
    nuclearFluxCheck = (elements(currElemNumber, 5) == 1);
    radiationFluxCheck = (elements(currElemNumber, 5) == 2);

    if (nuclearFluxCheck)
        % Applied on the eta = -1 boundary. Eta = -1 corresponds with the line
        % of an element that faces the reactor. Similarily to above we are going to use Gauss
        % quadrature in order to intergrate the following:
        % $$\{r\}_{q}=\int_{-1}^{1} q_{in}^{*}\{N(\zeta, \eta=-1)\}\sqrt{ \left( \frac{\partial x}{\partial \zeta} 
        % \right)^{2}+\left( \frac{\partial y}{\partial \zeta} \right)^{2} } \, d \zeta $$

        % Look to Chapter 5 Pt 2 Pages 18 and 19 for more info

        for gaussPoint = 1:2 

            % zeta first iter = -1/sqrt(3) and zeta second iter = 1/sqrt(3)
            zeta = gaussTable1D(gaussPoint); 
            eta = -1;

            % Calculate the shape functions along the eta = -1 boundary
            NMat = shape_function_calc(zeta, eta);

            % Get the partial derivatives at the gauss points for each of the shape functions 
            % dN_dzeta is a 4x1 matrix containing \frac{\partial N_i(zeta, eta)}{\partial zeta}
            % dN_deta is a 4x1 matrix containing \frac{\partial N_i(zeta, eta)}{\partial eta}
            [dN_dzeta, dN_deta] = shape_function_derivatives(zeta, eta);

            % Get the Jacobian in order to get the partial derivatives with respect to zeta
            [J, ~] = jacobian(currElemCoords, dN_dzeta, dN_deta);

            % Calculate the local load vector and preform gauss-quadrature intergration
            dx_dzeta_sq = J(1, 1) .^ 2;
            dy_dzeta_sq = J(2, 1) .^ 2;
            Re = Re + (qIn * NMat * sqrt(dx_dzeta_sq + dy_dzeta_sq) * t);
        end
    end

    if (radiationFluxCheck)
        % Applied on the eta = 1 boundary. Eta = 1 corresponds with the line
        % of an element that faces outer space. Similarily to above we are going to use Gauss
        % quadrature in order to intergrate the following:
        % $$\{r\}_{q}=\int_{-1}^{1} q_{out}^{*}\{N(\zeta, \eta=1)\}\sqrt{ \left( \frac{\partial x}{\partial \zeta} 
        % \right)^{2}+\left( \frac{\partial y}{\partial \zeta} \right)^{2} } \, d \zeta $$

        % Applied on the zeta = 1 boundary
        for gaussPoint = 1:2
            % zeta first iter = -1/sqrt(3) and zeta second iter = 1/sqrt(3)
            zeta = gaussTable1D(gaussPoint); 
            eta = 1;

            % Calculate the shape functions along the eta = -1 boundary
            NMat = shape_function_calc(zeta, eta);

            % Get the partial derivatives at the gauss points for each of the shape functions 
            % dN_dzeta is a 4x1 matrix containing \frac{\partial N_i(zeta, eta)}{\partial zeta}
            % dN_deta is a 4x1 matrix containing \frac{\partial N_i(zeta, eta)}{\partial eta}
            [dN_dzeta, dN_deta] = shape_function_derivatives(zeta, eta);

            % Get the Jacobian in order to get the partial derivatives with respect to zeta
            [J, ~] = jacobian(currElemCoords, dN_dzeta, dN_deta);

            % Calculate the local load vector and preform gauss-quadrature intergration
            dx_dzeta_sq = J(1, 1) .^ 2;
            dy_dzeta_sq = J(2, 1) .^ 2;
            Re = Re + (qOut * NMat * sqrt(dx_dzeta_sq + dy_dzeta_sq) * t);
        end
    end

    % Get the number of each node for this element in the global system 
    % Add its contribution to the global stiffness matrix for its particular nodes
    % Add its contribution to the global load vector for its particular nodes
    globalNodes = elements(currElemNumber,1:4);
    K(globalNodes, globalNodes) = K(globalNodes, globalNodes) + Ke;
    R(globalNodes) = R(globalNodes) + Re;
end

%% Enforcing Temperature Boundary Conditions



%% Plotting 

% Element Numbering 
figure;
hold on;
axis equal;
grid on;

numElements = size(elements,1);

for e = 1:numElements
    
    elemNodes = elements(e, 1:4);
    coords = nodes(elemNodes, :);   % 4×2 matrix
    
    % Close the quad for plotting
    coords_closed = [coords; coords(1,:)];
    
    % Draw element
    plot(coords_closed(:,1), coords_closed(:,2), 'k-');
    
    % Compute centroid
    centroid = mean(coords, 1);
    
    % Label element number
    text(centroid(1), centroid(2), sprintf('%d', e), ...
        'Color','r', 'FontSize',10, 'HorizontalAlignment','center');
end

for n = 1:size(nodes,1)
    x = nodes(n,1);
    y = nodes(n,2);

    text(x + offset, y + offset, sprintf('%d', n), ...
        'Color','b', ...
        'FontSize',8, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
end

xlabel('X');
ylabel('Y');
title('Mesh with Global Element and Nodes');
hold off

%% Functions

function result = shape_function_calc(zeta, eta)

    % This function should be pretty self explanatory as it simply just calculates the shape
    % functions at a particular zeta and eta.

    % N_1(zeta, eta) = 1/4 * (1 - zeta) * (1 - eta) = 1/4 * (1 - zeta - eta + zeta*eta)
    % N_2(zeta, eta) = 1/4 * (1 - zeta) * (1 + eta) = 1/4 * (1 - zeta + eta - zeta*eta)
    % N_3(zeta, eta) = 1/4 * (1 + zeta) * (1 + eta) = 1/4 * (1 + zeta + eta + zeta*eta)
    % N_4(zeta, eta) = 1/4 * (1 + zeta) * (1 - eta) = 1/4 * (1 + zeta - eta - zeta*eta)

    N1 = 0.25 * (1 - zeta) * (1 - eta);
    N2 = 0.25 * (1 - zeta) * (1 + eta);
    N3 = 0.25 * (1 + zeta) * (1 + eta);
    N4 = 0.25 * (1 + zeta) * (1 - eta);
    
    result = [N1; N2; N3; N4];
end

function [dN_dzeta, dN_deta] = shape_function_derivatives(zeta, eta)

    % Above there is a section called shape functions which are the shape functions for each node in
    % the 4-node quad. All I did in this section was take the derivative of each of those shape
    % functions with respect to zeta and eta. Each index of dN_dzeta and dN_deta corresponds to the
    % partial derivatives with respect to zeta and eta respectively. The first row is for node 1,
    % the second row is for node 2, and so on and so forth. Look to Figure 1 to determine what nodes
    % correspond with what location on the quad.

    dN_dzeta = 1/4 * [
        -(1 - eta)
        -(1 + eta)
         (1 + eta)
         (1 - eta)
    ];

    dN_deta = 1/4 * [
        -(1 - zeta)
         (1 - zeta)
         (1 + zeta)
        -(1 + zeta)
    ];
end

function [J, detJ] = jacobian(coords, dN_dzeta, dN_deta)
    % coords: 4x2 matrix [x y]
    J = zeros(2,2);
    
    %{
    The Jacobian matrix is equal to:
    J = [dx/dzeta, dx/deta; dy/dzeta, dy/deta] note that these are
    all partial derivatives. If you remember, we can get the x and y loacations from
    a given zeta and eta by using the mapping equation. 

    x(zeta, eta) = \Sum_{i=1}^{M}x_i * N_i(zeta,eta) and y(zeta, eta)
    y(zeta, eta) = \Sum_{i=1}^{M}y_i * N_i(zeta,eta)

    Where x_i and y_i are the the global x and y locations of the node that N_i corresponds to and M
    = 4 due to us using a 4-node quad.
   
    We can also get the partial derivates (like dx/dzeta and dy/deta). We do this by multiplying not
    by the normal shape functions but their derivatives instead. For example, \frac{\partial
    x}{\partial zeta} = \Sigma_{i=1}^{M}x_i * \frac{\partial N_i(zeta, eta)}{\partial zeta}. We get
    these partial derivatives of these shape functions from the shape_function_derivatives function
    above.
    %}

    for i = 1:4
        J(1,1) = J(1,1) + dN_dzeta(i) * coords(i,1); % dx/dzeta
        J(1,2) = J(1,2) + dN_deta(i)  * coords(i,1); % dx/deta
        J(2,1) = J(2,1) + dN_dzeta(i) * coords(i,2); % dy/dzeta
        J(2,2) = J(2,2) + dN_deta(i)  * coords(i,2); % dy/deta
    end
    
    detJ = det(J);
end