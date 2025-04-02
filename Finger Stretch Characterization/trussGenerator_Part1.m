%Truss generator
%SampleCode for Part 1: Solve for reactions of truss with pin and roller
%support. Load of 1 force unit along the x-direction at node 2
%CIV 216 2023
%Project Part B

% nodes_xy: table of x, y positrx
% ions of nodes in the following format
%   [n1_x n1_y
%    n2_x n2_y
%    n3_x n3_y ...]

nodes_xy = [0 0;cos(pi*60/180) sin(pi*60/180);1 0];
%node 1, x_coordinate = 0, y_coordinate =0
% node 2, x_coordinate = lcos(60),y_coordinate = lsin(60)
%node 3, x_coordinate = l, y_coordinate = 0;

% el: table of truss elements specifying the element type and the
%     numbers of the nodes it is connected to. 'type' is either 1 (for a
%     strut) or 2 (for a cable).
% A strut is a compressive or tensile member. Cables only support tensile
% forces
%   [type_el_1 node_1_el_1 node_2_el_1
%    type_el_2 node_1_el_2 node_2_el_2
%    type_el_3 node_1_el_3 node_2_el_3 ...]

el = [1 1 2;1 2 3;1 1 3];
% element type = 1, node 1, node 2
% element type = 1, node 2, node 3
% element type = 1, node 1, node 3

% EA_el: vector of axial stiffnesses for each element type.
%   [EA_eltype_1, EA_eltype_2, ...]

EA_el = [1,1]; 
% [element axial stiffness EA = 1 for strut, element stiffness for cable = 1]. Second entry is not used here. 

% P: table of applied forces. 'direction' is either 1 for x or 2 for y.
%    'magnitude' may be either positive or negative.
%   [node_1 direction_1 magnitude_1
%    node_2 direction_2 magnitude_2
%    node_3 direction_3 magnitude_3 ...]

P=[2,1,1]; 
%[External force at node 2, x-component, magnitude of force = 1]

% n_fix: table of displacement boundary conditions.
%   [node_1 direction_1 displacement_1
%    node_2 direction_2 displacement_2
%    node_3 direction_3 displacement_3 ...]

n_fix=[1,1,0;1,2,0;3 2 0]; 

%[node 1, x-component, u_x = 0; node 1, y-component, u_y = 0; node 3;
%y-component, u_y = 0];

% Description of outputs:
% 
% ux: : vector of x-displacements for each node
% uy: : vector of y-displacements for each node
% Rx: : vector of reaction forces in the x-direction for each node (entries
%       are zero for unconstrained nodes)
% Ry  : vector of reaction forces in the y-direction
% Fint: vector of internal forces for each element

% For part 3 of the project, use your this truss generator code to generate
% the truss by modifying, nodes_xy,el,EA_el,P, and n_fix, according to the problem statement 

[ux, uy, Rx, Ry, Fint]=solveTruss(nodes_xy,el,EA_el,P,n_fix);

drawTruss(nodes_xy, el, Fint);
