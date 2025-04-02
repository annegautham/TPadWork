
function [ux, uy, Rx, Ry, Fint] = solveTruss(nodes_xy, el,...
                                             EA_el, P, n_fix)
% Calculate the nodal displacements and reaction forces for a truss
% composed of two different types of elements.
% 
% Description of input arguments:
% 
% nodes_xy: table of x, y positions of nodes in the following format
%   [n1_x n1_y
%    n2_x n2_y
%    n3_x n3_y ...]
% 
% el: table of truss elements specifying the element type and the
%     numbers of the nodes it is connected to. 'type' is either 1 (for a
%     strut) or 2 (for a cable).
%   [type_el_1 node_1_el_1 node_2_el_1
%    type_el_2 node_1_el_2 node_2_el_2
%    type_el_3 node_1_el_3 node_2_el_3 ...]
% 
% EA_el: vector of axial stiffnesses for each element type.
%   [EA_eltype_1, EA_eltype_2, ...]
% 
% P: table of applied forces. 'direction' is either 1 for x or 2 for y.
%    'magnitude' may be either positive or negative.
%   [node_1 direction_1 magnitude_1
%    node_2 direction_2 magnitude_2
%    node_3 direction_3 magnitude_3 ...]
% 
% n_fix: table of displacement boundary conditions.
%   [node_1 direction_1 displacement_1
%    node_2 direction_2 displacement_2
%    node_3 direction_3 displacement_3 ...]
% 
% Description of outputs:
% 
% ux: : vector of x-displacements for each node
% uy: : vector of y-displacements for each node
% Rx: : vector of reaction forces in the x-direction for each node (entries
%       are zero for unconstrained nodes)
% Ry  : vector of reaction forces in the y-direction
% Fint: vector of internal forces for each element

[n_n, ~] = size(nodes_xy);  % Number of nodes
[n_el, ~] = size(el);       % Number of struts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Element stiffness matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create global stiffness matrix
K_global = zeros(2*n_n);

% Iterate over all elements and calculate stiffness matrix for each
for e = 1:n_el

    n1 = el(e, 2);  % index of first element node
    n2 = el(e, 3);  % index of second element node
    
    % Degrees of freedom associated with this element
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    
    % Get element angle and length
    [theta, l] = cart2pol(nodes_xy(n2, 1) - nodes_xy(n1, 1), ...
                      nodes_xy(n2, 2) - nodes_xy(n1, 2));

    % Determine the element type to determine the axial stiffness
    if el(e, 1) > length(EA_el)
        error('Unrecognized element type: element number: %d', e);
    end
    
    % set stiffness based on element type
    EA = EA_el(el(e, 1));
    
    %%% FINISH THE CODE BELOW:
    %%% Create a 4x4 stiffness matrix for the current element
    %%% the angle and length have already been calculated for you above
    % directional cosines
    c = cos(theta);
    s = sin(theta);
    k = EA/l;
    k_el = k * [ c^2,  c*s, -c^2, -c*s;
             c*s,  s^2, -c*s, -s^2;
            -c^2, -c*s,  c^2,  c*s;
            -c*s, -s^2,  c*s,  s^2 ];
    
    % Add the stiffness matrix to the appropriate location within the
    % global stiffness matrix
    K_global(dofs, dofs) = K_global(dofs, dofs) + k_el;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forces and boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% force vector
f = zeros(2*n_n, 1);
f(2*(P(:,1) - 1) + P(:,2)) = P(:,3);

% Apply boundary conditions
K_aug = K_global;

[n_bc, ~] = size(n_fix);

for i = 1:n_bc
    
    dof = 2*(n_fix(i, 1) - 1) + n_fix(i, 2);
    
    K_aug(dof, :) = 0;
    K_aug(dof, dof) = 1;
    f(dof) = n_fix(i, 3);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate solution from matrix equation K*u = f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOLVE
u = K_aug \ f;

ux = u(1:2:end);
uy = u(2:2:end);

% Get reaction forces
f_r = K_global * u - f;
Rx = f_r(1:2:end);
Ry = f_r(2:2:end);

% Get member internal forces
Fint = zeros(n_el, 1);
for e = 1:n_el
    
    n1 = el(e, 2);
    n2 = el(e, 3);
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    
    % Get element angle and length
    [theta, l] = cart2pol(nodes_xy(n2, 1) - nodes_xy(n1, 1), ...
                      nodes_xy(n2, 2) - nodes_xy(n1, 2));
    
    % Determine the element type to determine the axial stiffness
    if el(e, 1) > length(EA_el)
        error('Unrecognized element type: element number: %d', e);
    end
    
    % set stiffness based on element type
    EA = EA_el(el(e, 1));

    k = EA/l;  % axial stiffness
    
    Fint(e) = k * ((u(dofs(3)) - u(dofs(1)))*cos(theta) + ...
                   (u(dofs(4)) - u(dofs(2)))*sin(theta));
    
end

end
