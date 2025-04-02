% === runBridge.m: Driver Script ===
clear; clc;

% Parameters
n = 6;          % number of diamond cells
L = 0.2;        % horizontal spacing between columns
EA_bar = 10;    % axial stiffness for bars (type=1)
EA_cable = 1;   % axial stiffness for cables (type=2)

% Generate bridge geometry using our modified pattern.
[nodes_xy, el, n_fix] = genBridge(n, L);
EA_el = [EA_bar, EA_cable];

% Choose the load application point.
% In our design, the full interior column is at physical column 2.
% In MATLAB indexing, column 3 corresponds to physical column 2 (since col 1 = physical 0).
centerCol = 3;  
% For a full column, the bottom node is given by pos = 2.
% (Using our helper, nodeIndex(centerCol,2) would give the global node; here we know
% that in column 3 the starting index is computed in genBridge. In this driver we do not
% have direct access to nodeIndex, so we determine it from the known pattern.)
%
% Let’s reassemble the node numbering from the pattern:
%   Column 0 (physical 0, full): nodes 1–3.
%   Column 1 (physical 1, reduced): nodes 4–5.
%   Column 2 (physical 2, full): nodes 6–8.
%   Column 3 (physical 3, reduced): nodes 9–10.
%   Column 4 (physical 4, reduced): nodes 11–12.
%   Column 5 (physical 5, full): nodes 13–15.
%
% We choose the full interior column (physical column 2), which in our numbering is Column 2 (MATLAB index 3).
% Its bottom node is: for a full column, node = colStart(3) + 2.
% From genBridge, column 1 started at 1 (3 nodes), column 2 (MATLAB index 2) is reduced (2 nodes),
% so column 3 (MATLAB index 3) is full. Its starting index can be computed as:
%   colStart(1)=1; (full: 3 nodes) → colStart(2)=4;
%   col 2 is reduced → colStart(3)=4+2 = 6.
% Thus, the bottom node in column 3 is node 6 + 2 = 8.
bottomCenterNode = 8;
P = [bottomCenterNode, 2, -1];  % Apply a load of –1 in the y–direction

% Solve the truss (using your provided solveTruss function)
[ux, uy, Rx, Ry, Fint] = solveTruss(nodes_xy, el, EA_el, P, n_fix);

% Draw the color-coded truss.
figure; 
drawTruss(nodes_xy, el, Fint);
title('5-Diamond Bridge with Alternating Middle Nodes Removed');

% Identify which element carries the highest tension and compression.
[maxTension, eTension] = max(Fint);
[minCompression, eCompression] = min(Fint);
fprintf('Highest Tension: Element %d, Force = %.3f\n', eTension, maxTension);
fprintf('Highest Compression: Element %d, Force = %.3f\n', eCompression, minCompression);
