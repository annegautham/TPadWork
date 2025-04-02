function [nodes_xy, el, n_fix] = genBridgeAlt(n, L)
% genBridgeAlt  Generate a "bridge" truss with n diamond cells,
%   but some columns have no middle node (omitting them in an alternating way).
%
% In this example, we do the specific pattern:
%   - Columns = 0..n (so n+1 total columns => n "bays")
%   - For n=5, that means columns 0..5.
%   - We pin (fix x,y) all nodes at column 0 and column n.
%   - We choose which columns have a middle node, e.g.:
%       col0 => yes, col1 => no, col2 => yes, col3 => no,
%       col4 => no, col5 => yes.
%   - This matches a figure where node 5, 11, 14 do not exist,
%     yet columns 0,2,5 still have the middle node.
%
% INPUTS:
%   n : number of diamond cells (so there are n+1 columns)
%   L : horizontal spacing between adjacent columns
%
% OUTPUTS:
%   nodes_xy : (#nodes x 2) array of nodal coordinates
%   el       : element connectivity [type, node1, node2]
%              (type=1 => "bar", type=2 => "cable")
%   n_fix    : boundary conditions [node, dir, prescribed_disp]
%              pinned in x & y for the nodes in columns 0 and n.
%
%   You can adapt the pattern below if you want to change which
%   columns have a middle node vs. which do not.

    % --- 1) Decide which columns have a "middle" node ---
    % Here we hard-code a pattern for demonstration. For n=5, that
    % means columns i=0..5. Suppose we want:
    %   col0 => has middle
    %   col1 => no middle
    %   col2 => has middle
    %   col3 => no middle
    %   col4 => no middle
    %   col5 => has middle
    %
    % Of course, you can generalize this to any pattern you like.
    hasMiddle = false(1, n+1);  % default false
    % Always give column 0 a middle node:
    hasMiddle(1) = true;
    % Always give column n a middle node:
    hasMiddle(n+1) = true;
    % Manually specify the interior pattern (for n=5):
    % col1 => false, col2 => true, col3 => false, col4 => false
    if n==5
        hasMiddle(3) = true;  % col2 (index=2+1=3)
        % (col1, col3, col4 remain false)
    end
    
    % --- 2) Assign node numbers and coordinates ---
    % We will scan from column 0..n and build:
    %   top(i) = always
    %   mid(i) = only if hasMiddle(i)
    %   bot(i) = always
    %
    % This yields a variable total number of nodes. We'll store them
    % in arrays topID(i+1), midID(i+1), botID(i+1), with 0-based i
    % but 1-based indexing in the arrays. If a middle node does not
    % exist, midID(i+1) = 0.
    
    topID = zeros(1, n+1);
    midID = zeros(1, n+1);
    botID = zeros(1, n+1);
    
    % We'll store the (x,y) coordinates in a growing list
    nodeCount = 0;
    nodes_xy = [];  % will grow row by row
    
    for i = 0:n
        x_i = i*L;
        
        % Top node
        nodeCount = nodeCount + 1;
        topID(i+1) = nodeCount;
        nodes_xy(nodeCount, :) = [x_i, 2*L];
        
        % Middle node if hasMiddle(i+1)
        if hasMiddle(i+1)
            nodeCount = nodeCount + 1;
            midID(i+1) = nodeCount;
            nodes_xy(nodeCount, :) = [x_i, L];
        else
            midID(i+1) = 0;  % no middle node in this column
        end
        
        % Bottom node
        nodeCount = nodeCount + 1;
        botID(i+1) = nodeCount;
        nodes_xy(nodeCount, :) = [x_i, 0];
    end
    
    % --- 3) Build the element connectivity ---
    % We'll define:
    %   - Vertical cables in columns that have a middle node (top->mid->bot).
    %   - Horizontal cables (top chord, bottom chord) between columns i-1 and i.
    %   - Diagonals (bars) that form each diamond cell in each bay i=1..n.
    %     Typically: top_{i-1}->bottom_i, bottom_{i-1}->top_i, etc.
    %     We can also connect to the middle node if it exists on either side.
    
    el = [];
    addEl = @(type, n1, n2) [type, n1, n2];
    
    % 3a) Vertical cables in columns that have a middle node
    for i = 0:n
        if midID(i+1) ~= 0
            % top->mid
            el = [el; addEl(2, topID(i+1), midID(i+1))];
            % mid->bottom
            el = [el; addEl(2, midID(i+1), botID(i+1))];
        end
    end
    
    % 3b) Horizontal cables along top & bottom chords
    for i = 1:n
        el = [el;
              addEl(2, topID(i),   topID(i+1));  % top chord
              addEl(2, botID(i),   botID(i+1))]; % bottom chord
    end
    
    % 3c) Diagonals (bars) for each bay i=1..n
    %     Usually we do an "X" between columns (i-1) and i:
    %       top_{i-1} -> bottom_i
    %       bottom_{i-1} -> top_i
    %     plus connections to the middle node(s) if present.
    for i = 1:n
        iL = i - 1;   % left column index
        iR = i;       % right column index
        
        % Basic "X"
        el = [el;
              addEl(1, topID(iL+1), botID(iR+1));  % top_{i-1} -> bottom_i
              addEl(1, botID(iL+1), topID(iR+1))]; % bottom_{i-1} -> top_i
        
        % If left column has mid node, connect it to top_i and bottom_i
        if midID(iL+1) ~= 0
            el = [el;
                  addEl(1, midID(iL+1), topID(iR+1));
                  addEl(1, midID(iL+1), botID(iR+1))];
        end
        
        % If right column has mid node, connect it to top_{i-1} and bottom_{i-1}
        if midID(iR+1) ~= 0
            el = [el;
                  addEl(1, topID(iL+1), midID(iR+1));
                  addEl(1, botID(iL+1), midID(iR+1))];
        end
    end
    
    % --- 4) Boundary conditions: pinned at columns 0 and n ---
    % That is, fix x and y for top/mid/bottom at col0 and colN (if mid exists).
    n_fix = [];
    % pin col0
    n_fix = addPin(n_fix, topID(1));
    if midID(1)~=0, n_fix = addPin(n_fix, midID(1)); end
    n_fix = addPin(n_fix, botID(1));
    % pin colN
    n_fix = addPin(n_fix, topID(n+1));
    if midID(n+1)~=0, n_fix = addPin(n_fix, midID(n+1)); end
    n_fix = addPin(n_fix, botID(n+1));
end

function n_fix = addPin(n_fix, nodeID)
% addPin  appends x=0 and y=0 constraints for the given node
    n_fix = [n_fix;
             nodeID, 1, 0;  % fix x
             nodeID, 2, 0]; % fix y
end
