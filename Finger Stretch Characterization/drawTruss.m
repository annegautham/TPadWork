function [] = drawTruss(nodes_xy, el, Fint)
% Draw the truss structure and colorcode elements according to internal
% force, if provided
% 
% Format of input arguments are identical to those in solveTruss.m

[n_el, ~] = size(el);        % Number of elements
[n_n, ~] = size(nodes_xy);   % Number of nodes

% line spec and width for each element type
lspec = {'ko-', 'r-'};
lwidth = [3, 1];

% Setup options for colorcoding
if nargin < 3
    colorcoding = 0;
else
    colorcoding = 1;

    % Create color scale
    Fscale = max(Fint) - min(Fint);
    cmap = colormap('jet');
    n_colors = length(cmap);
end


figure(1)
hold on

% Draw each element as a solid line
for e = 1:n_el
    
    n1 = el(e, 2);
    n2 = el(e, 3);
    el_type = el(e, 1);
    
    if el_type > length(lspec)
        error('Unrecognized element type: element number: %d', e);
    else
        if colorcoding
            % Calculate color based on internal force
            f_norm = (Fint(e) - min(Fint)) / Fscale;
            c = cmap(round(f_norm * (n_colors-1)) + 1, :);
    
            plot(nodes_xy([n1, n2], 1), nodes_xy([n1, n2], 2),...
            lspec{el_type},'linewidth', lwidth(el_type),...
            'Color', c);
        else
            plot(nodes_xy([n1, n2], 1), nodes_xy([n1, n2], 2),...
            lspec{el_type},'linewidth', lwidth(el_type));
        end
    end
    
end

% Label nodes
for n = 1:n_n
    text(nodes_xy(n, 1), nodes_xy(n, 2), ['    ' num2str(n)],...
         'HorizontalAlignment', 'left')
end

% Create a color legend
colorbar('Ticks', [0, 0.5, 1], 'TickLabels', {'min', '', 'max'})