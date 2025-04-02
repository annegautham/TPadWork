function [nodes_xy, el, n_fix, P, EA_el] = genBridge2(n, L)
    nodes_xy = [];
    %generate nodes
    for i = 0:n+1
        if i == 0
            f1 = [0, 0];
            f2 = [0, L/2];
            f3 = [0, L];
            nodes_xy = [nodes_xy; f1; f2; f3];
    
        elseif i == n+1
            f4 = [(i-1)*L, 0];
            f5 = [(i-1)*L, L/2];
            f6 = [(i-1)*L, L];
            nodes_xy = [nodes_xy; f4; f5; f6];
        else
            ni1 = [(i-1)*L + L/2, 0];
            ni2 = [(i-1)*L + L/2, L];
            ni3 = [i*L, L/2];
            nodes_xy = [nodes_xy; ni1; ni2; ni3];
        end
    end
    nodes_xy = unique(nodes_xy, 'rows');
    
    disp(nodes_xy)
    scatter(nodes_xy(:,1), nodes_xy(:,2), 100, 'filled'); % Plot nodes
    axis equal
    a = 6+(n-1)*3;
    %generate fixed constraints
    n_fix = [1, 1, 0;
             1, 2, 0;
             2, 1, 0;
             2, 2, 0;
             3, 1, 0;
             3, 2, 0;
             a, 1, 0;
             a, 2, 0;
             a+1, 1, 0;
             a+1, 2, 0;
             a+2, 1, 0;
             a+2, 2, 0];
    
    %add cables (vertical + horizontal)
    el_vert =[];
    for i = 1:n
        el_vert = [el_vert; 2, 1+3*i, 2+3*i];
    end
    el_horz =[2, 1, 4;
              2, 3, 5;
              ];
    
    for i = 1:n-1
        q = 1+3*i;
        r = q+3;
        s = 2+3*i;
        t = s+3;
        el_horz = [el_horz; 2,q,r];
        el_horz = [el_horz; 2,s,t];
    end
    q = 1+3*n;
    
    el_horz = [el_horz;2, q, q+2];
    el_horz = [el_horz;2, q+1, q+4];
    
    el = [el_vert;el_horz];
    a = a+2;

    %add struts
    el_strut = [];
    el_strut = [1, 2, 4;
                1, 2, 5;
                1, a-1, a-3;
                1, a-1, a-4];
    
    for i = 1:n-1
        el_strut = [el_strut; 1, 6+(i-1)*3,6+(i-1)*3+1];
        el_strut = [el_strut; 1, 6+(i-1)*3,6+(i-1)*3+2];
    end
    
    for i = 1:n-1
        el_strut = [el_strut; 1, 6+(i-1)*3,6+(i-1)*3-1];
        el_strut = [el_strut; 1, 6+(i-1)*3,6+(i-1)*3-2];
    end
    
    el = [el; el_strut];
    P = [1.5*n+2.5, 2, -1];
    EA_el = [10, 1];

end