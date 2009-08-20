function mask = adapt_grid(fvec, jmax, jmin, porder, eps)
% Grid adaptation is performed for 1D field given by values {fvec}, which was 
% previously transformed from level {jmax} to level {jmin}. The new adapted
% grid is returned via {mask} and contains: significant nodes (as prescribed 
% by {eps}), adjacent to significant nodes, nodes at all levels which are 
% necessary to compute these significant and adjacent nodes (reconstruction 
% check criterion), and all nodes at coarsest level (c-coefficients).
% Array {porder} determines polynomial order of transform:
% [p_po - order for predict stage
%  u_po - order for update stage]
%
% For details see: 
% O.V.Vasilyev and C.Bowman, "Second Generation Wavelet Collocation Method for 
% the Solution of Partial Differential Equations",J.Comp.Phys.,165,660-693,2000.

p_po = porder(1); % predict polynomial order

nx = length(fvec);

mask = zeros(nx,1);

% add significant nodes to the mask -
% d-coefficients at all levels with values above eps
for j = jmax:-1:(jmin+1)
    
    % step
    s = 2^(jmax-j);
    
    % loop over d-coefficients
    for i = (s+1):2*s:nx
        % add significant node to the mask
        mask(i) = abs(fvec(i)) > eps;
    end
    
end

mask_adj = mask;
% add to the mask nodes adjacent to significant nodes
for j = jmax:-1:(jmin+1)
    
    % step
    s = 2^(jmax-j);
    % step for adjacent nodes at current level
    s_adj1 = s;
    % step for adjacent nodes at one above level
    s_adj2 = 2^(jmax-min(j+1,jmax));
    
    % loop over d-coefficients
    for i = (s+1):2*s:nx
        % check if it's a significant node
        if (~mask(i))
            continue;
        end
        % add to the mask adjacent nodes at current level
        L = max(1,i-s_adj1);
        R = min(i+s_adj1,nx);
        mask_adj(L:s_adj1:R) = true;
        % add to the mask adjacent nodes at one above level
        L = max(1,i-s_adj2);
        R = min(i+s_adj2,nx);
        mask_adj(L:s_adj2:R) = true;
    end
    
end
mask = mask_adj;

% reconstruction check - add to the mask nodes at all
% levels which are necessary to compute significant and adjacent nodes
for j = jmax:-1:(jmin+1)
    
    % step
    s = 2^(jmax-j);
    
    % loop over d-coefficients
    for i = (s+1):2*s:nx
        % check if current node is in the mask
        if (~mask(i))
            continue;
        end
        % add to the mask nodes which are 
        % necessary to compute current d-coefficient
        L = max(1,i-p_po*s);
        R = min(L+2*p_po*s,nx);
        L = max(1,R-2*p_po*s);
        mask(L:2*s:R) = true;
    end
    
end

% add c-coefficients at coarsest level to the mask
for i = 1:2*s:nx
    mask(i) = true;
end

end