function mask = adapt_grid_2d(fmat, jmax, jmin, porder, eps)
% Grid adaptation is performed for 2D field given by values {fmat}, which was 
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
%
% $Id$

p_po = porder(1); % predict polynomial order

ny = size(fmat,1);
nx = size(fmat,2);

mask = zeros(size(fmat));

% add significant nodes to the mask -
% d-coefficients at all levels with values above eps
for j = jmax:-1:(jmin+1)
    
    % step
    s = 2^(jmax-j);
    
    % loop over d-coefficients on X-slices
    for iy = 1:s:ny
        for ix = (s+1):2*s:nx
            % add significant node to the mask
            mask(iy,ix) = abs(fmat(iy,ix)) > eps;
        end
    end
    
    % loop over d-coefficients on Y-slices
    for ix = 1:s:nx
        for iy = (s+1):2*s:ny
            % add significant node to the mask
            mask(iy,ix) = abs(fmat(iy,ix)) > eps;
        end
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
    
    % loop over d-coefficients on X-slices
    for iy = 1:s:ny
        for ix = (s+1):2*s:nx
            % check if it's a significant node
            if (~mask(iy,ix))
                continue;
            end
            % add to the mask adjacent nodes at current level
            Lx = max(1,ix-s_adj1);
            Rx = min(ix+s_adj1,nx);
            Ly = max(1,iy-s_adj1);
            Ry = min(iy+s_adj1,ny);
            mask_adj(Ly:s_adj1:Ry,Lx:s_adj1:Rx) = true;
            % add to the mask adjacent nodes at one above level
            Lx = max(1,ix-s_adj2);
            Rx = min(ix+s_adj2,nx);
            Ly = max(1,iy-s_adj2);
            Ry = min(iy+s_adj2,ny);
            mask_adj(Ly:s_adj2:Ry,Lx:s_adj2:Rx) = true;
        end
    end
    
    % loop over d-coefficients on Y-slices
    for ix = 1:s:nx
        for iy = (s+1):2*s:ny
            % check if it's a significant node
            if (~mask(iy,ix))
                continue;
            end
            % add to the mask adjacent nodes at current level
            Lx = max(1,ix-s_adj1);
            Rx = min(ix+s_adj1,nx);
            Ly = max(1,iy-s_adj1);
            Ry = min(iy+s_adj1,ny);
            mask_adj(Ly:s_adj1:Ry,Lx:s_adj1:Rx) = true;
            % add to the mask adjacent nodes at one above level
            Lx = max(1,ix-s_adj2);
            Rx = min(ix+s_adj2,nx);
            Ly = max(1,iy-s_adj2);
            Ry = min(iy+s_adj2,ny);
            mask_adj(Ly:s_adj2:Ry,Lx:s_adj2:Rx) = true;
        end
    end
    
end
mask = mask_adj;

% reconstruction check - add to the mask nodes at all
% levels which are necessary to compute significant and adjacent nodes
for j = jmax:-1:(jmin+1)
    
    % step
    s = 2^(jmax-j);
    
    % loop over d-coefficients on X-slices
    for iy = 1:s:ny
        for ix = (s+1):2*s:nx
            % check if current node is in the mask
            if (~mask(iy,ix))
                continue;
            end
            % add to the mask nodes which are
            % necessary to compute current d-coefficient
            L = max(1,ix-p_po*s);
            R = min(L+2*p_po*s,nx);
            L = max(1,R-2*p_po*s);
            mask(iy,L:2*s:R) = true;
        end
    end
    
    % loop over d-coefficients on Y-slices
    for ix = 1:s:nx
        for iy = (s+1):2*s:ny
            % check if current node is in the mask
            if (~mask(iy,ix))
                continue;
            end
            % add to the mask nodes which are
            % necessary to compute current d-coefficient
            L = max(1,iy-p_po*s);
            R = min(L+2*p_po*s,ny);
            L = max(1,R-2*p_po*s);
            mask(L:2*s:R,ix) = true;
        end
    end
    
end

% add c-coefficients at coarsest level to the mask
for iy = 1:2*s:ny
    for ix = 1:2*s:nx
        mask(iy,ix) = true;
    end
end

end