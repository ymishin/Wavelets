function fvec = forward_transform(xvec, fvec, jmax, jmin, porder, mask)
% Forward interpolating wavelet transform is performed on 1D field given
% by values {fvec} at points {xvec}. Only points marked by mask {mask} are 
% considered, if {mask} is equal -1 then all points are considered. Number of 
% given values of {fvec} should be M*2^(jmax-1)+1 (M is arbitrary). Transform 
% is performed from level {jmax} to level {jmin}. Array {porder} determines 
% polynomial order of transform:
% [p_po - order for predict stage
%  u_po - order for update stage]
%
% $Id$

% check mask
if (mask == -1)
    mask = ones(length(fvec),1); % mask will include all points
end

% loop over levels
for j = jmax:-1:(jmin+1)
    
    % step
    s = 2^(jmax-j);
    
    % transform from current level to lower one
    fvec = forward_transform_step(xvec, fvec, s, porder, mask);
    
end

end