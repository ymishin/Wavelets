function fvec = forward_transform_2d(xvec, yvec, fvec, jmax, jmin, porder)
% Forward interpolating wavelet transform is performed on 2D field given by
% values {fvec} at points {xvec},{yvec}. Number of given values of {fvec}
% should be Mx*2^(jmax-1)+1 and My*2^(jmax-1)+1 in x- and y-directions (Mx and
% My are arbitrary). Transform is performed from level {jmax} to level {jmin}.
% Array {porder} determines polynomial order of transform:
% [p_po - order for predict stage
%  u_po - order for update stage]

% loop over levels
for j = jmax:-1:(jmin+1)
    
    % step
    s = 2^(jmax-j);
    
    % transform from current level to lower one on X-slices
    for iy = 1:s:size(fvec,1)
        fvec(iy,:) = forward_transform_step(xvec, fvec(iy,:), s, porder);
    end
    
    % transform from current level to lower one on Y-slices
    for ix = 1:s:size(fvec,2)
        fvec(:,ix) = forward_transform_step(yvec, fvec(:,ix), s, porder);
    end
    
end

end