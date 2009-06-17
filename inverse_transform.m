function fvec = inverse_transform(xvec, fvec, jmax, jmin, porder)
% Inverse interpolating wavelet transform is performed on 1D field given
% by values {fvec} at points {xvec}. Number of given values of {fvec} should
% be M*2^(jmax-1)+1 (M is arbitrary). Transform is performed from level {jmin}
% to level {jmax}. Array {porder} determines polynomial order of transform:
% [p_po - order for predict stage
%  u_po - order for update stage]

% loop over levels
for j = (jmin+1):1:jmax
    
    % step
    s = 2^(jmax-j);
    
    % transform from current level to higher one
    fvec = inverse_transform_step(xvec, fvec, s, porder);
    
end

end