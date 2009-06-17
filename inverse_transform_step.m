function fvec = inverse_transform_step(xvec, fvec, s, porder)
% Inverse interpolating wavelet transform is performed on 1D field given
% by values {fvec} at points {xvec}. Only values at positions with step {s}
% are considered. Array {porder} determines polynomial order of transform.

p_po = porder(1); % predict polynomial order
u_po = porder(2); % update polynomial order

nx = length(xvec);

% UN-UPDATE
% loop over values of function which were updated
for i = 1:2*s:nx
    
    % left and right boundaries of the interval to use for interpolation
    L = max(s+1,i-u_po*s);
    R = min(L+2*u_po*s,nx-s);
    L = max(s+1,R-2*u_po*s);
    
    % interpolate value at point i
    f = lagrange_interp(xvec(i), xvec(L:2*s:R), fvec(L:2*s:R));
    
    % un-update value of function
    fvec(i) = fvec(i) - f;
    
end

% UN-PREDICT
% loop over values of function which were predicted
for i = (s+1):2*s:nx
    
    % left and right boundaries of the interval to use for interpolation
    L = max(1,i-p_po*s);
    R = min(L+2*p_po*s,nx);
    L = max(1,R-2*p_po*s);
    
    % interpolate value at point i
    f = lagrange_interp(xvec(i), xvec(L:2*s:R), fvec(L:2*s:R));
    
    % un-predict value of function
    fvec(i) = 2.0 * fvec(i) + f;
    
end

end