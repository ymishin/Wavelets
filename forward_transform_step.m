function fvec = forward_transform_step(xvec, fvec, s, porder)
% Forward interpolating wavelet transform is performed on 1D field given
% by values {fvec} at points {xvec}. Only values at positions with step {s}
% are considered. Array {porder} determines polynomial order of transform.

p_po = porder(1); % predict polynomial order
u_po = porder(2); % update polynomial order

nx = length(xvec);

% PREDICT
% loop over values of function which should be predicted
for i = (s+1):2*s:nx
    
    % left and right boundaries of the interval to use for interpolation
    L = max(1,i-p_po*s);
    R = min(L+2*p_po*s,nx);
    L = max(1,R-2*p_po*s);
    
    % interpolate value at point i
    f = lagrange_interp(xvec(i), xvec(L:2*s:R), fvec(L:2*s:R));
    
    % compute and store difference -> d coefficient
    fvec(i) = 0.5 * (fvec(i) - f);
    
end

% UPDATE
% loop over values of function which should be updated
for i = 1:2*s:nx
    
    % left and right boundaries of the interval to use for interpolation
    L = max(s+1,i-u_po*s);
    R = min(L+2*u_po*s,nx-s);
    L = max(s+1,R-2*u_po*s);
    
    % interpolate value at point i
    f = lagrange_interp(xvec(i), xvec(L:2*s:R), fvec(L:2*s:R));
    
    % update (lift) value of function -> c coefficient
    fvec(i) = fvec(i) + f;
    
end

end