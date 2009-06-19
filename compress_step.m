function fvec = compress_step(fvec, s, eps)
% Compress 1D field given by {fvec} at one level. As the result of compression
% all d-coefficients at positions with step {s} below {eps} are zeroized.

nx = length(fvec);

% loop over d-coefficients
for i = (s+1):2*s:nx
    
    % zeroize
    if (abs(fvec(i)) < eps)
        fvec(i) = 0;
    end
    
end

end