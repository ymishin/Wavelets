function fvec = compress(fvec, jmax, jmin, eps)
% Compress 1D field given by {fvec} which underwent forward transform from
% level {jmax} to level {jmin}. As the result of compression all d-coefficients
% at all levels below {eps} are zeroized.

if (eps < 0)
    return;
end

% loop over levels
for j = jmax:-1:(jmin+1)
    
    % step
    s = 2^(jmax-j);
    
    % compress at current level
    fvec = compress_step(fvec, s, eps);
    
end

end