function fvec = compress_2d(fvec, jmax, jmin, eps)
% Compress 2D field given by {fvec} which underwent forward transform from
% level {jmax} to level {jmin}. As the result of compression all d-coefficients
% at all levels below {eps} are zeroized.
%
% $Id$

if (eps < 0)
    return;
end

% loop over levels
for j = jmax:-1:(jmin+1)
    
    % step
    s = 2^(jmax-j);
    
    % compress at current level on X-slices
    for iy = 1:s:size(fvec,1)
        fvec(iy,:) = compress_step(fvec(iy,:), s, eps);
    end
    
    % compress at current level on Y-slices
    for ix = 1:s:size(fvec,2)
        fvec(:,ix) = compress_step(fvec(:,ix), s, eps);
    end
    
end

end