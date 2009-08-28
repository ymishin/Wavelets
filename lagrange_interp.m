function f = lagrange_interp(x, xvec, fvec)
% Perform Lagrange polynomial interpolation - interpolating polynomial
% is built to have values {fvec} at points {xvec}, and value of this
% polynomial at point {x} is returned. Polynomial order is (length(xvec)-1).
% see e.g. E.Kreyszig "Advanced Engineering Mathematics, 9th Ed."
%
% $Id$

nx = length(xvec);
f = 0;

% loop over given points
for k = 1:nx
    
    % construct and evaluate polynomial Lk
    % Lk has value 1 at point k, and 0 at all others given points
    Lk = (prod(x - xvec(1:k-1)) * prod(x - xvec(k+1:nx))) / ...
        (prod(xvec(k) - xvec(1:k-1)) * prod(xvec(k) - xvec(k+1:nx)));
    
    % update value of interpolating polynomial at point x
    f = f + fvec(k) * Lk;
    
end

end