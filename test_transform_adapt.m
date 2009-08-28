function test_transform_adapt
%
% $Id$

close all;
clear all;

M = 8;
jmax = 6;
jmin = 1;
%porder = [1 1]; % linear
porder = [3 3]; % cubic
%porder = [5 5]; % 5th order
eps = 1e-5;

% sample function
nx = M*2^(jmax-1)+1;
xvec = linspace(-0.5,0.5,nx);
x0 = -0.3; v = 1; t = 0; dt = 0.05; nu = 1e-2;
fvec = -tanh((xvec - x0 - v * t) ./ (2 * nu));
fvecp = fvec;
enorm = max(fvec) - min(fvec);
eps = eps * enorm; % normalize
mask = ones(nx,1); % initially all nodes are marked by mask

for k = 1:3
    
    % perform forward transform
    fvec1 = forward_transform(xvec, fvec, jmax, jmin, porder, mask);
    
    % build adapted grid
    mask = adapt_grid(fvec1, jmax, jmin, porder, eps);
    
    % plot
    figure;
    plot(xvec,fvecp); % sample function
    figure;
    axis([xvec(1) xvec(end) jmin-0.1 jmax+0.1]);
    grid on;
    hold on;
    % loop over levels
    for j = jmax:-1:(jmin+1)
        s = 2^(jmax-j);
        % plot d coefficients at current level
        for i = (s+1):2*s:nx
            plot(xvec(i), j, 'o', 'MarkerEdgeColor', 'r', 'MarkerSize', 5);
            if (mask(i))
                % fill in if it's in the mask
                plot(xvec(i), j, 'o', 'MarkerEdgeColor', 'r', 'MarkerSize', 5, ...
                    'MarkerFaceColor', 'r');
            end
        end
    end
    % plot c coefficients at lowest level
    for i = 1:2*s:nx
        plot(xvec(i), jmin, 'o', 'MarkerEdgeColor', 'b', 'MarkerSize', 5, ...
            'MarkerFaceColor', 'b');
    end
    hold off;
    
    % sample function for next time step at adapted grid
    t = t + dt;
    fvec = -tanh((xvec - x0 - v * t) ./ (2 * nu));
    fvecp = fvec;
    fvec(find(~mask)) = 0; % leave only values at adapted grid
    
end

end