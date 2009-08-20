function test_transform_adapt_2d

close all;
clear all;

Mx = 8;
My = 8;
jmax = 3;
jmin = 1;
porder = [1 1]; % linear
%porder = [3 3]; % cubic
%porder = [5 5]; % 5th order
eps = 1e-5;

% sample function
nx = Mx * 2^(jmax-1)+1;
ny = My * 2^(jmax-1)+1;
xvec = linspace(-0.5,0.5,nx);
yvec = linspace(-0.5,0.5,ny);
[X Y] = meshgrid(xvec,yvec);
A = 1.0; sigma = 0.07; dsigma = 0.05;
fmat = A * exp(-0.5*(X.^2/sigma^2 + Y.^2/sigma^2)); % gaussian
fmatp = fmat;
enorm = max(fmat(:)) - min(fmat(:));
eps = eps * enorm; % normalize
mask = ones(ny,nx); % initially all nodes are marked by mask

for k = 1:2
    
    % perform forward transform
    fmat1 = forward_transform_2d(xvec, yvec, fmat, jmax, jmin, porder, mask);
    
    % build adapted grid
    mask = adapt_grid_2d(fmat1, jmax, jmin, porder, eps);
    
    % plot
    figure;
    surf(xvec, yvec, fmatp); % sample function
    % loop over levels
    for j = jmax:-1:(jmin+1)
        s = 2^(jmax-j);
        figure;
        hold on;
        % plot all sample points
        plot(X(:), Y(:), 'o', 'MarkerEdgeColor', 'k', 'MarkerSize', 3);
        % plot d coefficients at current level
        for iy = 1:s:ny
            for ix = (s+1):2*s:nx
                plot(xvec(ix), yvec(iy), 'o', 'MarkerEdgeColor', 'r', ...
                    'MarkerSize', 5);
                if (mask(iy,ix))
                    % fill in if it's in the mask
                    plot(xvec(ix), yvec(iy), 'o', 'MarkerEdgeColor', 'r', ...
                        'MarkerSize', 5, 'MarkerFaceColor', 'r');
                end
            end
        end
        for ix = 1:s:nx
            for iy = (s+1):2*s:ny
                plot(xvec(ix), yvec(iy), 'o', 'MarkerEdgeColor', 'r', ...
                    'MarkerSize', 5);
                if (mask(iy,ix))
                    % fill in if it's in the mask
                    plot(xvec(ix), yvec(iy), 'o', 'MarkerEdgeColor', 'r', ...
                        'MarkerSize', 5, 'MarkerFaceColor', 'r');
                end
            end
        end
        hold off;
    end
    figure;
    hold on;
    % plot all sample points
    plot(X(:), Y(:), 'o', 'MarkerEdgeColor', 'k', 'MarkerSize', 3);
    % plot c coefficients at lowest level
    for iy = 1:2*s:ny
        for ix = 1:2*s:nx
            plot(xvec(ix), yvec(iy), 'o', 'MarkerEdgeColor', 'b', ...
                'MarkerSize', 5, 'MarkerFaceColor', 'b');
        end
    end
    hold off;
    
    % sample function for next time step at adapted grid
    sigma = sigma + dsigma;
    fmat = A * exp(-0.5*(X.^2/sigma^2 + Y.^2/sigma^2));
    fmatp = fmat;
    fmat(find(~mask)) = 0; % leave only values at adapted grid
    
end

end