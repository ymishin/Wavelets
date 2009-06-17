function test_transform_2d

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
A = 1.0; sigmax = 0.07; sigmay = 0.07;
fmat = A * exp(-0.5*(X.^2/sigmax^2 + Y.^2/sigmay^2)); % gaussian

% perform forward transform
fmat1 = forward_transform_2d(xvec, yvec, fmat, jmax, jmin, porder);

% get rid of values below eps
fmat1(find(abs(fmat1) < eps)) = 0;

% perform inverse transform
fmat2 = inverse_transform_2d(xvec, yvec, fmat1, jmax, jmin, porder);

% compute max error
err = max(max(abs(fmat - fmat2)));
disp(err);

% plot
figure;
surf(xvec, yvec, fmat); % sample function
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
            if (abs(fmat1(iy,ix)) > 0)
                % fill in if value > 0
                plot(xvec(ix), yvec(iy), 'o', 'MarkerEdgeColor', 'r', ...
                    'MarkerSize', 5, 'MarkerFaceColor', 'r');
            end
        end
    end
    for ix = 1:s:nx
        for iy = (s+1):2*s:ny
            plot(xvec(ix), yvec(iy), 'o', 'MarkerEdgeColor', 'r', ...
                'MarkerSize', 5);
            if (abs(fmat1(iy,ix)) > 0)
                % fill in if value > 0
                plot(xvec(ix), yvec(iy), 'o', 'MarkerEdgeColor', 'r', ...
                    'MarkerSize', 5, 'MarkerFaceColor', 'r');
            end
        end
    end
    hold off;
end
% plot c coefficients at lowest level
figure;
hold on;
% plot all sample points
plot(X(:), Y(:), 'o', 'MarkerEdgeColor', 'k', 'MarkerSize', 3);
for iy = 1:2*s:ny
    for ix = 1:2*s:nx
        plot(xvec(ix), yvec(iy), 'o', 'MarkerEdgeColor', 'b', 'MarkerSize', 5);
        if (abs(fmat1(iy,ix)) > 0)
            % fill in if value > 0
            plot(xvec(ix), yvec(iy), 'o', 'MarkerEdgeColor', 'b', ...
                'MarkerSize', 5, 'MarkerFaceColor', 'b');
        end
    end
end
hold off;

end