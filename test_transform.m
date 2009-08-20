function test_transform

close all;
clear all;

M = 8;
jmax = 6;
jmin = 1;
%porder = [1 1]; % linear
porder = [3 3]; % cubic
%porder = [5 5]; % 5th order
eps = 1e-3;

% sample function
nx = M*2^(jmax-1)+1;
xvec = linspace(-0.5,0.5,nx);
fvec = cos(80*pi*xvec).*exp(-64*xvec.^2);
enorm = max(fvec) - min(fvec);
eps = eps * enorm; % normalize

% perform forward transform
fvec1 = forward_transform(xvec, fvec, jmax, jmin, porder, -1);

% compress - get rid of d-coefficients below eps
fvec1 = compress(fvec1, jmax, jmin, eps);

% perform inverse transform
fvec2 = inverse_transform(xvec, fvec1, jmax, jmin, porder);

% max error and compression ratio
err = max(abs(fvec - fvec2));
err = err / enorm; % normalize
disp(err);
comp_ratio = 100 * (1.0 - nnz(fvec1) / nnz(fvec));
disp(comp_ratio);

% plot
figure;
plot(xvec,fvec); % sample function
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
        if (abs(fvec1(i)) > 0)
            % fill in if value > 0
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

end