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

% perform forward transform
fvec1 = forward_transform(xvec, fvec, jmax, jmin, porder);

% get rid of values below eps
fvec1(find(abs(fvec1) < eps)) = 0;

% perform inverse transform
fvec2 = inverse_transform(xvec, fvec1, jmax, jmin, porder);

% compute max error
err = max(abs(fvec - fvec2));
disp(err);

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
    plot(xvec(i), jmin, 'o', 'MarkerEdgeColor', 'b', 'MarkerSize', 5);
    if (abs(fvec1(i)) > 0)
        % fill in if value > 0
        plot(xvec(i), jmin, 'o', 'MarkerEdgeColor', 'b', 'MarkerSize', 5, ...
            'MarkerFaceColor', 'b');
    end
end
hold off;

end