clear;
options = optimoptions('fsolve','Display','none');

iter = 2;
eiter = 10; %max extrapolation iterations
for j = 1:eiter  
    t = 0.00001;
    iterj = iter*(2^(j - 1));
    h = (3*pi - t)/iterj; %stepsize
    y = zeros(iter*2^eiter, 2);
    y(1, :) = [.000005 .5]; %initial values
    for i = 2:iterj + 1
        x = fsolve(@(x)trap(x, y(i - 1, :), t, h), y(i - 1,:), options);
        y(i,:) = x;
        t = t + h;
    end
    z(j) = x(1); %proposed solution before extrapolating
    ex = [];
    ex(:,1) = z'; %extrapolate. Note that this could be improved slightly efficiency-wise (some things being recomputed unnecessarily) 
    if j > 1 
        for m = 2:j
            for k = 2:m
                ex(m,k) = ex(m, k - 1) + (ex(m, k - 1) - ex(m - 1, k - 1))/(2^k - 1);
            end
        end
        if abs(ex(j,j) - ex(j, j - 1)) < 10^-7 %stop iteration if change is sufficiently small
            break
        end
    end
end

function [f] = fun(t, y) %the system to solve. Converted from the original equation using the standard substitutions u = y, v = y'
x1 = y(1,2);
x2 = -(t * y(1,2) + (t^2 - 1)*y(1,1))/(t^2);
f = [x1 x2];
end

function [tra] = trap(y1, y2, t, h) %trapezoid rule step
tra = y1 - y2 - 1/2 * h * (fun(t, y2) + fun(t + h, y1));
end
