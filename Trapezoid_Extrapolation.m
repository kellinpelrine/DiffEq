clear;
options = optimoptions('fsolve','Display','none');
lambda = 6.7;
iter = 2;
eiter = 10; %max extrapolation iterations
sol = [];
for l = 1:1000 %check a grid with 1000 values of lambda between 6.7 and 6.8
    lambda = lambda + .0001;
for j = 1:eiter  
    t = 0;
    iterj = iter*(2^(j - 1));
    h = 1/iterj; %stepsize
    y = zeros(iter*2^eiter, 2);
    y(1, :) = [0 1]; %initial values
    for i = 2:iterj + 1
        x = fsolve(@(x)trap(x, y(i - 1, :), t, h, lambda), y(i - 1,:), options);
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
    sol(l) = ex(j,j);
end
end
[m, i] = min(abs(sol)) %display solution [value of y(1), lambda = 6.7 + .0001*i]

function [f] = fun(t, y, lambda) %the system to solve. Converted from the original equation using the standard substitutions u = y, v = y'
x1 = y(1,2);
x2 = y(1,2)/(1 + t) - (1 + t)*lambda*y(1,1);
f = [x1 x2];
end

function [tra] = trap(y1, y2, t, h, lambda) %trapezoid rule step
tra = y1 - y2 - 1/2 * h * (fun(t, y2, lambda) + fun(t + h, y1, lambda));
end
