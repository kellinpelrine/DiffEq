clear;
endt = 1;
tstepcount = 400;
tstep = endt/tstepcount;
xstepcount = 20;
xstep = 1/xstepcount;

problem = 1; %problem to solve
sol = CN(tstepcount,tstep,xstepcount,xstep,problem);
solexp = explicit(tstepcount,tstep,xstepcount,xstep,problem);
exact = truegrid(tstepcount,tstep,xstepcount,xstep,problem);

errorCN = mean(mean((sol - exact).^2))
errorexp = mean(mean((solexp - exact).^2))

%Check stability (lack thereof for explicit)
xstepcount = 50;
xstep = 1/xstepcount;

sol = CN(tstepcount,tstep,xstepcount,xstep,problem);
solexp = explicit(tstepcount,tstep,xstepcount,xstep,problem);
exact = truegrid(tstepcount,tstep,xstepcount,xstep,problem);

errorCN = mean(mean((sol - exact).^2))
errorexp = mean(mean((solexp - exact).^2))

%next problem
problem = 2; 
xstepcount = 10;
xstep = 1/xstepcount;

sol = CN(tstepcount,tstep,xstepcount,xstep,problem);
solexp = explicit(tstepcount,tstep,xstepcount,xstep,problem);
exact = truegrid(tstepcount,tstep,xstepcount,xstep,problem);

errorCN = mean(mean((sol - exact).^2))
errorexp = mean(mean((solexp - exact).^2))

%stability check
xstepcount = 30;
xstep = 1/xstepcount;

sol = CN(tstepcount,tstep,xstepcount,xstep,problem);
solexp = explicit(tstepcount,tstep,xstepcount,xstep,problem);
exact = truegrid(tstepcount,tstep,xstepcount,xstep,problem);

errorCN = mean(mean((sol - exact).^2))
errorexp = mean(mean((solexp - exact).^2))


%reversing time
problem = 1; 
endt = -1;
tstepcount = 400;
tstep = endt/tstepcount;
xstepcount = 20;
xstep = 1/xstepcount;

sol = CN(tstepcount,tstep,xstepcount,xstep,problem);
reversetime = min(abs(sol(tstepcount,2:xstepcount))) %showing solution explodes at every x (except the fixed boundary).

%wave
problem = 3;
endt = 1;

tstepcount = 1000; %we need h_t < h_x for stability (so this is overkill from that perspective)
tstep = endt/tstepcount;
xstepcount = 30;
xstep = 1/xstepcount;

solwave = wave(tstepcount,tstep,xstepcount,xstep,problem);
exact = truegrid(tstepcount,tstep,xstepcount,xstep,problem);
errorwave = mean(mean((solwave(2:end) - exact(1:end - 1)).^2))

%ill-posed elliptic
problem = 4;

tstepcount = 1000;
tstep = endt/tstepcount;
xstepcount = 30;
xstep = 1/xstepcount;

solwave = wave(tstepcount,tstep,xstepcount,xstep,problem);
reverseRHSsign = min(abs(solwave(tstepcount,2:xstepcount))) %like above, reversing the RHS sign though instead of time. Solution explodes. 

%wave 2
problem = 5;
endt = 1;

tstepcount = 1000;
tstep = endt/tstepcount;
xstepcount = 30;
xstep = 1/xstepcount;

solwave = wave(tstepcount,tstep,xstepcount,xstep,problem);
exact = truegrid(tstepcount,tstep,xstepcount,xstep,problem);
errorwave = mean(mean((solwave(2:end) - exact(1:end - 1)).^2))

function [ao] = a(x,i) %i = problem index
if i == 1
    ao = 1/pi^2; 
end
if i == 2
    ao = 1;
end
if i == 3
   ao = 1; 
end
if i == 4
   ao = -1; 
end
if i == 5;
    ao = x*(1 - x);
end
end

function [adero] = ader(x,i) %derivative of a
if i == 1
    adero = 0;
end
if i == 2
    adero = 0;
end
if i == 3
    adero = 0;
end
if i == 4
    adero = 0;
end
if i == 5
    adero = 1 - 2*x;
end
end

function [fo] = f(x,t,i)
if i == 1
    fo = 0;
end
if i == 2
    fo = 0;
end
if i == 3
    fo = (2 - x*(1 - x))*sin(t);
end
if i == 4
    fo = -(2 - x*(1 - x))*sin(t);
end
if i == 5
    fo = -(x - 1)*x*(4*x - 3)*sin(t);
end
end

function [inito] = init(x,i) %initial condition
if i == 1
    inito = 1/(pi^2) * sin(pi*x);
end
if i == 2
    inito = 12 * sin(9*pi*x) - 7*sin(4*pi*x);
end
if i == 3
    inito = 0;
end
if i == 4
    inito = 0;
end
if i == 5
    inito = 0;
end
end

function [init2o] = init2(x,i) %derivative initial condition for wave equation
if i == 3
    init2o = x*(1 - x);
end
if i == 4
    init2o = x*(1 - x);
end
if i == 5
    init2o = x*(1 - x);
end
end

function [trueo] = true(x,t,i) %for computing error
if i == 1
    trueo = 1/(pi^2) .* exp(-t) .* sin(pi*x);
end
if i == 2
   trueo = 12*sin(9*pi*x)*exp( -(9*pi)^2 * t) - 7 * sin(4 * pi * x) * exp(-(4 * pi)^2 * t); 
end
if i == 3
   trueo = x*(1 - x) * sin(t);
end
if i == 5
   trueo = x*(1 - x) * sin(t);
end
end

function [exact] = truegrid(tstepcount,tstep,xstepcount,xstep,problem) %compute values of exact function
    for i = 1:tstepcount %computational efficiency of this (and probably some stuff below) could be improved by vectorization
        for j = 0:xstepcount
            exact(i,j + 1) = true(j*xstep,i*tstep,problem);
        end
    end
end

function [sol] = CN(tstepcount,tstep,xstepcount,xstep,problem) %Crank-Nicolson
    A = zeros(xstepcount - 1,xstepcount - 1);
    for i = 1:xstepcount - 1
        tempx = i*xstep;
        s = ader(tempx,problem)*tstep/(4*xstep);
        r = a(tempx,problem)*tstep/(2*(xstep^2));    

        b(i) = (r + s)*init(tempx + xstep,problem) + (1 - 2*r)*init(tempx,problem) ...
            + (r - s)*init(tempx - xstep,problem) + (tstep/2) * (f(tempx,tstep,problem) + f(tempx,0,problem)); 
        if i ~= 1 && i ~= xstepcount - 1
            A(i,i - 1) = s - r;
            A(i,i) = 1 + 2*r;
            A(i, i + 1) = -r - s;
        elseif i == 1
            A(1,1) = 1 + 2*r;
            A(1,2) = -r - s;
        else %i = xstepcount - 1
            A(xstepcount - 1,xstepcount - 1) = 1 + 2*r;
            A(xstepcount - 1,xstepcount - 2) = s - r;
        end
    end
    sol(1,:) = [0 (A\b')' 0]; %first time step
    for j = 2:tstepcount
        for i = 1:xstepcount - 1
            tempx = i*xstep;
            s = ader(tempx,problem)*tstep/(4*xstep);
            r = a(tempx,problem)*tstep/(2*xstep^2);    

            b(i) = (r + s)*sol(j - 1,i + 2) + (1 - 2*r)*sol(j - 1,i + 1) ...
                + (r - s)*sol(j - 1, i) + (tstep/2) * (f(tempx,j*tstep,problem) + f(tempx,(j - 1)*tstep,problem)); 
            if i ~= 1 && i ~= xstepcount - 1
                A(i,i - 1) = s - r;
                A(i,i) = 1 + 2*r;
                A(i, i + 1) = -r - s;
            elseif i == 1
                A(1,1) = 1 + 2*r;
                A(1,2) = -r - s;
            else %i = xstepcount - 1
                A(xstepcount - 1,xstepcount - 1) = 1 + 2*r;
                A(xstepcount - 1,xstepcount - 2) = s - r;
            end
        end
        sol(j,:) = [0 (A\b')' 0];
    end
end

function [solexp] = explicit(tstepcount,tstep,xstepcount,xstep,problem) %Explicit finite difference
    for i = 1:xstepcount - 1
        tempx = i*xstep;
        v = ader(tempx,problem)*tstep/(2*xstep);
        w = a(tempx,problem)*tstep/(xstep^2);   
        solexp(1,i) = (w + v)*init(tempx + xstep,problem) + (1 - 2*w)*init(tempx,problem) ...
            + (w - v)*init(tempx - xstep,problem) + tstep * f(tempx,0,problem);
    end
    solexp = [0 solexp 0];
    for j = 2:tstepcount
        for i = 1:xstepcount - 1
            tempx = i*xstep;
            v = ader(tempx,problem)*tstep/(2*xstep);
            w = a(tempx,problem)*tstep/(xstep^2);   
            solhold(j,i) = (w + v)*solexp(j - 1,i + 2) + (1 - 2*w)*solexp(j - 1,i + 1) ...
                    + (w - v)*solexp(j - 1, i) + tstep * f(tempx,(j - 1)*tstep,problem); 
        end
        solexp(j,:) = [0 solhold(j,:) 0];
    end
end

function[solwave] = wave(tstepcount,tstep,xstepcount,xstep,problem) %explicit wave equation
    for i = 1:xstepcount - 1 %compute t = 0
        tempx = i*xstep;
        solwave(1,i) = init(tempx,problem);
    end
    solwave = [0 solwave 0];
    for i = 1:xstepcount - 1 %compute next time step
        tempx = i*xstep;
        v = ader(tempx,problem)*tstep^2/(2*xstep);
        w = a(tempx,problem)*tstep^2/(xstep^2);   
        hold(1,i) = (1/2)*(   (w + v)*solwave(1,i + 2) + (2 - 2*w)*solwave(1, i + 1) ...
            + (w - v)*solwave(1, i) + tstep^2 * f(tempx,0,problem) + 2*tstep*init2(tempx,problem)    );
    end    
    solwave(2,:) = [0 hold 0];
    for j = 3:tstepcount %compute the rest
        for i = 1:xstepcount - 1
            tempx = i*xstep;
            v = ader(tempx,problem)*tstep^2/(2*xstep);
            w = a(tempx,problem)*tstep^2/(xstep^2);   
            hold(1,i) =  (w + v)*solwave(j - 1,i + 2) + (2 - 2*w)*solwave(j - 1,i + 1) ...
                + (w - v)*solwave(j - 1,i) + tstep^2 * f(tempx,(j - 1)*tstep,problem) - solwave(j - 2,i + 1);
        end 
         solwave(j,:) = [0 hold 0];
    end
end