clear;
load('Poisson_DST_square.mat') %load square (for boundary, when comparing with matlab solution)
close all;
problem = 1;
M = 200;
N = 200; %gridsizes
gridx = linspace(0,1,M);
gridy = linspace(0,1,N);
cgridx = linspace(1,M,M).^2;
cgridy = linspace(1,N,N).^2;
cgrid = cgridx + cgridy'; %for converting coefficients for f to ones for u

fgrid = f(gridx,gridy,1);
fhat = sintrans(sintrans(fgrid,1),2);

%confirming (with known values) that inverse is returning to the original function
%note that DST is its own inverse up to scaling
funhat = 4/(M*N) * sintrans(sintrans(fhat,2),1); 
invconfirm = sum(sum(abs(fgrid - funhat))) %0 => inverse = original;

uhat = -fhat ./ (pi^2 * cgrid);
ugrid = 8/(M * N) * sintrans(sintrans(uhat,2),1);
su = surf(ugrid);
view(2)
colorbar
colormap gray
set(su, 'edgecolor','none')
xlim([1 M]);
ylim([1 N]);
figure

%compare with matlab results
 model = createpde();
 geometryFromEdges(model, decsg(gd,sf,ns)); %pre-specified as square
 applyBoundaryCondition(model,'dirichlet','edge',1:model.Geometry.NumEdges,'u',0);
 specifyCoefficients(model,'m',0,...
                          'd',0,...
                          'c',-1,...
                          'a',0,...
                          'f',@f1);
generateMesh(model,'Hmax',0.02);
results = solvepde(model);
pdeplot(model,'XYData',results.NodalSolution)
colormap gray

%try another f
problem = 2;

fgrid = f(gridx,gridy,problem);
fhat = sintrans(sintrans(fgrid,1),2);
funhat = 4/(M*N) * sintrans(sintrans(fhat,2),1); 
invconfirm2 = sum(sum(abs(fgrid - funhat))) %0 => inverse = original;
uhat = -fhat ./ (pi^2 * cgrid);
ugrid = 16/(M * N) * sintrans(sintrans(uhat,1),2); %extra factor needed to match built in solver: 4
figure
su = surf(ugrid);
view(2)
colorbar
colormap gray
set(su, 'edgecolor','none')
xlim([1 M]);
ylim([1 N]);
figure

%compare with matlab results
model = createpde();
geometryFromEdges(model, decsg(gd,sf,ns));
applyBoundaryCondition(model,'dirichlet','edge',1:model.Geometry.NumEdges,'u',0);
specifyCoefficients(model,'m',0,...
                      'd',0,...
                      'c',-1,...
                      'a',0,...
                      'f',@f2);
generateMesh(model,'Hmax',0.02);
results = solvepde(model);
pdeplot(model,'XYData',results.NodalSolution)
colormap gray


function [s] = sintrans(X,dir) %dir = direction, down row or across column
    s1 = -imag(fft(X,size(X,dir)*2,dir));
    if dir == 1
        s = s1(1:size(X,1),:);
    else
        s = s1(:,1:size(X,2));
    end
end

function [fo] = f(X,Y,i) %i = problem index
    if i == 1
        fo = sin(2*pi*X) .* sin(2*pi*Y');
    end
    if i == 2
        fo = X.*(1 - X).*Y'.*(1 - Y');
    end
end

function [fo] = f1(region,state)
    fo = sin (2*pi*region.x).*sin(2*pi*region.y);    
end

function [fo] = f2(region,state)
   fo = region.x .* (1 - region.x) .* region.y .* (1 - region.y);
end