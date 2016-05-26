%% 1D Multi-layer Additive Manufacturing with locally refined Mesh
% evaluate a 1D poisson transient problem using h-FEM

clear all;
clc;

writerObj = VideoWriter('AdditiveManufacturingRefined.avi');
writerObj.Quality = 100;
writerObj.FrameRate = 5;
open(writerObj);

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 1000.0;                                          % density [kg/m^3]
c = 10.0;                                              % specific heat [J/(kg°C)]
k = 10.0;                                              % thermal conductivity [W/(m°C)]
T0 = 20.0;                                             % Initial temperature [°C]
heatCapacity= rho*c;                                   % heat capacity [kJ / kg °C]
Tsource = 2000.0;                                      % source temperature [°C]


tEnd = 500.0;
xEnd = 1.0;

dirichletLeftBC = @(t) T0;
dirichletRightBC = @(t) T0 + Tsource;
rhs = @(x, t) 0.0;

timeSteps = 20;
numberOfElementsInX = 20;
refinementDepth = 3;

t = linspace(0, tEnd, timeSteps + 1);                                       % time discretization
x = linspace(0.0, xEnd, numberOfElementsInX + 1);                           % spatial discretization X
x_ref = linspace(0.0, 1.0, 2^refinementDepth + 1);                          % refinement discretization X_ref
x_PostProcess = linspace(0.0, xEnd, 3*(numberOfElementsInX + 1));           % post-processing coordinates             


[X, T] = meshgrid(x_PostProcess, t);
[X_ref, T_ref] = meshgrid(x_ref, t);


%% Analysis and plot

[temperatureSolution, heatFlux, temperatureSolutionRefined, heatFluxRefined, internalEnergy]...
    = backwardEulerRefined(x, x_PostProcess, rhs, dirichletLeftBC, dirichletRightBC, k, heatCapacity, t, refinementDepth);

figure(202)
surf(X, T, temperatureSolution')

figure(303)
surf(X, T, heatFlux')

figure(404)
F(size(t,2)) = struct('cdata',[],'colormap',[]);

for i=1:size(t,2)
    plot(x_PostProcess',temperatureSolution(:,i))
    drawnow
    F(i) = getframe;
    writeVideo(writerObj, getframe(gcf, [ 0 0 560 420 ]));
end

close(writerObj);