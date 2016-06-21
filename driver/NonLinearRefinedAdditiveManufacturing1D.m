%% Non-linear Refined Additive Manufacturing 1D
% evaluate a 1D poisson non-linear transient problem using h-FEM and a
% refinement strategy toward the right-end-side of the bar

clear all;
clc;

writerObj = VideoWriter('NonLinearAdditiveManufacturing1D.avi');
writerObj.Quality = 100;
writerObj.FrameRate = 5;
open(writerObj);

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 7580.0;                                                       % density [kg/m^3]
c = 440.0;                                                          % specific heat [J/(kg°C)]
k =@(T) 54.0 - 26.7 * tanh( ( T - 800.0 ) / 800.0 + 1.0 );          % thermal conductivity [W/(m°C)]
T0 = 120.0;                                                         % Initial temperature [°C]
heatCapacity= rho*c;                                                % heat capacity [kJ / kg °C]
Tsource = 3000.0;                                                   % source temperature [°C]

tEnd = 500;
xEnd = 0.1;

dirichletLeftBC = @(t) T0;
dirichletRightBC = @(t) T0 + Tsource;
rhs = @(x, t) 0.0;

timeSteps = 10;
numberOfElementsInX = 10;
refinementDepth = 2;

tolerance = 1.0e-05;
maxNumberOfIterations = 15;

t = linspace(0, tEnd, timeSteps + 1);                                       % time discretization
x = linspace(0.0, xEnd, numberOfElementsInX + 1);                           % spatial discretization X
x_postProcessing = linspace(0.0, xEnd, 5*numberOfElementsInX*2^refinementDepth + 1);          % spatial discretization X for Post-Processing

[X, T] = meshgrid(x_postProcessing, t);

%% Analysis and plot
[temperatureSolution, heatFlux] = nonLinearBackwardEulerRefined(x, x_postProcessing, rhs, T0, dirichletLeftBC,...
    dirichletRightBC, k, heatCapacity, t, refinementDepth, maxNumberOfIterations, tolerance);

figure(2001)
surf(X, T, temperatureSolution')

figure(3001)
surf(X, T, heatFlux')

figure(4001)
F(size(t,2)) = struct('cdata',[],'colormap',[]);

for i=1:size(t,2)
    plot(x_postProcessing',temperatureSolution(:,i))
    drawnow
    F(i) = getframe;
    writeVideo(writerObj, getframe(gcf, [ 0 0 560 420 ]));
end

close(writerObj);