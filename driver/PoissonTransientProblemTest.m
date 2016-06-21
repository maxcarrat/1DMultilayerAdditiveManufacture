%% Poisson Transient Problem Test 1D
% evaluate a 1D poisson transient problem using h-FEM

clear all;
clc;

writerObj = VideoWriter('PoissonTest.avi');
writerObj.Quality = 100;
writerObj.FrameRate = 5;
open(writerObj);

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 7200.0 ;                                      % density [kg/m^3]
c = 440.5;                                          % specific heat [J/(kg°C)]
k = 35.0;                                           % thermal conductivity [W/(m°C)]
T0 = 0.0;                                           %Initial temperature [°C]
heatCapacity= rho*c;                                % heat capacity [kJ / kg °C]

tEnd = 32.0;
xEnd = 0.1e-00;

dirichletLeftBC = @(t) 0.0;
dirichletRightBC = @(t) T0 + 100.0 * sin( pi * t / 40 );
rhs = @(x, t) 10.0e+06 * cos( pi * x / 0.05 );

timeSteps = 32;
numberOfElementsInX = 40;

t = linspace(0, tEnd, timeSteps + 1);                                       % time discretization
x = linspace(0.0, xEnd, numberOfElementsInX + 1);                           % spatial discretization X
xPosProcessing = linspace(0.0, xEnd, 100);                                  % spatial discretization X post-processing

[X, T] = meshgrid(xPosProcessing, t);

%% Analysis and plot
[temperatureSolution, heatFlux, internalEnergy] = backwardEulerSolver(x, xPosProcessing, rhs, dirichletLeftBC, dirichletRightBC, k, heatCapacity, t);

figure(4)
surf(X, T, temperatureSolution')

figure(5)
F(size(t,2)) = struct('cdata',[],'colormap',[]);
for i=1:size(t,2)
    plot(xPosProcessing',temperatureSolution(:,i))
    drawnow
    F(i) = getframe;
    writeVideo(writerObj, getframe(gcf, [ 0 0 560 420 ]));
end
close(writerObj);