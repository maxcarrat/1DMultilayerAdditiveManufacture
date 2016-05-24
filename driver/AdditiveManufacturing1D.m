%% 1D Multi-layer Additive Manufacturing
% evaluate a 1D poisson transient problem using h-FEM

clear all;
clc;

writerObj = VideoWriter('AdditiveManufacturing1D.avi');
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

tEnd = 500.0;
xEnd = 1.0;

dirichletLeftBC = @(t) T0;
dirichletRightBC = @(t) T0 + 200.0;
rhs = @(x, t) 0.0;

timeSteps = 50;
numberOfElementsInX = 50;

t = linspace(0, tEnd, timeSteps + 1);                                       % time discretization
x = linspace(0.0, xEnd, numberOfElementsInX + 1);                           % spatial discretization X
x_postProcessing = linspace(0.0, xEnd, 2*numberOfElementsInX + 1);            % spatial discretization X for Post-Processing

[X, T] = meshgrid(x_postProcessing, t);

%% Analysis and plot
[temperatureSolution, heatFlux, internalEnergy] = backwardEulerMultilayers(x, x_postProcessing, rhs, dirichletLeftBC, dirichletRightBC, k, heatCapacity, t);

figure(20)
surf(X, T, temperatureSolution')

figure(30)
surf(X, T, heatFlux')

figure(40)
F(size(t,2)) = struct('cdata',[],'colormap',[]);

for i=1:size(t,2)
    plot(x_postProcessing',temperatureSolution(:,i))
    drawnow
    F(i) = getframe;
    writeVideo(writerObj, getframe(gcf, [ 0 0 560 420 ]));
end

close(writerObj);