%% 1D Multi-layer Additive Manufacturing with locally refined Mesh
% evaluate a 1D poisson transient problem using h-FEM

clear all;
clc;

writerObj = VideoWriter('TravellingSource.avi');
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

timeSteps = 20;
numberOfElementsInX = 25;
refinementDepth = 5;

t = linspace(0, tEnd, timeSteps + 1);                                       % time discretization
x = linspace(0.0, xEnd, numberOfElementsInX + 1);                           % spatial discretization X
x_ref = linspace(0.0, 1.0, 2^refinementDepth + 1);                          % refinement discretization X_ref

[X, T] = meshgrid(x, t);
[X_ref, T_ref] = meshgrid(x_ref, t);


%% Analysis and plot

[temperatureSolution, heatFlux, temperatureSolutionRefined, heatFluxRefined] = backwardEulerRefined(x, rhs, dirichletLeftBC, dirichletRightBC, k, heatCapacity, t, refinementDepth);

figure(2)
surf(X, T, temperatureSolution')

figure(3)
surf(X, T, heatFlux')

figure(4)
surf(X_ref, T_ref, temperatureSolutionRefined')

figure(5)
surf(X_ref, T_ref, heatFluxRefined')

% figure(4)
% F(size(t,2)) = struct('cdata',[],'colormap',[]);
% 
% for i=1:size(t,2)
%     plot(x',temperatureSolution(:,i))
%     drawnow
%     F(i) = getframe;
%     writeVideo(writerObj, getframe(gcf, [ 0 0 560 420 ]));
% end

close(writerObj);