%% 1D Multi-layer Linear Additive Manufacturing with enriched local functions
% evaluate a 1D poisson transient linear problem using X-FEM

clear all;
clc;

writerObj = VideoWriter('AdditiveManufacturingPODBasis.avi');
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
xEnd = 0.10;

dirichletLeftBC = @(t) T0;
dirichletRightBC = @(t) T0 + Tsource;
rhs = @(x, t) 0.0;

timeSteps = 20;
maxTrainingTimeSteps = timeSteps*0.5;
relErrorEnergyNorm = zeros(maxTrainingTimeSteps-3,1);
modes = zeros(maxTrainingTimeSteps-3,1);
tainingVector = linspace(3,maxTrainingTimeSteps, maxTrainingTimeSteps-2);

trainingTimeSteps = 6;

% for trainingTimeSteps = 3:maxTrainingTimeSteps

numberOfElementsInX = timeSteps;
refinementDepth = 3;
numberOfPODModes = 0;

t = linspace(0, tEnd, timeSteps + 1);                                       % time discretization
x = linspace(0.0, xEnd, numberOfElementsInX + 1);                           % spatial discretization X
x_ref = linspace(0.0, 1.0, 2^refinementDepth + 1);                          % refinement discretization X_ref
x_PostProcess = linspace(0.0, xEnd, 5*(numberOfElementsInX +1));            % post-processing coordinates

[X, T] = meshgrid(x_PostProcess, t);
[X_ref, T_ref] = meshgrid(x_ref, t);


%% Analysis
[temperatureSolution, heatFlux, internalEnergy, modes(trainingTimeSteps-2)] = backwardEulerXFEM(x, x_PostProcess, rhs, T0,...
    dirichletLeftBC, dirichletRightBC, k, heatCapacity, t, refinementDepth, trainingTimeSteps, numberOfPODModes);

%% Post-Process
overkilledInternalEnergy = 4.221615373269894e+07;

relErrorEnergyNorm(trainingTimeSteps-2) = sqrt((internalEnergy(end)-overkilledInternalEnergy)^2)/sqrt(overkilledInternalEnergy^2);

figure(2)

surf(X, T, temperatureSolution')

figure(3)
surf(X, T, heatFlux')

% end
%
% figure(1)
% semilogy( tainingVector, relErrorEnergyNorm)
% figure(101)
% plot( tainingVector, modes)



figure(4)
F(size(t,2)) = struct('cdata',[],'colormap',[]);

for i=1:size(t,2)
    plot(x_PostProcess',temperatureSolution(:,i))
    drawnow
    F(i) = getframe;
    writeVideo(writerObj, getframe(gcf, [ 0 0 560 420 ]));
end

close(writerObj);