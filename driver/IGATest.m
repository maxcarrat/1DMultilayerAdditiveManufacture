%% IGA Test
% evaluate a 1D poisson transient problem using IGA
clear all;
clc;

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 7500.0;                                               % density [kg/m^3]
c = 600.0;                                                  % specific heat [J/(kg°C)]
T0 = 20.0;                                                  % Initial temperature [°C]
entalphyJump = 261e+03;                                     % entalphy [J/Kg]
Tsource = 1500.0;                                           % source temperature [°C]
Tmelt = 1450.0;                                             % melting temperature [°C]

tEnd = 0.5;                                                  % total time [sec]
xEnd = 0.1;                                                % length of the bar 0.001 [m]

dirichletLeftBC = @(t) 0.0;
dirichletRightBCValue = 0.0;
nuemannRightBC = 0.0e+04;
bodySource = 1.0e+02;

tolerance = 1.0e-03;
maxIteration = 100;

p = 2;                                                      % polynomial degree
integrationOrder = p+1;                                     % integration order

maxRefinementDepth = 1;                                     % max refinement depth
DOFs = zeros(maxRefinementDepth, 1);                        % DOFs vector to print

for depth = 1:maxRefinementDepth
    
    % bar specifics
    numberOfLayers = 1;
    trainingTimeSteps = numberOfLayers;
    numberOfTimeStepsPerLayer = 10;              % total time per layer 0.225 sec
    numberOfHeatingTimeSteps = 1;               % heating laser time per layer 0.090 sec
    
    numberOfElementsInX = numberOfLayers;       % one element per layer
    timeSteps = numberOfLayers;
    refinementDepth = depth;                    % refinement depth on the last layer
    
    numberOfPODModes = 0;
    PODRefinementDepth = 0;
    numberOfRefinedElementsToBeKept = 1;
    
    % Dirichlet condition on the bar tip
    dirichletRightBC = @(t) dirichletRightBCValue;
    
    % body heat source in the bar
    bodyLoad = @(x, t) bodySource;
    
    % conductivity function                                                     % thermal conductivity [W/(m K)]
    k = @(x, t, T) 1.0;
    
    % heat capacity function                                                    % heat capacity [kJ / kg °C]
    heatCapacity = @(x, T, T_last) 100.0;
    
    % discretization
    t = linspace(0, tEnd, numberOfTimeStepsPerLayer*numberOfLayers + 1);        % time discretization
    layerThickness = xEnd/numberOfLayers;                                       % thickness of the single layer
    
    TotalNumberOfControlPoints =  numberOfElementsInX + p;
    
    x_PostProcess = linspace(0.0, xEnd, 1000);        % post-processing coordinates
    
    [X, T] = meshgrid(x_PostProcess, t);
    
    
    %% Analysis
    [temperatureSolution, heatFlux, internalEnergy, CPUTime, DOFs(depth)] = ...
        IGAMultiPhaseBackwardEulerSolver(p, x_PostProcess, bodyLoad, T0,...
        dirichletLeftBC, dirichletRightBC, nuemannRightBC, k, heatCapacity, t, tolerance,...
        maxIteration,numberOfRefinedElementsToBeKept, TotalNumberOfControlPoints,...
        refinementDepth, PODRefinementDepth, trainingTimeSteps, numberOfTimeStepsPerLayer,...
        numberOfLayers, numberOfPODModes, integrationOrder, integrationOrder, layerThickness);
    
    %% Post-Process
    figure(depth+6)
    
    % Create axes
    axes1 = axes;
    
    for i=1:size(temperatureSolution, 2)
        
        % Create multiple lines using matrix input to plot
        plot(x_PostProcess,temperatureSolution(:,i)');
        
        % Title
        title('Temperature distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Temperature [°C]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes1,[0 xEnd]);
        ylim(axes1,[-1 1]);
        
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    % Create ylabel
    
    
    hold off
    
    figure(depth+1006)
    % Create axes
    axes2 = axes;
    
    for i=1:size(heatFlux, 2)
        
        % Create multiple lines using matrix input to plot
        plot(x_PostProcess,heatFlux(:,i)');
        
        % Title
        title('Temperature distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Temperature flux [°C/m^2]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes2,[0 xEnd]);
        ylim(axes2,[-1.0e+1 1.0e+1]);
        
        box(axes2,'on');
        % Set the remaining axes properties
        set(axes2,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end

end % depth loop


