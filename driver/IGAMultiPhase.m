%% 1D Multi-layer Non-Linear Multi-phase Additive Manufacturing
% evaluate a 1D poisson transient non-linear problem using IGA

clear all;
clc;

writerObj1 = VideoWriter('IGAMultiPhaseTemperature.avi');
writerObj1.Quality = 100;
writerObj1.FrameRate = 5;
open(writerObj1);

writerObj2 = VideoWriter('IGAMultiPhaseHeatFlux.avi');
writerObj2.Quality = 100;
writerObj2.FrameRate = 5;
open(writerObj2);

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 7500.0;                                               % density [kg/m^3]
c = 600.0;                                                  % specific heat [J/(kg°C)]
T0 = 20.0;                                                  % Initial temperature [°C]
entalphyJump = 261e+03;                                     % entalphy [J/Kg]
Tsource = 1500.0;                                           % source temperature [°C]
Tmelt = 1450.0;                                             % melting temperature [°C]

tEnd = 4.5;                                                 % total time [sec]
xEnd = 0.001;                                               % length of the bar 0.001 [m]

dirichletLeftBC = @(t) T0;
dirichletRightBCValue = T0 + Tsource;
nuemannRightBC = 0.0e+04;
bodySource = 0.0e+00;

tolerance = 1.0e-03;
maxIteration = 100;

p = 2;                                                      % polynomial degree
integrationOrder = p+1;                                     % integration order

maxRefinementDepth = 1;                                     % max refinement depth
DOFs = zeros(maxRefinementDepth, 1);                        % DOFs vector to print

for depth = 1:maxRefinementDepth
    
    % bar specifics
    numberOfLayers = 20;
    trainingTimeSteps = numberOfLayers;
    numberOfTimeStepsPerLayer = 10;     % total time per layer 0.225 sec
    numberOfHeatingTimeSteps = 4;       % heating laser time per layer 0.090 sec
    
    numberOfElementsInX = numberOfLayers;   % one element per layer
    timeSteps = numberOfLayers;
    refinementDepth = depth;                % refinement depth on the last layer
    
    if depth == 1
        PODRefinementDepth = 0;
    else
        PODRefinementDepth = floor( depth / 2 );
    end
    
    numberOfPODModes = 0;
    
    numberOfRefinedElementsToBeKept = 1;
    
    % Dirichlet condition on the bar tip
    dirichletRightBC = @(t) heatingCoolingBoundary(t, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, numberOfHeatingTimeSteps, dirichletRightBCValue);
    %     dirichletRightBC = @(t) dirichletRightBCValue;
    
    % body heat source in the bar
    bodyLoad = @(x, t) externalSource(x, t, xEnd, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, refinementDepth, bodySource);
    %     bodyLoad = @(x, t) bodySource;
    
    % conductivity function                                                     % thermal conductivity [W/(m K)]
    k = @(x, t, T) steelThermalConductivity(x, t, xEnd, numberOfLayers, tEnd,...
        numberOfTimeStepsPerLayer, T);
    %     k = @(x, t, T) 27.0;
    
    % heat capacity function                                                    % heat capacity [kJ / kg °C]
        heatCapacity= @(x, T, T_last)  capacityPhaseTransition(x, T, T_last, c, rho,...
            entalphyJump, Tsource, Tmelt);
%     heatCapacity = @(x, T, T_last) rho * c;
    
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
        ylim(axes1,[-100 3000]);
        
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
        ylabel('\fontname{Latin Modern Math} Heat flux [W/m^2]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes2,[0 xEnd]);
        ylim(axes2,[-1.0e+07 4.0e+08]);
        
        box(axes2,'on');
        % Set the remaining axes properties
        set(axes2,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    % Create ylabel
    
    
    hold off
    % Write results to a file
    formatSpec = 'myIGAMultiPhaseResultsFile_%d.txt';
    filename = sprintf(formatSpec,depth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(temperatureSolution, 1)
        fprintf(resultFile, '%d, %\t', x_PostProcess(i));
        for j=1:numberOfTimeStepsPerLayer
            fprintf(resultFile, '%d, %\t', temperatureSolution(i,end-(j-1)));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    % Write fluxes to a file
    formatSpec = 'myIGAMultiPhaseFluxesFile_%d.txt';
    filename = sprintf(formatSpec,depth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:size(heatFlux, 1)
        fprintf(resultFile, '%d, %\t', x_PostProcess(i));
        for j=1:numberOfTimeStepsPerLayer
            fprintf(resultFile, '%d, %\t', heatFlux(i,end-(j-1)));
        end
        fprintf(resultFile, '\n');
    end
    fclose(resultFile);
    
    % Write CPU time to a file
    formatSpec = 'myIGAMultiPhaseTimeFile_%d.txt';
    filename = sprintf(formatSpec,depth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:numel(CPUTime)
        fprintf(resultFile, '%d, %\t', CPUTime(i));
    end
    
    fclose(resultFile);
    
    % Write DOFs to a file
    formatSpec = 'myIGAMultiPhaseDOFsFile_%d.txt';
    filename = sprintf(formatSpec,depth);
    resultFile = fopen(filename, 'wt'); % Open for writing
    for i=1:numel(DOFs)
        fprintf(resultFile, '%d, %\t', DOFs(i));
    end
    
    fclose(resultFile);
    
    
end % depth loop


% figure(10101)
% 
% % Create axes
% axes1 = axes;
% 
% F(size(t,2)) = struct('cdata',[],'colormap',[]);
% 
% for i=1:size(t,2)
%     plot(x_PostProcess,temperatureSolution(:,i)')
% 
%     % Title
%     title('Temperature evolution of the bar', 'Interpreter','latex');
% 
%     % Create ylabel
%     ylabel('\fontname{Latin Modern Math} Temperature [°C]');
%     xlabel('\fontname{Latin Modern Math} length [m]');
% 
%     % Axis limit
%     xlim(axes1,[0 xEnd]);
%     ylim(axes1,[-100 3000]);
% 
%     box(axes1,'on');
%     % Set the remaining axes properties
%     set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
% 
%     grid on
%     drawnow
%     F(i) = getframe;
%     writeVideo(writerObj1, getframe(gcf, [ 0 0 500 400 ]));
% end
% close(writerObj1);
% 
% figure(10102)
% 
% % Create axes
% axes1 = axes;
% 
% F(size(t,2)) = struct('cdata',[],'colormap',[]);
% 
% for i=1:size(t,2)
%     plot(x_PostProcess,heatFlux(:,i)')
% 
%     % Title
%     title('Temperature evolution of the bar', 'Interpreter','latex');
% 
%     % Create ylabel
%     ylabel('\fontname{Latin Modern Math} Heat Flux [W/m^2]');
%     xlabel('\fontname{Latin Modern Math} length [m]');
% 
%     % Axis limit
%     xlim(axes1,[0 xEnd]);
%     ylim(axes1,[-1.0e+07 4.0e+08]);
% 
%     box(axes1,'on');
%     % Set the remaining axes properties
%     set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
% 
%     grid on
%     drawnow
%     F(i) = getframe;
%     writeVideo(writerObj2, getframe(gcf, [ 0 0 500 400 ]));
% end
% 
% close(writerObj2);