%% Phase-Change Validation Benchmark
% test the code wrt an analytical solution for a non-linear phase-change
% problem

clear all;
clc;

%% Problem SetUp
% Define the problem parameter, the boundary conditions
% and the discretization.

rho = 4.51e+03;                                               % density [kg/m^3]
c = 520;                                                      % specific heat [J/(kg°C)]
Tsol = 1000.0;                                                % Initial temperature [°C]
entalphyJump = 325e+03;                                       % entalphy [J/Kg]
Tliq = 2000.0;
Tmel = 1670;

tEnd = 100.0;                                              % total time [sec]
xEnd = 0.1;                                                % length of the bar 0.1 [m]

dirichletLeftBC = @(t) Tliq;
nuemannRightBC = 0.0e+00;
bodySource = 0.0e+13;

tolerance = 1.0e-03;
maxIteration = 20;

maxRefinementDepth = 6;
DOFs = zeros(maxRefinementDepth, 1);

% allocate plots 
p_analytical = zeros(1,1);                                  
p_ABAQUS = zeros(1,1);
p_numerical = zeros(1,1);
pFlux_ABAQUS = zeros(1,1);
pFlux_numerical = zeros(1,1);

for depth = 6:maxRefinementDepth
    
    integrationOrder = 2;
    numberOfElementsInX = 10;
    
    timeSteps = 10;
    refinementDepth = depth;
    
    dirichletRightBC = @(t) [];                     % No RHS Dirichlet condition
    
    bodyLoad = @(x, t) 0.0;
    
    k = @(x, t, T) 16;                            % constant thermal conductivity [W/(m °C)]
    
    heatCapacity= @(x, T, Tlast) capacityPhaseTransition(x, T, Tlast, c,...
        rho, entalphyJump, Tliq, Tmel);                                        % heat capacity [kJ / kg °C]
    % heat capacity [J / kg °C]
    
    t = linspace(0, tEnd, timeSteps + 1);        % time discretization
    x = linspace(0.0, xEnd, numberOfElementsInX + 1);                           % spatial discretization X

    numberOfElementsPostProcess = 150;
    x_PostProcess = linspace(0.0, xEnd, numberOfElementsPostProcess);          % post-processing coordinates
    
    [X, T] = meshgrid(x_PostProcess, t);
    
    
    %% Analysis
    [temperatureSolution, heatFlux, analyticalSolution, CPUTime, DOFs(depth)] =...
        multiPhaseBackwardEulerSolverNoAM(x, x_PostProcess, bodyLoad, Tsol,...
        dirichletLeftBC, dirichletRightBC, nuemannRightBC, k, heatCapacity, t, tolerance,...
        refinementDepth, integrationOrder, maxIteration);
    
    
    %% Post-Process
    
    figure(depth)
    % Create axes
    axes1 = axes;
    
    for i=1:size(temperatureSolution, 2)
        
        % Create multiple lines using matrix input to plot
        p_analytical(i) = plot(x_PostProcess,analyticalSolution(:,i)','r');
        
        % Title
        title('Temperature distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Temperature [°C]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes1,[0 xEnd]);
        ylim(axes1,[1000 2000]);
        
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
abaqusSolution = dlmread('abaqusFluxes.txt');

    for i=1:size(temperatureSolution, 2)-2
        
        % Create multiple lines using matrix input to plot
        p_ABAQUS(i) = plot(linspace(0.0, xEnd, size(abaqusSolution,1)),...
            abaqusSolution(1:end,i+11)','--b');
        
        % Title
        title('Temperature distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Temperature flux [°C/m^2]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes1,[0 xEnd]);
        ylim(axes1,[1000 2000]);
        
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    
    for i=1:size(temperatureSolution, 2)
        
        % Create multiple lines using matrix input to plot
        p_numerical(i) = plot(x_PostProcess,temperatureSolution(:,i)','-.k');
        
        % Title
        title('Temperature distibution at each time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Temperature [°C]');
        xlabel('\fontname{Latin Modern Math} length [m]');

        % Axis limit
        xlim(axes1,[0 xEnd]);
        ylim(axes1,[1000 2000]);
        
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    legend([p_analytical(end),p_ABAQUS(end),p_numerical(end)],'Analytical','ABAQUS', 'Numerical')
    
    hold off
    
    
    figure(depth+100)
    % Create axes
    axes2 = axes;
    
    for i=1:1%size(heatFlux, 2)
        
        % Create multiple lines using matrix input to plot
        pFlux_numerical(1)= plot(x_PostProcess,heatFlux(:,size(heatFlux, 2)-1)','--r');
        
        % Title
        title('Heat flux distibution at the last time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Heat flux [W/m^2]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes2,[0 xEnd]);
        ylim(axes2,[-8.0e+05 0.0e+01]);
        
        box(axes2,'on');
        % Set the remaining axes properties
        set(axes2,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    abaqusFluxSolution = dlmread('abaqusFluxes.txt');

    for i=1:1%size(heatFlux, 2)
        
        % Create multiple lines using matrix input to plot
        pFlux_ABAQUS(1) = plot(linspace(0.0, xEnd, size(abaqusFluxSolution,1)),...
            abaqusFluxSolution(1:end,6)','-.b');
        
        % Title
        title('Heat flux at the last time step', 'Interpreter','latex');
        
        % Axis labels
        ylabel('\fontname{Latin Modern Math} Heat flux [W/m^2]');
        xlabel('\fontname{Latin Modern Math} length [m]');
        
        % Axis limit
        xlim(axes2,[0 xEnd]);
        ylim(axes2,[-8.0e+05 0.0e+01]);
        
        box(axes2,'on');
        % Set the remaining axes properties
        set(axes2,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
        
        hold on
        
    end
    
    legend([pFlux_ABAQUS(end),pFlux_numerical(end)],'ABAQUS', 'Numerical')
   
    hold off
 
    
%     
%     hold off
%     % Write results to a file
%     formatSpec = 'myFEMMultiPhaseResultsFile_%d.txt';
%     filename = sprintf(formatSpec,depth);
%     resultFile = fopen(filename, 'wt'); % Open for writing
%     for i=1:size(temperatureSolution, 1)
%         fprintf(resultFile, '%d, %\t', x_PostProcess(i));
%         for j=1:timeSteps
%             fprintf(resultFile, '%d, %\t', temperatureSolution(i,end-(j-1)));
%         end
%         fprintf(resultFile, '\n');
%     end
%     fclose(resultFile);
%     
%     % Write fluxes to a file
%     formatSpec = 'myFEMMultiPhaseFluxesFile_%d.txt';
%     filename = sprintf(formatSpec,depth);
%     resultFile = fopen(filename, 'wt'); % Open for writing
%     for i=1:size(heatFlux, 1)
%         fprintf(resultFile, '%d, %\t', x_PostProcess(i));
%         for j=1:timeSteps
%             fprintf(resultFile, '%d, %\t', heatFlux(i,end-(j-1)));
%         end
%         fprintf(resultFile, '\n');
%     end
%     fclose(resultFile);
%     
%     % Write CPU time to a file
%     formatSpec = 'myFEMMultiPhaseTimeFile_%d.txt';
%     filename = sprintf(formatSpec,depth);
%     resultFile = fopen(filename, 'wt'); % Open for writing
%     for i=1:numel(CPUTime)
%         fprintf(resultFile, '%d, %\t', CPUTime(i));
%     end
%     
%     fclose(resultFile);
%     
%     % Write DOFs to a file
%     formatSpec = 'myXFEMMultiPhaseDOFsFile_%d.txt';
%     filename = sprintf(formatSpec,depth);
%     resultFile = fopen(filename, 'wt'); % Open for writing
%     for i=1:numel(DOFs)
%         fprintf(resultFile, '%d, %\t', DOFs(i));
%     end
%     
%     fclose(resultFile);
%     
%     
% end % depth loop
% 
% 
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
%     ylim(axes1,[1000 2000]);
%     
%     box(axes1,'on');
%     % Set the remaining axes properties
%     set(axes1,'FontSize',15,'FontWeight','normal', 'TickLabelInterpreter','latex');
%     
%     grid on
%     
end

