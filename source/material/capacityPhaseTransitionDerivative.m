function [ hDerivative ] = capacityPhaseTransitionDerivative( x, T, TlastConv, specificHeat, rho, entalphyJump,...
    Tliq, Tmel )
%CAPACITYPHASETRANSITIONDERIVATIVE  returns the thermal capacity derivative for the
%liquid and the solid part. K. C. Mills, "Recommended Values of
%Thermophysic Properties for selected commercial alloys", 2002.

hDerivative = zeros(numel(x),1);

for i=1:numel(x)
    
    hDerivative(i) =  rho * entalphyJump * phaseChangeFuncSecondDer(T(i), TlastConv, Tmel) * (tanh( (T(i) - Tmel) / 1) + 1);

end

end

function f_pcDer2 = phaseChangeFuncSecondDer(T, TlastConvergent, Tmel)

if T == TlastConvergent
    f_pcDer2 = 0.0;
else
    TIncrement = abs(T-TlastConvergent);
    f_pcDer2 = (phaseChangeFunc(TlastConvergent+TIncrement, Tmel) -...
        2*phaseChangeFunc(TlastConvergent, Tmel) + phaseChangeFunc(TlastConvergent-TIncrement, Tmel))...
        / (T - TlastConvergent)^2;
end

end

% function f_pcDer = phaseChangeFuncDer(T, TlastConvergent, Tmel)
% 
% if T == TlastConvergent
%     f_pcDer = 0.0;
% else
%     f_pcDer = (phaseChangeFunc(T, Tmel) - phaseChangeFunc(TlastConvergent, Tmel))...
%         / (T - TlastConvergent);
% end
% 
% end

function f_pc = phaseChangeFunc(T, Tmel)

if T <= Tmel
    f_pc = 0.0;
else
    f_pc = 1.0;
end

end

