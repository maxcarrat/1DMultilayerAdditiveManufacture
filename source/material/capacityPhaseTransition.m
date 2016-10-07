function [ h ] = capacityPhaseTransition( x, T, TlastConv, specificHeat, rho, entalphyJump,...
    Tliq, Tmel)
%CAPACITYPHASETRANSITION returns the thermal capacity for the
%liquid and the solid part. K. C. Mills, "Recommended Values of
%Thermophysic Properties for selected commercial alloys", 2002.

h = zeros(numel(x),1);

for i=1:numel(x)
    
    h(i) = rho * specificHeat + rho * entalphyJump * ...
        phaseChangeFuncDer(T, TlastConv, Tmel) * (tanh( (T - Tmel) / 1) + 1);

end

end

function f_pcDer = phaseChangeFuncDer(T, TlastConvergent, Tmel)

if T == TlastConvergent
    f_pcDer = 0.0;
else
    f_pcDer = (phaseChangeFunc(T, Tmel) - phaseChangeFunc(TlastConvergent, Tmel))...
        / (T - TlastConvergent);
end

end

function f_pc = phaseChangeFunc(T, Tmel)

if T <= Tmel
    f_pc = 0.0;
else
    f_pc = 1.0;
end

end