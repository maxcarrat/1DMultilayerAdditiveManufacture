function [ h ] = capacityPhaseTransition( x, T, TlastConv, specificHeat, rho, entalphyJump,...
    Tliq, Tmel)
%CAPACITYPHASETRANSITION returns the thermal capacity for the
%liquid and the solid part. K. C. Mills, "Recommended Values of
%Thermophysic Properties for selected commercial alloys", 2002.

h = zeros(numel(x),1);

for i=1:numel(x)
    
    %    if T > 0.0 &&  T < Tmel
    %        h(i) = rho * specificHeat;
    %    elseif T >= Tmel  &&  T <= Tliq
    %        h(i) = rho * specificHeat + rho * entalphyJump / (10);
    %    else
    %        h(i) = rho * specificHeat + rho * entalphyJump;
    %    end
    
    %     h(i) = rho * specificHeat + rho * entalphyJump / 2 / (Tliq-Tmel) *...
    %         (tanh( (T - Tmel) / 10) + 1);
    
    h(i) = rho * specificHeat + rho * entalphyJump * ...
        phaseChangeFuncDer(T, TlastConv, Tmel);% * (tanh( (T - Tmel) / 1) + 1);

    
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