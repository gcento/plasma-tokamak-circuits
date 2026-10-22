function [names, Ic] = circuit_currents(Ia, coil_names)
% The 12 circuit currents as the controller measures them, from the 14 coil currents
% Ia (one coil per column, named coil_names): CS1 = (CS1U + CS1L)/2,
% VS3 = (VS3U - VS3L)/2, the other circuits one to one.
coil = @(c) Ia(:, strcmp(coil_names, c));
names = {'CS3U','CS2U','CS1','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5','PF6','VS3'};
Ic = [coil('CS3U'), coil('CS2U'), (coil('CS1U') + coil('CS1L'))/2, coil('CS2L'), ...
    coil('CS3L'), coil('PF1'), coil('PF2'), coil('PF3'), coil('PF4'), coil('PF5'), ...
    coil('PF6'), (coil('VS3U') - coil('VS3L'))/2];
end
