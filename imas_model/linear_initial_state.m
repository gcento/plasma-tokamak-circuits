function x0 = linear_initial_state(model, i_coil)
% initial state of an absolute linear model (configure_linear_model, coordinates
% 'absolute') from the coil currents i_coil (model.n_coils, machine-description order):
% least squares on the coil (circuit) states, passive states 0. With the ITER circuits
% CS1 is the mean of CS1U and CS1L and VS3 half the difference of VS3U and VS3L (with
% the sign of the descriptor). Currents that violate the circuits (CS1U ~= CS1L,
% VS3U ~= -VS3L) are not representable: error plasmaless:circuits when the mismatch
% exceeds 1e-6 of the largest current, as linear_models/utils/vertical_circuit/
% map_coils_to_circuits.m.

i_coil = i_coil(:);
if numel(i_coil) ~= model.n_coils
    error('plasmaless:coilCurrents','%d coil currents for %d coils.', numel(i_coil), model.n_coils);
end
is_coil_state = any(model.x_to_coil ~= 0,1);
x0 = zeros(model.n_states,1);
x0(is_coil_state) = model.x_to_coil(:,is_coil_state) \ i_coil;
mismatch = max(abs(model.x_to_coil*x0 - i_coil));
if mismatch > 1e-6*max(1,max(abs(i_coil)))
    error('plasmaless:circuits', ['Coil currents not consistent with the circuits of the ', ...
        'model (series CS1U = CS1L, antiseries VS3U = -VS3L): max mismatch %.3g A.'], mismatch);
end
end
