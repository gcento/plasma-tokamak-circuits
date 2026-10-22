function [same, msg] = same_model(dir_new, dir_ref)
% Check that a model folder reproduces a reference folder, as a generator must: the
% descriptor fields of linear_models/README.md equal and the matrices A, B, C equal
% within 1e-12 of the largest reference entry. Prints one line (msg).
%   same = same_model(dir_new, dir_ref)
%   [same, msg] = same_model(dir_new, dir_ref)
% Needs plant/ on the path (read_npy).

dn = jsondecode(fileread(fullfile(dir_new,'linear_model.json')));
dr = jsondecode(fileread(fullfile(dir_ref,'linear_model.json')));
differ = {};
for f = {'states','inputs','coil_states','passive_states','outputs', ...
        'coordinates','convention','operating_point','shot','time'}
    if isfield(dn,f{1}) ~= isfield(dr,f{1}) || ...
            (isfield(dn,f{1}) && ~isequaln(dn.(f{1}), dr.(f{1})))
        differ{end+1} = f{1}; %#ok<AGROW>
    end
end
for m = {'A','B','C'}
    Mn = read_npy(fullfile(dir_new, dn.files.(m{1})));
    Mr = read_npy(fullfile(dir_ref, dr.files.(m{1})));
    tol = 1e-12*max([abs(Mr(:)); realmin]);
    if ~isequal(size(Mn), size(Mr)) || any(abs(Mn - Mr) > tol, 'all')
        differ{end+1} = m{1}; %#ok<AGROW>
    end
end
same = isempty(differ);
if same
    msg = sprintf('%s reproduces %s', dir_new, dir_ref);
else
    msg = sprintf('%s differs from %s in: %s', dir_new, dir_ref, strjoin(differ, ', '));
end
fprintf('%s\n', msg);
end
