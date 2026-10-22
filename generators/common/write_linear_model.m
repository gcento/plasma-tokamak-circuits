function write_linear_model(out_dir, A, B, C, d, overwrite)
% Write a model folder for the plant (plant/configure_linear_model.m): the continuous
% matrices A, B, C as A.npy, B.npy, C.npy and the descriptor struct d as
% linear_model.json (d.files is set here). The matrix files are read back with
% plant/read_npy.m and must equal the matrices exactly.
%   write_linear_model(out_dir, A, B, C, d)             % out_dir must not hold a model
%   write_linear_model(out_dir, A, B, C, d, overwrite)  % overwrite = true replaces it
% Used by the generators of generators/plasmaless, generators/vertical and
% generators/nice.

if nargin < 6 || isempty(overwrite)
    overwrite = false;
end
out_dir = char(out_dir);
desc_file = fullfile(out_dir,'linear_model.json');
if isfile(desc_file) && ~overwrite
    error('write_linear_model:exists', ['%s already holds a model: choose another ', ...
        'out_dir or set overwrite = true.'], out_dir);
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

d.files = struct('A','A.npy','B','B.npy','C','C.npy');
write_npy(fullfile(out_dir,'A.npy'), A);
write_npy(fullfile(out_dir,'B.npy'), B);
write_npy(fullfile(out_dir,'C.npy'), C);
fid = fopen(desc_file,'w');
if fid < 0
    error('write_linear_model:file','Cannot write %s.', desc_file);
end
fprintf(fid,'%s\n',jsonencode(d,'PrettyPrint',true));
fclose(fid);

% check: the matrix files read back as the plant reads them
names = {'A','B','C'};
mats = {A, B, C};
for k = 1:3
    back = read_npy(fullfile(out_dir,[names{k} '.npy']));
    if ~isequal(size(back), size(mats{k})) || ~isequal(back, full(mats{k}))
        error('write_linear_model:check','%s.npy read back differs from the matrix.', names{k});
    end
end
end
