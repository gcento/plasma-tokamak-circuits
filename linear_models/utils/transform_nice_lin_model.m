close all
clear 
clc
% Unlike the raw plasmaless model, NICE already enforces circuit constraints in the
% underlying physics.However, NICE exports VS3 and -VS3 as two redundant
% states. Because the physical constraint is already into A and B the model
% can be transforemd via pure kinematic transformation:

repo      = "/home/ITER/centong/Desktop/PDSs/IMAS-PDS/local_installs/plasmaless-tokamak-circuits";
addpath("/home/ITER/centong/Desktop/PDSs/IMAS-PDS/local_installs/plasmaless-tokamak-circuits/imas_model");

% raw NICE models (linear_models/sources/<name>) -> circuit models (linear_models/<name>)
for name = ["niceLin_105073_t84.4176", "niceSimpleLin_105073_t84.4176"]

    model_dir = fullfile(repo, "linear_models", "sources", name);
    out_dir   = fullfile(repo, "linear_models", name);
    d = jsondecode(fileread(fullfile(model_dir, "linear_model.json")));

    A = read_npy(fullfile(model_dir, d.files.A));
    B = read_npy(fullfile(model_dir, d.files.B));
    C = read_npy(fullfile(model_dir, d.files.C));

    % nx size of the original state (13 + np)
    % np size of the passive structure
    % na size of the original active coils
    % nc size of the new state
    [~,nx] = size(A);
    np = nx - 13;
    na = 13;
    nc = 12;

    Ta = [  eye(11)  , zeros(11,2);
          zeros(1,11), 0.5, -0.5 ];

    Wa = [     eye(12)       ;
          zeros(1,11), -1   ];


    T = [     Ta     , zeros(nc,np);
         zeros(np,na),    eye(np)  ];

    W = [     Wa      , zeros(na,np);
          zeros(np,nc),    eye(np) ]; % pinv(T)

    Anew = T*A*W;
    Bnew = T*B;
    Cnew = C*W;

    % new model folder: A, B, C (.npy) and the descriptor
    if ~isfolder(out_dir), mkdir(out_dir); end
    write_npy(fullfile(out_dir, "A.npy"), Anew);
    write_npy(fullfile(out_dir, "B.npy"), Bnew);
    write_npy(fullfile(out_dir, "C.npy"), Cnew);

    d.states(strcmp(d.states, '-VS3')) = [];    % coil_states unchanged: VS3L "-VS3" = minus VS3
    d.files = struct('A', 'A.npy', 'B', 'B.npy', 'C', 'C.npy');
    d.notes{end+1} = ['Reduced to 12 circuit states (2026-10-21) from linear_models/sources/' ...
        char(name) ' by linear_models/utils/transform_nice_lin_model.m: VS3 = (VS3 - (-VS3))/2, ' ...
        'the VS3U current with VS3L = -VS3U; A = T*A0*pinv(T), B = T*B0, C = C0*pinv(T).'];
    fid = fopen(fullfile(out_dir, "linear_model.json"), "w");
    fprintf(fid, "%s\n", jsonencode(d, "PrettyPrint", true));
    fclose(fid);
    fprintf('%s: %d -> %d states\n', out_dir, nx, size(Anew,1));

end

function write_npy(filename, M)
% minimal .npy writer: float64, 2-D, Fortran order (MATLAB's column-major storage)
header = sprintf('{''descr'': ''<f8'', ''fortran_order'': True, ''shape'': (%d, %d), }', ...
    size(M,1), size(M,2));
pad = mod(-(10 + numel(header) + 1), 64);       % total header length a multiple of 64
header = [header repmat(' ',1,pad) newline];
fid = fopen(filename, 'w', 'ieee-le');
fwrite(fid, [147 double('NUMPY') 1 0], 'uint8');
fwrite(fid, numel(header), 'uint16');
fwrite(fid, header, 'char');
fwrite(fid, full(M(:)), 'float64');
fclose(fid);
end