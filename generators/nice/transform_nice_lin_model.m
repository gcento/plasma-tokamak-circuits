function transform_nice_lin_model(source_dir, out_dir, overwrite)
% Reduce a raw NICE linear model (generators/nice/sources/<name>: 13 coil states, VS3 and
% -VS3 separate) to the 12 ITER circuit states and write the model folder read by the
% plant (linear_models/<name>).
%   transform_nice_lin_model(source_dir, out_dir)
%   transform_nice_lin_model(source_dir, out_dir, overwrite)   % overwrite out_dir
% Launcher: utils/make_nice_models.m (both NICE models).
%
% Unlike the raw plasmaless model, NICE already enforces circuit constraints in the
% underlying physics.However, NICE exports VS3 and -VS3 as two redundant
% states. Because the physical constraint is already into A and B the model
% can be transforemd via pure kinematic transformation:

if nargin < 3 || isempty(overwrite)
    overwrite = false;
end
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(here));                 % generators/nice
addpath(fullfile(repo_root,'plant'), fullfile(repo_root,'generators','common'));

source_dir = char(source_dir);
[~, name, ext] = fileparts(source_dir);
name = [name ext];                                      % model names contain a dot
d = jsondecode(fileread(fullfile(source_dir, "linear_model.json")));

A = read_npy(fullfile(source_dir, d.files.A));
B = read_npy(fullfile(source_dir, d.files.B));
C = read_npy(fullfile(source_dir, d.files.C));

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
d.states(strcmp(d.states, '-VS3')) = [];    % coil_states unchanged: VS3L "-VS3" = minus VS3
d.notes{end+1} = ['Reduced to 12 circuit states (2026-10-21) from generators/nice/sources/' ...
    name ' by generators/nice/transform_nice_lin_model.m: VS3 = (VS3 - (-VS3))/2, ' ...
    'the VS3U current with VS3L = -VS3U; A = T*A0*pinv(T), B = T*B0, C = C0*pinv(T).'];
write_linear_model(out_dir, Anew, Bnew, Cnew, d, overwrite);
fprintf('transform_nice_lin_model: %s: %d -> %d states\n', out_dir, nx, size(Anew,1));
end
