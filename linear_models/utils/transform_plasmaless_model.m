close all
clear 
clc
%The plasmaless model treats every coil as independent so one need to
%modify the model to take in account for circuits. consider 
% x : original state (14 coils)
% x_new : wanted state ( 12 circuits)
% we can find a map x = W x_new    M W x_newdot + R W x_new = [eye(14);zeros(np,14)] V u
% W' M W x_newdot + W' R W x_new = W' [eye(14);zeros(np,14)] V u
% from this expression we can derive the new model or find the way to
% modify the original one

repo      = "/home/ITER/centong/Desktop/PDSs/IMAS-PDS/local_installs/plasmaless-tokamak-circuits";

% raw coil model (written by export_plasmaless_model.m) -> circuit model read by the plant
model_dir = fullfile(repo, "linear_models", "sources", "plasmaless_model");
out_dir   = fullfile(repo, "linear_models", "plasmaless_model");
addpath(fullfile(repo, "imas_model"));      % read_npy, load_plasmaless_md, get_coupling_matrix
d = jsondecode(fileread(fullfile(model_dir, "linear_model.json")));

A = read_npy(fullfile(model_dir,'A.npy'));
B = read_npy(fullfile(model_dir,'B.npy'));
C = read_npy(fullfile(model_dir,'C.npy'));

% M of the 14 coils + 102 loops from the machine description (data/md_dd4, needs IMAS-MATLAB)
em = load_plasmaless_md();
Maa = get_coupling_matrix(em, 'mutual_active_active',   'pf_active/coil',  'pf_active/coil');
Mea = get_coupling_matrix(em, 'mutual_passive_active',  'pf_passive/loop', 'pf_active/coil');
Mee = get_coupling_matrix(em, 'mutual_passive_passive', 'pf_passive/loop', 'pf_passive/loop');
M = [Maa Mea'; Mea Mee];

% nx size of the original state (na + np)
% np size of the passive structure
% na size of the original active coils
% nc size of the new state
na = 14;          
nc = 12;           
nx = size(A, 1);
np = nx - na;      

% x = W x_new
Wa = zeros(na, nc);
Wa(1:2, 1:2)   = eye(2);             % CS3U, CS2U
Wa(3:4, 3)     = [1; 1];             % CS1U = CS1, CS1L = CS1
Wa(5:12, 4:11) = eye(8);             % CS2L, CS3L, PF1-PF6
Wa(13:14, 12)  = [1; -1];            % VS3U = VS3, VS3L = -VS3

W = blkdiag(Wa, eye(np));

Wm = (W'*M*W)\(W'*M);

% mapping u_old = V u_new to feed B
% u_old = CS3U, CS2U,CS1U, CS1L, CS2L, CS3L, PF1, PF2, PF3, PF4, PF5, PF6, VS3U,VS3L
% u_new = CS3U, CS2U, CS1U, CS1L,CS2L, CS3L, PF1, PF2, PF3, PF4, PF5, VS1, PF6, VS3

V = zeros(14, 14);
V(1:11, 1:11) = eye(11);
V(8:9, 12)   = -1;
V(10:11, 12) =  1;
V(12, 13) = 1;
V(13, 14) = -0.5;                    % VS3: V(VS3L) - V(VS3U) = Vps (positive VS3 -> VS3L current up)
V(14, 14) =  0.5;


Anew = Wm*A*W;
Bnew = Wm*B*V;
Cnew = C * W;


if ~isfolder(out_dir), mkdir(out_dir); end
write_npy(fullfile(out_dir, "A.npy"), Anew);
write_npy(fullfile(out_dir, "B.npy"), Bnew);
write_npy(fullfile(out_dir, "C.npy"), Cnew);

new_coils = {'CS3U', 'CS2U', 'CS1', 'CS2L', 'CS3L', ...
             'PF1', 'PF2', 'PF3', 'PF4', 'PF5', 'PF6', 'VS3'}';

d.states = [new_coils; d.passive_states(:)];

d.coil_states.CS1U = 'CS1';
d.coil_states.CS1L = 'CS1';
d.coil_states.VS3U = 'VS3';
d.coil_states.VS3L = '-VS3';

d.inputs = {'CS3U', 'CS2U', 'CS1U', 'CS1L', 'CS2L', 'CS3L', ...
            'PF1', 'PF2', 'PF3', 'PF4', 'PF5', 'VS1', 'PF6', 'VS3'}';

d.description = strrep(d.description, '14 coils', '12 coil circuits (14 voltage inputs)');
d.files = struct('A', 'A.npy', 'B', 'B.npy', 'C', 'C.npy');
d.provenance = ['linear_models/sources/plasmaless_model reduced to 12 circuits by linear_models/utils/transform_plasmaless_model.m: ', ...
    'M from the machine description (data/md_dd4), Wm = (W''*M*W)\(W''*M), A = Wm*A0*W, B = Wm*B0*V, C = C0*W'];
d.notes{end+1} = ['Circuit states (2026-10-21): CS1U = CS1L = CS1 (series), VS3U = VS3, VS3L = -VS3 ', ...
    '(antiseries, VS3 is the VS3U current); inputs: the 14 power supplies, VS1 = -PF2 -PF3 +PF4 +PF5, ', ...
    'VS3 = V(VS3L) - V(VS3U).'];

fid = fopen(fullfile(out_dir, "linear_model.json"), "w");
fprintf(fid, "%s\n", jsonencode(d, "PrettyPrint", true));
fclose(fid);

function write_npy(filename, M)
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