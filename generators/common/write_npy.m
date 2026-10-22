function write_npy(filename, M)
% minimal .npy writer: float64, 2-D, Fortran order (MATLAB's column-major storage);
% plant/read_npy.m reads it back
header = sprintf('{''descr'': ''<f8'', ''fortran_order'': True, ''shape'': (%d, %d), }', ...
    size(M,1), size(M,2));
pad = mod(-(10 + numel(header) + 1), 64);       % total header length a multiple of 64
header = [header repmat(' ',1,pad) newline];
fid = fopen(filename, 'w', 'ieee-le');
if fid < 0
    error('write_npy:file', 'Cannot write %s.', filename);
end
fwrite(fid, [147 double('NUMPY') 1 0], 'uint8');
fwrite(fid, numel(header), 'uint16');
fwrite(fid, header, 'char');
fwrite(fid, full(M(:)), 'float64');
fclose(fid);
end
