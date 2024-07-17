function [f, uuid] = tmatrix_hdf5(f, tmatrix, wavelength, epsilon, geometry, computation, comments)
%tmatrix_hdf5 Write T-matrix into standard HDF5 format
%   Arguments somewhat tailored for SMARTIES
%
%  f: filename
%  tmatrix: 3D array, first dimension is wavelength index
%  wavelength: scalar or vector [in nm]
%  epsilon: struct with fields: embedding, particle
%  geometry: struct with fields: description, shape, radiusxy, radiusz
%  computation: struct with fields: Lmax, Ntheta
%  comments: struct with fields: script


% check that sizes make sense
[a,b,c] = size(tmatrix);

a == length(wavelength);
b == c;
Lmax = computation.Lmax;
qmax = 2*(Lmax*(Lmax+1)+Lmax); % T-matrix size
b == qmax;

i=1;
for l=1:Lmax
    for m=-l:l
        for s=1:2
            modes.l(i) = int64(l);
            modes.m(i) = int64(m);
            modes.s(i) = int64(s);
            i=i+1;
        end
    end
end

polars = ["electric", "magnetic"];
polarization = polars(modes.s);
modes = rmfield(modes,'s');


% smarties does not do magnetic materials, and has only 2 regions
embedding = struct('relative_permeability', 1.0,'relative_permittivity', epsilon.embedding);
particle = struct('relative_permeability', 1.0,'relative_permittivity', epsilon.particle);

scatterer = struct('material', particle, ...
                  'geometry', geometry);


method_parameters = struct('Lmax', int64(computation.Lmax), ...
                            'Ntheta', int64(computation.Ntheta));
s = struct('tmatrix', tmatrix, ...
    'embedding', embedding,...
    'scatterer', scatterer, ...
    'computation', struct('method_parameters', method_parameters, ...
                          'analytical_zeros', computation.analytical_zeros));

saveh5(s, f, 'ComplexFormat', {'r','i'}, 'rootname', '', 'Compression', 'deflate'); 


% deal with objects manually
% character vectors, and also 1D vectors, 
% as easyh5 would store them as 1xN arrays for back-compatibility

h5create(f,'/vacuum_wavelength', length(wavelength), 'Datatype', 'double')
h5write(f,'/vacuum_wavelength', wavelength)


% modes
h5create(f,'/modes/l', length(modes.l), 'Datatype', 'int64')
h5write(f,'/modes/l', modes.l)

h5create(f,'/modes/m', length(modes.m), 'Datatype', 'int64')
h5write(f,'/modes/m', modes.m)

h5create(f,'/modes/polarization', length(polarization), 'Datatype', 'string')
h5write(f,'/modes/polarization', polarization)

script = convertCharsToStrings(fileread(comments.script));
h5create(f,'/computation/files/script', size(script), 'Datatype', 'string')
h5write(f,'/computation/files/script', script)


% root attributes
h5writeatt(f, '/', 'name', comments.name);
h5writeatt(f, '/', 'storage_format_version', 'v1'); 
h5writeatt(f, '/','description', comments.description);
h5writeatt(f, '/','keywords', comments.keywords);

% object and group attributes
h5writeatt(f, '/vacuum_wavelength', 'unit', 'nm');

h5writeatt(f, '/embedding', 'name', epsilon.embedding_name);
h5writeatt(f, '/embedding', 'reference', epsilon.embedding_reference);
h5writeatt(f, '/embedding', 'keywords', epsilon.embedding_keywords);

h5writeatt(f, '/scatterer/material', 'name', epsilon.material_name);
h5writeatt(f, '/scatterer/material', 'reference', epsilon.material_reference);
h5writeatt(f, '/scatterer/material', 'keywords', epsilon.material_keywords);

h5writeatt(f, '/scatterer/geometry', 'name', 'homogeneous spheroid with symmetry axis z');
h5writeatt(f, '/scatterer/geometry', 'unit', 'nm');
h5writeatt(f, '/scatterer/geometry', 'shape', 'spheroid')

[h5major,h5minor,h5rel] = H5.get_libversion(); % HDF5 version
matlabv = version ; % Matlab version
software = sprintf('SMARTIES=1.1, matlab=%s, HDF5=%d.%d.%d',matlabv,h5major,h5minor,h5rel);

h5writeatt(f, '/computation', 'description', computation.description);
h5writeatt(f, '/computation', 'accuracy', computation.accuracy);
h5writeatt(f, '/computation', 'name', 'SMARTIES');
h5writeatt(f, '/computation', 'method', 'EBCM, Extended Boundary Condition Method');
h5writeatt(f, '/computation', 'software', software);

% no longer needed, the database will add a uuid
% uuid = char(matlab.lang.internal.uuid());
% h5writeatt(f, '/uuid', 'version', '4'); % not needed


end