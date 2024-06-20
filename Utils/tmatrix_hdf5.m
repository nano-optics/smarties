function [f, uuid] = tmatrix_hdf5(f, tmatrix, modes, wavelength, epsilon, geometry, computation, comments)
%tmatrix_hdf5 Write T-matrix into standard HDF5 format
%   Arguments somewhat tailored for SMARTIES
%
%  f: filename
%  tmatrix: array (wide format)
%  modes: struct with fields: l, m, polarization: note must be ordered as
%  per standard
%  wavelength: scalar or vector [in nm]
%  epsilon: struct with fields: embedding, [particle material]
%  geometry: struct with fields: description, [particle material]
%  computation: struct with fields: embedding, [particle material]
%  comments: struct with fields: script,
%
%

uuid = char(matlab.lang.internal.uuid());

% append file script
% computation.files = struct('script', fileread(comments.script));

% remove script as it's not handled correctly by easyh5
script = convertCharsToStrings(fileread(comments.script));

% remove polarization as it's not handled correctly by easyh5
polarization = modes.polarization;
modes = rmfield(modes,'polarization');


embedding = struct('relative_permeability', 1.0,'relative_permittivity', epsilon.embedding);
matname = setdiff(fieldnames(epsilon),'embedding');
particle = struct('relative_permeability', 1.0,'relative_permittivity', epsilon.(matname{1}));
% materials = struct('embedding', embedding, matname{1}, particle);

scatterer = struct('material', particle, ...
                  'geometry', geometry);

s = struct('tmatrix', tmatrix, ...
    'vacuum_wavelength', wavelength, ...
    'embedding', embedding,...
    'scatterer', scatterer, ...
    'modes', modes, ...
    'computation', computation);

    % 'materials', materials, ...
    % 'geometry', geometry, ...

saveh5(s, f, 'ComplexFormat', {'r','i'}, 'rootname', '', 'Compression', 'deflate'); 


% deal with string objects manually

h5create(f,'/computation/files/script', size(script), 'Datatype', 'string')
h5write(f,'/computation/files/script', script)

h5create(f,'/modes/polarization', size(polarization), 'Datatype', 'string')
h5write(f,'/modes/polarization', polarization)


% root attributes
h5writeatt(f, '/', 'name', comments.name);
h5writeatt(f, '/', 'created_with', 'Matlab easyh5');
% [major,minor,rel] = H5.get_libversion(); %sprintf('%d.%d.%d',major,minor,rel)
h5writeatt(f, '/', 'storage_format_version', 'v0.01'); 
h5writeatt(f, '/','description', comments.description);
h5writeatt(f, '/','keywords', comments.keywords);

h5writeatt(f, '/computation', 'description', comments.description);
h5writeatt(f, '/computation', 'name', 'SMARTIES');
h5writeatt(f, '/computation', 'method', 'EBCM');
h5writeatt(f, '/computation', 'software', 'SMARTIES');
h5writeatt(f, '/computation', 'version', '1.1');

h5writeatt(f, '/embedding', 'description', 'constant refractive index');
h5writeatt(f, '/embedding', 'keywords', 'non-dispersive');
h5writeatt(f, '/embedding', 'name', comments.material_embedding);

h5writeatt(f, '/scatterer/material', 'name', comments.material_spheroid);
h5writeatt(f, '/scatterer/material', 'reference', comments.material_reference);

h5writeatt(f, '/scatterer/geometry', 'name', 'homogeneous spheroid with symmetry axis z');
h5writeatt(f, '/scatterer/geometry', 'unit', 'nm');
h5writeatt(f, '/scatterer/geometry', 'shape', 'spheroid')


h5writeatt(f, '/vacuum_wavelength', 'unit', 'nm');

% h5writeatt(f, '/uuid', 'version', '4'); % not needed


end