function [f, uuid] = tmatrix_hdf5(f, tmatrix, modes, wavelength, epsilon, geometry, computation, comments)
%tmatrix_hdf5 Write T-matrix into standard HDF5 format
%   Arguments somewhat tailored for SMARTIES
%
%  f: filename
%  tmatrix: array (wide format)
%  modes: struct with fields: l, m, polarization
%  wavelength: scalar or vector [in nm]
%  epsilon: struct with fields: embedding, [particle material]
%  geometry: struct with fields: description, [particle material]
%  computation: struct with fields: embedding, [particle material]
%  comments: struct with fields: script,
%
%

uuid = char(matlab.lang.internal.uuid());

% append file script
computation.script = fileread(comments.script);

% grab polarization as it's not handled correctly by easyh5
polarization = modes.polarization;
modes = rmfield(modes,'polarization');

% strip description (written as attribute, not value)
geo_description = geometry.description;
geometry = rmfield(geometry,'description');

embedding = struct('relative_permeability', 1.0,'relative_permittivity', epsilon.embedding);
matname = setdiff(fieldnames(epsilon),'embedding');
particle = struct('relative_permeability', 1.0,'relative_permittivity', epsilon.(matname{1}));
materials = struct('embedding', embedding, matname{1}, particle);


s = struct('tmatrix', tmatrix, ...
    'vacuum_wavelength', wavelength, ...
    'embedding', embedding,...
    'materials', materials, ...
    'geometry', geometry, ...
    'modes', modes, ...
    'computation', computation, ...
    'uuid', uuid);

saveh5(s, f, 'ComplexFormat', {'r','i'}, 'rootname', '', 'Compression', 'deflate'); 

% deal with polarization manually
h5create(f,'/modes/polarization', size(polarization), 'Datatype', 'string')
h5write(f,'/modes/polarization', polarization)

% attributes
h5writeatt(f, '/', 'name', comments.name);
h5writeatt(f, '/', 'created_with', 'Matlab easyh5');
[major,minor,rel] = H5.get_libversion();
h5writeatt(f, '/', 'storage_format_version', sprintf('%d.%d.%d',major,minor,rel));
h5writeatt(f, '/','description', comments.description);
h5writeatt(f, '/','keywords', comments.keywords);
h5writeatt(f, '/vacuum_wavelength', 'unit', 'nm');
h5writeatt(f, '/uuid', 'version', '4');
h5writeatt(f, '/geometry', 'name', geo_description);

end