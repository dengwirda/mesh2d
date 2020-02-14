function savemsh(name,mesh)
%SAVEMSH save a *.MSH file for JIGSAW.
%
%   SAVEMSH(NAME,MESH);
%
%   The following are optionally written to "NAME.MSH". Ent-
%   ities are written if they are present in MESH:
%
%   .IF. MESH.MSHID == 'EUCLIDEAN-MESH':
%   -----------------------------------
%
%   MESH.POINT.COORD - [NPxND+1] array of point coordinates,
%       where ND is the number of spatial dimenions.
%       COORD(K,ND+1) is an ID tag for the K-TH point.
%
%   MESH.POINT.POWER - [NPx 1] array of vertex "weights",
%       associated with the dual "power" tessellation.
%
%   MESH.EDGE2.INDEX - [N2x 3] array of indexing for EDGE-2
%       elements, where INDEX(K,1:2) is an array of
%       "point-indices" associated with the K-TH edge, and
%       INDEX(K,3) is an ID tag for the K-TH edge.
%
%   MESH.TRIA3.INDEX - [N3x 4] array of indexing for TRIA-3
%       elements, where INDEX(K,1:3) is an array of
%       "point-indices" associated with the K-TH tria, and
%       INDEX(K,4) is an ID tag for the K-TH tria.
%
%   MESH.QUAD4.INDEX - [N4x 5] array of indexing for QUAD-4
%       elements, where INDEX(K,1:4) is an array of
%       "point-indices" associated with the K-TH quad, and
%       INDEX(K,5) is an ID tag for the K-TH quad.
%
%   MESH.TRIA4.INDEX - [M4x 5] array of indexing for TRIA-4
%       elements, where INDEX(K,1:4) is an array of
%       "point-indices" associated with the K-TH tria, and
%       INDEX(K,5) is an ID tag for the K-TH tria.
%
%   MESH.HEXA8.INDEX - [M8x 9] array of indexing for HEXA-8
%       elements, where INDEX(K,1:8) is an array of
%       "point-indices" associated with the K-TH hexa, and
%       INDEX(K,9) is an ID tag for the K-TH hexa.
%
%   MESH.WEDG6.INDEX - [M6x 7] array of indexing for WEDG-6
%       elements, where INDEX(K,1:6) is an array of
%       "point-indices" associated with the K-TH wedg, and
%       INDEX(K,7) is an ID tag for the K-TH wedg.
%
%   MESH.PYRA5.INDEX - [M5x 6] array of indexing for PYRA-5
%       elements, where INDEX(K,1:5) is an array of
%       "point-indices" associated with the K-TH pyra, and
%       INDEX(K,6) is an ID tag for the K-TH pyra.
%
%   MESH.BOUND.INDEX - [NBx 3] array of "boundary" indexing
%       in the domain, indicating how elements in the
%       geometry are associated with various enclosed areas
%       /volumes, herein known as "parts". INDEX(:,1) is an
%       array of "part" ID's, INDEX(:,2) is an array of
%       element numbering and INDEX(:,3) is an array of
%       element "tags", describing which element "kind" is
%       numbered via INDEX(:,2). Element tags are defined
%       via a series of constants instantiated in LIBDATA.
%       In the default case, where BOUND is not specified,
%       all elements in the geometry are assumed to define
%       the boundaries of enclosed "parts".
%
%   MESH.VALUE - [NPxNV] array of "values" associated with
%       the vertices of the mesh. NV values are associated
%       with each vertex.
%
%   MESH.SLOPE - [NPx 1] array of "slopes" associated with
%       the vertices of the mesh. Slope values define the
%       gradient-limits ||dh/dx|| used by the Eikonal solver
%       MARCHE.
%
%   .IF. MESH.MSHID == 'ELLIPSOID-MESH':
%   -----------------------------------
%
%   MESH.RADII - [ 3x 1] array of principle ellipsoid radii.
%
%   Additionally, entities described in the 'EUCLIDEAN-MESH'
%   data-type are optionally written.
%
%   .IF. MESH.MSHID == 'EUCLIDEAN-GRID':
%   .OR. MESH.MSHID == 'ELLIPSOID-GRID':
%   -----------------------------------
%
%   MESH.POINT.COORD - [NDx1] cell array of grid coordinates
%       where ND is the number of spatial dimenions. Each
%       array COORD{ID} should be a vector of grid coord.'s,
%       increasing or decreasing monotonically.
%
%   MESH.VALUE - [NMxNV] array of "values" associated with
%       the vertices of the grid, where NM is the product of
%       the dimensions of the grid. NV values are associated
%       with each vertex.
%
%   MESH.SLOPE - [NMx 1] array of "slopes" associated with
%       the vertices of the grid, where NM is the product of
%       the dimensions of the grid. Slope values define the
%       gradient-limits ||dh/dx|| used by the Eikonal solver
%       MARCHE.
%
%   See also JIGSAW, LOADMSH
%

%-----------------------------------------------------------
%   Darren Engwirda
%   github.com/dengwirda/jigsaw-matlab
%   20-Aug-2019
%   darren.engwirda@columbia.edu
%-----------------------------------------------------------
%

   [ok] = certify(mesh);

    if (~ischar  (name))
        error('NAME must be a valid file-name!') ;
    end
    if (~isstruct(mesh))
        error('MESH must be a valid structure!') ;
    end

   [path,file,fext] = fileparts(name) ;

    if(~strcmp(lower(fext),'.msh'))
        name = [name,'.msh'];
    end

    try
%-- try to write data to file

    ffid = fopen(name, 'w') ;

    nver = +3;

    if (exist('OCTAVE_VERSION','builtin')  > 0)
    fprintf(ffid,[ ...
    '# %s.msh; created by JIGSAW''s OCTAVE interface\n'],file) ;
    else
    fprintf(ffid,[ ...
    '# %s.msh; created by JIGSAW''s MATLAB interface\n'],file) ;
    end

    if (isfield(mesh,'mshID'))
        mshID =  mesh.mshID ;
    else
        mshID = 'EUCLIDEAN-MESH';
    end

    switch (upper(mshID))

    case 'EUCLIDEAN-MESH'
        save_mesh_format( ...
            ffid,nver,mesh,'EUCLIDEAN-MESH') ;
    case 'EUCLIDEAN-GRID'
        save_grid_format( ...
            ffid,nver,mesh,'EUCLIDEAN-GRID') ;
    case 'EUCLIDEAN-DUAL'
       %save_dual_format( ...
       %    ffid,nver,mesh,'EUCLIDEAN-DUAL') ;

    case 'ELLIPSOID-MESH'
        save_mesh_format( ...
            ffid,nver,mesh,'ELLIPSOID-MESH') ;
    case 'ELLIPSOID-GRID'
        save_grid_format( ...
            ffid,nver,mesh,'ELLIPSOID-GRID') ;
    case 'ELLIPSOID-DUAL'
       %save_dual_format( ...
       %    ffid,nver,mesh,'ELLIPSOID-DUAL') ;

    otherwise
        error('Invalid mshID!') ;

    end

    fclose(ffid);

    catch err

%-- ensure that we close the file regardless!
    if (ffid>-1)
    fclose(ffid) ;
    end
    rethrow(err) ;

    end

end

function save_mesh_format(ffid,nver,mesh,kind)
%SAVE-MESH-FORMAT save mesh data in unstructured-mesh format

    switch (upper(kind))
    case 'EUCLIDEAN-MESH'
            fprintf( ...
        ffid,'MSHID=%u;EUCLIDEAN-MESH\n',nver) ;

    case 'ELLIPSOID-MESH'
            fprintf( ...
        ffid,'MSHID=%u;ELLIPSOID-MESH\n',nver) ;

    end

    npts = +0 ;

    if (isfield(mesh,'radii') && ...
            ~isempty(mesh.radii) )

    %------------------------------------ write "RADII" data

        if (~isnumeric(mesh.radii))
            error('Incorrect input types');
        end
        if (ndims(mesh.radii) ~= 2)
            error('Incorrect dimensions!');
        end
        if (numel(mesh.radii) ~= 3)
            mesh.radii = mesh.radii(1) * ones(+3,+1);
        end

        fprintf(ffid,'RADII=%f;%f;%f\n',mesh.radii');

    end

    if (isfield(mesh,'point') && ...
            isfield(mesh.point,'coord') && ...
                ~isempty(mesh.point.coord) )

    %------------------------------------ write "POINT" data

        if (~isnumeric(mesh.point.coord))
            error('Incorrect input type!');
        end
        if (ndims(mesh.point.coord) ~= 2)
            error('Incorrect dimensions!');
        end

        ndim = size(mesh.point.coord,2)-1 ;
        npts = size(mesh.point.coord,1)-0 ;
        fprintf(ffid,['NDIMS=%u','\n'],ndim);

        fprintf(ffid, ...
            ['POINT=%u','\n'],size(mesh.point.coord,1));

        if (isa(mesh.point.coord,'double'))
            vstr = sprintf('%%1.%ug;',+16);
        else
            vstr = sprintf('%%1.%ug;',+ 8);
        end

        fprintf(ffid, ...
        [repmat(vstr,1,ndim),'%i\n'],mesh.point.coord');

    end

    if (isfield(mesh,'point') && ...
            isfield(mesh.point,'power') && ...
                ~isempty(mesh.point.power) )

    %------------------------------------ write "POWER" data

        if (~isnumeric(mesh.point.power))
            error('Incorrect input type!');
        end
        if (ndims(mesh.point.power) ~= 2)
            error('Incorrect dimensions!');
        end

        npwr = size(mesh.point.power,2) - 0 ;
        nrow = size(mesh.point.power,1) - 0 ;

        if (isa(mesh.point.power,'double'))
        vstr = sprintf('%%1.%ug;',+16) ;
        else
        vstr = sprintf('%%1.%ug;',+ 8) ;
        end
        vstr = repmat(vstr,+1,npwr) ;

        fprintf(ffid,['POWER=%u;%u','\n'],[nrow,npwr]);
        fprintf(ffid, ...
            [vstr(+1:end-1), '\n'], mesh.point.power');

    end

    if (isfield(mesh,'value'))

    %------------------------------------ write "VALUE" data

        if (~isnumeric(mesh.value))
            error('Incorrect input type!');
        end
        if (ndims(mesh.value) ~= 2)
            error('Incorrect dimensions!');
        end
        if (size(mesh.value,1) ~= npts)
            error('Incorrect dimensions!');
        end

        nrow = size(mesh.value,1);
        nval = size(mesh.value,2);

        if     (isa(mesh.value, 'double'))
        vstr = sprintf('%%1.%ug;',+16) ;
        elseif (isa(mesh.value, 'single'))
        vstr = sprintf('%%1.%ug;',+ 8) ;
        elseif (isa(mesh.value,'integer'))
        vstr = '%d;' ;
        else
            error('Incorrect input type!');
        end
        vstr = repmat(vstr,+1,nval) ;

        fprintf(ffid,['VALUE=%u;%u','\n'],[nrow,nval]);
        fprintf(ffid,[vstr(1:end-1),'\n'],mesh.value');

    end

    if (isfield(mesh,'slope'))

    %------------------------------------ write "SLOPE" data

        if (~isnumeric(mesh.slope))
            error('Incorrect input type!');
        end
        if (ndims(mesh.slope) ~= 2)
            error('Incorrect dimensions!');
        end
        if (size(mesh.slope,1) ~= npts && ...
            size(mesh.slope,1) ~= +1 )
            error('Incorrect dimensions!');
        end

        nrow = size(mesh.slope,1);
        nval = size(mesh.slope,2);

        if     (isa(mesh.slope, 'double'))
        vstr = sprintf('%%1.%ug;',+16) ;
        elseif (isa(mesh.slope, 'single'))
        vstr = sprintf('%%1.%ug;',+ 8) ;
        elseif (isa(mesh.slope,'integer'))
        vstr = '%d;' ;
        else
            error('Incorrect input type!');
        end
        vstr = repmat(vstr,+1,nval) ;

        fprintf(ffid,['SLOPE=%u;%u','\n'],[nrow,nval]);
        fprintf(ffid,[vstr(1:end-1),'\n'],mesh.slope');

    end

    if (isfield(mesh,'edge2') && ...
            isfield(mesh.edge2,'index') && ...
                ~isempty(mesh.edge2.index) )

    %------------------------------------ write "EDGE2" data

        if (~isnumeric(mesh.edge2.index))
            error('Incorrect input type!');
        end
        if (ndims(mesh.edge2.index) ~= 2)
            error('Incorrect dimensions!');
        end

        index = mesh.edge2.index;

        if (min(min(index(:,1:2))) < +1 || ...
            max(max(index(:,1:2))) > npts)
            error('Invalid EDGE-2 indexing!') ;
        end

        index(:,1:2) = ...
        index(:,1:2)-1 ; % zero-indexing!

        fprintf( ...
        ffid,['EDGE2=%u','\n'],size(index,1)) ;
        fprintf(ffid, ...
        [repmat('%u;',1,2),'%i','\n'],index') ;

    end

    if (isfield(mesh,'tria3') && ...
            isfield(mesh.tria3,'index') && ...
                ~isempty(mesh.tria3.index) )

    %------------------------------------ write "TRIA3" data

        if (~isnumeric(mesh.tria3.index))
            error('Incorrect input type!');
        end
        if (ndims(mesh.tria3.index) ~= 2)
            error('Incorrect dimensions!');
        end

        index = mesh.tria3.index;

        if (min(min(index(:,1:3))) < +1 || ...
            max(max(index(:,1:3))) > npts)
            error('Invalid TRIA-3 indexing!') ;
        end

        index(:,1:3) = ...
        index(:,1:3)-1 ; % zero-indexing!

        fprintf( ...
        ffid,['TRIA3=%u','\n'],size(index,1)) ;

        fprintf(ffid, ...
        [repmat('%u;',1,3),'%i','\n'],index') ;

    end

    if (isfield(mesh,'quad4') && ...
            isfield(mesh.quad4,'index') && ...
                ~isempty(mesh.quad4.index) )

    %------------------------------------ write "QUAD4" data

        if (~isnumeric(mesh.quad4.index))
            error('Incorrect input type!');
        end
        if (ndims(mesh.quad4.index) ~= 2)
            error('Incorrect dimensions!');
        end

        index = mesh.quad4.index;

        if (min(min(index(:,1:4))) < +1 || ...
            max(max(index(:,1:4))) > npts)
            error('Invalid QUAD-4 indexing!') ;
        end

        index(:,1:4) = ...
        index(:,1:4)-1 ; % zero-indexing!

        fprintf( ...
        ffid,['QUAD4=%u','\n'],size(index,1)) ;

        fprintf(ffid, ...
        [repmat('%u;',1,4),'%i','\n'],index') ;

    end

    if (isfield(mesh,'tria4') && ...
            isfield(mesh.tria4,'index') && ...
                ~isempty(mesh.tria4.index) )

    %------------------------------------ write "TRIA4" data

        if (~isnumeric(mesh.tria4.index))
            error('Incorrect input type!');
        end
        if (ndims(mesh.tria4.index) ~= 2)
            error('Incorrect dimensions!');
        end

        index = mesh.tria4.index;

        if (min(min(index(:,1:4))) < +1 || ...
            max(max(index(:,1:4))) > npts)
            error('Invalid TRIA-4 indexing!') ;
        end

        index(:,1:4) = ...
        index(:,1:4)-1 ; % zero-indexing!

        fprintf( ...
        ffid,['TRIA4=%u','\n'],size(index,1)) ;

        fprintf(ffid, ...
        [repmat('%u;',1,4),'%i','\n'],index') ;

    end

    if (isfield(mesh,'hexa8') && ...
            isfield(mesh.hexa8,'index') && ...
                ~isempty(mesh.hexa8.index) )

    %------------------------------------ write "HEXA8" data

        if (~isnumeric(mesh.hexa8.index))
            error('Incorrect input type!');
        end
        if (ndims(mesh.hexa8.index) ~= 2)
            error('Incorrect dimensions!');
        end

        index = mesh.hexa8.index;

        if (min(min(index(:,1:8))) < +1 || ...
            max(max(index(:,1:8))) > npts)
            error('Invalid HEXA-8 indexing!') ;
        end

        index(:,1:8) = ...
        index(:,1:8)-1 ; % zero-indexing!

        fprintf( ...
        ffid,['HEXA8=%u','\n'],size(index,1)) ;

        fprintf(ffid, ...
        [repmat('%u;',1,8),'%i','\n'],index') ;

    end

    if (isfield(mesh,'wedg6') && ...
            isfield(mesh.wedg6,'index') && ...
                ~isempty(mesh.wedg6.index) )

    %------------------------------------ write "WEDG6" data

        if (~isnumeric(mesh.wedg6.index))
            error('Incorrect input type!');
        end
        if (ndims(mesh.wedg6.index) ~= 2)
            error('Incorrect dimensions!');
        end

        index = mesh.wedg6.index ;

        if (min(min(index(:,1:6))) < +1 || ...
            max(max(index(:,1:6))) > npts)
            error('Invalid WEDG-6 indexing!') ;
        end

        index(:,1:6) = ...
        index(:,1:6)-1 ; % zero-indexing!

        fprintf( ...
        ffid,['WEDG6=%u','\n'],size(index,1)) ;

        fprintf(ffid, ...
        [repmat('%u;',1,6),'%i','\n'],index') ;

    end

    if (isfield(mesh,'pyra5') && ...
            isfield(mesh.pyra5,'index') && ...
                ~isempty(mesh.pyra5.index) )

    %------------------------------------ write "PYRA5" data

        if (~isnumeric(mesh.pyra5.index))
            error('Incorrect input type!');
        end
        if (ndims(mesh.pyra5.index) ~= 2)
            error('Incorrect dimensions!');
        end

        index = mesh.pyra5.index;

        if (min(min(index(:,1:5))) < +1 || ...
            max(max(index(:,1:5))) > npts)
            error('Invalid PYRA-5 indexing!') ;
        end

        index(:,1:5) = ...
        index(:,1:5)-1 ; % zero-indexing!

        fprintf( ...
        ffid,['PYRA5=%u','\n'],size(index,1)) ;

        fprintf(ffid, ...
        [repmat('%u;',1,6),'%i','\n'],index') ;

    end

    if (isfield(mesh,'bound') && ...
            isfield(mesh.bound,'index') && ...
                ~isempty(mesh.bound.index) )

    %------------------------------------ write "BOUND" data

        if (~isnumeric(mesh.bound.index))
            error('Incorrect input type!');
        end
        if (ndims(mesh.bound.index) ~= 2)
            error('Incorrect dimensions!');
        end

        index = mesh.bound.index ;
        index(:,2:2) = ...
        index(:,2:2)-1 ; % zero-indexing!

        fprintf( ...
        ffid,['BOUND=%u','\n'],size(index,1)) ;

        fprintf(ffid, ...
        [repmat('%u;',1,2),'%u','\n'],index') ;

    end

end

function save_grid_format(ffid,nver,mesh,kind)
%SAVE-GRID-FORMAT save mesh class in rectilinear-grid format

    switch (upper(kind))
    case 'EUCLIDEAN-GRID'
            fprintf( ...
        ffid,'MSHID=%u;EUCLIDEAN-GRID\n',nver) ;

    case 'ELLIPSOID-GRID'
            fprintf( ...
        ffid,'MSHID=%u;ELLIPSOID-GRID\n',nver) ;

    end

    dims = [] ;

    if (isfield(mesh,'radii') && ...
            ~isempty(mesh.radii) )

    %------------------------------------ write "RADII" data

        if (~isnumeric(mesh.radii))
            error('Incorrect input types');
        end
        if (ndims(mesh.radii) ~= 2)
            error('Incorrect dimensions!');
        end
        if (numel(mesh.radii) ~= 3)
            mesh.radii = mesh.radii(1) * ones(+3,+1);
        end

        fprintf(ffid,'RADII=%f;%f;%f\n',mesh.radii');

    end

    if (isfield(mesh,'point') && ...
            isfield(mesh.point,'coord') && ...
                ~isempty(mesh.point.coord) )

    %------------------------------------ write "COORD" data

        if(~iscell(mesh.point.coord) )
            error('Incorrect input types') ;
        end
        if ( numel(mesh.point.coord) ~= ...
            length(mesh.point.coord) )
            error('Incorrect dimensions!') ;
        end

        ndim = length(mesh.point.coord);
        dims = zeros(1,ndim);
        iord = [2,1,+3:ndim];

        fprintf(ffid, ...
        ['NDIMS=%u \n'],length(mesh.point.coord));

        for idim = +1 : length(mesh.point.coord)

        if ( numel(mesh.point.coord{idim}) ~= ...
            length(mesh.point.coord{idim}) )
            error('Incorrect dimensions!') ;
        end

        dims(iord(idim)) = ...
            length(mesh.point.coord{idim}) ;

        if (isa(mesh.point.coord{idim}, 'double'))
            vstr = sprintf('%%1.%ug\n',+16);
        else
            vstr = sprintf('%%1.%ug\n',+ 8);
        end

        fprintf(ffid,...
        'COORD=%u;%u\n', [idim,dims(iord(idim))]);

        fprintf(ffid,vstr,mesh.point.coord{idim});

    end

    end

    if (isfield(mesh,'value'))

    %------------------------------------ write "VALUE" data

        if (~isnumeric(mesh.value))
            error('Incorrect input types') ;
        end
        if (ndims(mesh.value) ~= length(dims)+0 && ...
            ndims(mesh.value) ~= length(dims)+1 )
            error('Incorrect dimensions!') ;
        end

        if (ndims(mesh.value) == length(dims))
            nval = size(mesh.value);
            nval = [nval, +1] ;
        else
            nval = size(mesh.value);
        end

        if (isvector(mesh.value))
        if (prod(nval(1:end-1)) ~= prod(dims))
            error('Incorrect dimensions!') ;
        end
        else
        if (~all(nval(1:end-1) == dims))
            error('Incorrect dimensions!') ;
        end
        end

        if     (isa(mesh.value, 'double'))
        vstr = sprintf('%%1.%ug;',+16) ;
        elseif (isa(mesh.value, 'single'))
        vstr = sprintf('%%1.%ug;',+ 8) ;
        elseif (isa(mesh.value,'integer'))
        vstr = '%d;' ;
        else
            error('Incorrect input type!');
        end
        vstr = repmat(vstr,+1,nval(end)) ;

        vals = ...
        reshape(mesh.value,[],nval(end)) ;

        fprintf(ffid, ...
          'VALUE=%u;%u\n',[prod(dims),nval(end)]);

        fprintf(ffid,[vstr(+1:end-1),'\n'],vals');

    end

    if (isfield(mesh,'slope'))

    %------------------------------------ write "SLOPE" data

        if (~isnumeric(mesh.slope))
            error('Incorrect input types') ;
        end
        if (ndims(mesh.slope) ~= length(dims)+0 && ...
            ndims(mesh.slope) ~= length(dims)+1 )
            error('Incorrect dimensions!') ;
        end

        if (ndims(mesh.slope) == length(dims))
            nval = size(mesh.slope);
            nval = [nval, +1] ;
        else
            nval = size(mesh.slope);
        end

        if (isvector(mesh.slope))
        if (prod(nval(1:end-1)) ~= prod(dims))
            error('Incorrect dimensions!') ;
        end
        else
        if (~all(nval(1:end-1) == dims))
            error('Incorrect dimensions!') ;
        end
        end

        if     (isa(mesh.slope, 'double'))
        vstr = sprintf('%%1.%ug;',+16) ;
        elseif (isa(mesh.slope, 'single'))
        vstr = sprintf('%%1.%ug;',+ 8) ;
        elseif (isa(mesh.slope,'integer'))
        vstr = '%d;' ;
        else
            error('Incorrect input type!');
        end
        vstr = repmat(vstr,+1,nval(end)) ;

        vals = ...
        reshape(mesh.slope,[],nval(end)) ;

        fprintf(ffid, ...
          'SLOPE=%u;%u\n',[prod(dims),nval(end)]);

        fprintf(ffid,[vstr(+1:end-1),'\n'],vals');

    end

end



