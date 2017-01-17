function makemsh(name,mesh)
%MAKEMSH make a *.MSH file for JIGSAW.
%
%   MAKEMSH(NAME,MESH);
%
%   The following entities are optionally written to "NAME.MSH". Ent-
%   ities are written if they are present in the sructure MESH:
%
%   MESH.POINT.COORD - [NPxND] array of point coordinates, where ND 
%       is the number of spatial dimenions.
%
%   MESH.EDGE2.INDEX - [N2x 3] array of indexing for edge-2 elements, 
%       where INDEX(K,1:2) is an array of "points" associated with 
%       the K-TH edge, and INDEX(K,3) is an ID tag for the K-TH edge.
%
%   MESH.TRIA3.INDEX - [N3x 4] array of indexing for tria-3 elements, 
%       where INDEX(K,1:3) is an array of "points" associated with 
%       the K-TH tria, and INDEX(K,4) is an ID tag for the K-TH tria.
%
%   MESH.QUAD4.INDEX - [N4x 5] array of indexing for quad-4 elements, 
%       where INDEX(K,1:4) is an array of "points" associated with 
%       the K-TH quad, and INDEX(K,5) is an ID tag for the K-TH quad.
%
%   MESH.TRIA4.INDEX - [M4x 5] array of indexing for tria-4 elements, 
%       where INDEX(K,1:4) is an array of "points" associated with 
%       the K-TH tria, and INDEX(K,5) is an ID tag for the K-TH tria.
%
%   MESH.HEXA8.INDEX - [M8x 9] array of indexing for hexa-8 elements, 
%       where INDEX(K,1:8) is an array of "points" associated with 
%       the K-TH hexa, and INDEX(K,9) is an ID tag for the K-TH hexa.
%
%   MESH.WEDG6.INDEX - [M6x 7] array of indexing for wedg-6 elements, 
%       where INDEX(K,1:6) is an array of "points" associated with 
%       the K-TH  wedg, and INDEX(K,7) is an ID tag for the K-TH wedg.
%
%   MESH.PYRA5.INDEX - [M5x 6] array of indexing for pyra-5 elements, 
%       where INDEX(K,1:5) is an array of "points" associated with 
%       the K-TH pyra, and INDEX(K,6) is an ID tag for the K-TH pyra.
%
%   See also READMSH, MAKEVTK, READVTK, MAKEMESH, READMESH, MAKEOFF, 
%            READOFF, MAKESTL, READSTL
%

%---------------------------------------------------------------------
%   Darren Engwirda
%   github.com/dengwirda/jigsaw-matlab
%   22-Mar-2016
%   d_engwirda@outlook.com
%---------------------------------------------------------------------
%

    if (~ischar  (name))
        error('NAME must be a valid file-name!') ;
    end
    if (~isstruct(mesh))
        error('MESH must be a valid structure!') ;
    end

   [path,file,fext] = fileparts(name);
   
    if(~strcmp(lower(fext),'.msh'))
        name = [name,'.msh'];
    end
 
    try
%-- try to write data to file
    
    ffid = fopen(name, 'w') ;
    
    nver = +1;
    
    fprintf(ffid,['# %s.msh file, created by JIGSAW','\n'],file);
    fprintf(ffid,['mshid=%u','\n'],nver);
    
    if (meshhas(mesh,'point'))
    
%-- write "POINT" data
        
    ndim = size(mesh.point.coord,2) - 1 ;    
    fprintf(ffid,['ndims=%u','\n'],ndim);
        
    fprintf(ffid,['point=%u','\n'],size(mesh.point.coord,1));        
    fprintf(ffid,[repmat('%1.16g;',1,ndim),'%i','\n'],...
        mesh.point.coord' ) ;

    end
    
    if (meshhas(mesh,'edge2'))
        
%-- write "EDGE2" data
        
    index = mesh.edge2.index;
    index(:,1:2) = index(:,1:2)-1 ; % file is zero-indexed!
       
    fprintf(ffid,['edge2=%u','\n'],size(mesh.edge2.index,1));        
    fprintf(ffid,[repmat('%u;',1,2),'%i','\n'],index');

    end
    
    if (meshhas(mesh,'tria3'))
       
%-- write "TRIA3" data
        
    index = mesh.tria3.index;
    index(:,1:3) = index(:,1:3)-1 ; % file is zero-indexed!
       
    fprintf(ffid,['tria3=%u','\n'],size(mesh.tria3.index,1));        
    fprintf(ffid,[repmat('%u;',1,3),'%i','\n'],index');

    end
    
    if (meshhas(mesh,'quad4'))
       
%-- write "QUAD4" data
        
    index = mesh.quad4.index;
    index(:,1:4) = index(:,1:4)-1 ; % file is zero-indexed!
       
    fprintf(ffid,['quad4=%u','\n'],size(mesh.quad4.index,1));        
    fprintf(ffid,[repmat('%u;',1,4),'%i','\n'],index');

    end
    
    if (meshhas(mesh,'tria4'))

%-- write "TRIA4" data
        
    index = mesh.tria4.index;
    index(:,1:4) = index(:,1:4)-1 ; % file is zero-indexed!
       
    fprintf(ffid,['tria4=%u','\n'],size(mesh.tria4.index,1));        
    fprintf(ffid,[repmat('%u;',1,4),'%i','\n'],index');

    end
    
    if (meshhas(mesh,'hexa8'))
    
%-- write "HEXA8" data
        
    index = mesh.hexa8.index;
    index(:,1:8) = index(:,1:8)-1 ; % file is zero-indexed!
       
    fprintf(ffid,['hexa8=%u','\n'],size(mesh.hexa8.index,1));        
    fprintf(ffid,[repmat('%u;',1,8),'%i','\n'],index');

    end
    
    if (meshhas(mesh,'wedg6'))
    
%-- write "WEDG6" data
        
    index = mesh.wedg6.index;
    index(:,1:6) = index(:,1:6)-1 ; % file is zero-indexed!
       
    fprintf(ffid,['wedg6=%u','\n'],size(mesh.wedg6.index,1));        
    fprintf(ffid,[repmat('%u;',1,6),'%i','\n'],index');

    end
    
    if (meshhas(mesh,'pyra5'))
    
%-- write "PYRA5" data
        
    index = mesh.pyra5.index;
    index(:,1:5) = index(:,1:5)-1 ; % file is zero-indexed!
       
    fprintf(ffid,['pyra5=%u','\n'],size(mesh.pyra5.index,1));        
    fprintf(ffid,[repmat('%u;',1,6),'%i','\n'],index');

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


