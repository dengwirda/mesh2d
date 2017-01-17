function [mesh] = readmsh(name)
%READMSH read a *.MSH file for JIGSAW.
%
%   MESH = READMSH(NAME);
%
%   The following entities are optionally read from "NAME.MSH". Ent-
%   ities are loaded if they are present in the file:
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
%   An additional set of "values" may also be optionally loaded, such 
%   that:
%
%       MESH.(BASE).(NAME).VALUE - [NFxNV] array of "values".
%       MESH.(BASE).(NAME).INDEX - [NFx 1] array of indexing.
%
%   Here BASE is a mesh entity: "POINT", "EDGE2", "TRIA3", "QUAD4", 
%   "TRIA4", "HEXA8", "WEDG6", "PYRA5", and NAME is the name assigned 
%   to the data field. This mechanism can be used to associate arbit-
%   rary "named" data fields with the primary mesh elements, where 
%   the data contained in MESH.(BASE).(NAME).VALUE(K,:) is associated 
%   with the MESH.(BASE).(NAME).INDEX(K)-TH element of MESH.(BASE).
%
%   See also MAKEMSH, MAKEVTK, READVTK, MAKEMESH, READMESH, MAKEOFF,
%            READOFF, MAKESTL, READSTL
%

%---------------------------------------------------------------------
%   Darren Engwirda
%   github.com/dengwirda/jigsaw-matlab
%   22-Mar-2016
%   d_engwirda@outlook.com
%---------------------------------------------------------------------
%

    try

    ffid = fopen(name,'r');
    
    real = '%f;';
    ints = '%i;';
    
    nver = +0;
    ndim = +0;
    
    while (true)
  
    %-- read next line from file
        lstr = fgetl(ffid);
        
        if (ischar(lstr) )
        
        if (length(lstr) > +0 && lstr(1) ~= '#')

        %-- tokenise line about '=' character
            tstr = regexp(lower(lstr),'=','split');
           
            switch (strtrim(tstr{1}))
            case 'mshid'

        %-- read "MSHID" data
        
                nver = str2double(tstr{2}) ;

            case 'ndims'

        %-- read "NDIMS" data
        
                ndim = str2double(tstr{2}) ;

            case 'point'

        %-- read "POINT" data

                nnum = str2double(tstr{2}) ;

                numr = nnum*(ndim+1);

                data = ...
            fscanf(ffid,[repmat(real,1,ndim),'%i'],numr);

                if (ndim == +2)
                mesh.point.coord = [ ...
                    data(1:3:end), ...
                    data(2:3:end), ...
                    data(3:3:end)] ;
                end
                if (ndim == +3)
                mesh.point.coord = [ ...
                    data(1:4:end), ...
                    data(2:4:end), ...
                    data(3:4:end), ...
                    data(4:4:end)] ;
                end

            case 'edge2'

        %-- read "EDGE2" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 3;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,2),'%i'],numr);
                
                mesh.edge2.index = [ ...
                    data(1:3:end), ...
                    data(2:3:end), ...
                    data(3:3:end)] ;
                    
                mesh.edge2.index(:,1:2) = ...
                mesh.edge2.index(:,1:2) + 1;
                
            case 'tria3'

        %-- read "TRIA3" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 4;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,3),'%i'],numr);
                
                mesh.tria3.index = [ ...
                    data(1:4:end), ...
                    data(2:4:end), ...
                    data(3:4:end), ...
                    data(4:4:end)] ;
                    
                mesh.tria3.index(:,1:3) = ...
                mesh.tria3.index(:,1:3) + 1;
            
            case 'quad4'

        %-- read "QUAD4" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 5;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,4),'%i'],numr);
                
                mesh.quad4.index = [ ...
                    data(1:5:end), ...
                    data(2:5:end), ...
                    data(3:5:end), ...
                    data(4:5:end), ...
                    data(5:5:end)] ;
        
                mesh.quad4.index(:,1:4) = ...
                mesh.quad4.index(:,1:4) + 1;
        
            case 'tria4'

        %-- read "TRIA4" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 5;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,4),'%i'],numr);
                
                mesh.tria4.index = [ ...
                    data(1:5:end), ...
                    data(2:5:end), ...
                    data(3:5:end), ...
                    data(4:5:end), ...
                    data(5:5:end)] ;
                    
                mesh.tria4.index(:,1:4) = ...
                mesh.tria4.index(:,1:4) + 1;
            
            case 'hexa8'

        %-- read "HEXA8" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 9;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,8),'%i'],numr);
                
                mesh.hexa8.index = [ ...
                    data(1:9:end), ...
                    data(2:9:end), ...
                    data(3:9:end), ...
                    data(4:9:end), ...
                    data(5:9:end), ...
                    data(6:9:end), ...
                    data(7:9:end), ...
                    data(8:9:end), ...
                    data(9:9:end)] ;
            
                mesh.hexa8.index(:,1:8) = ...
                mesh.hexa8.index(:,1:8) + 1;
            
            case 'wedg6'

        %-- read "WEDG6" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 7;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,6),'%i'],numr);
                
                mesh.wedg6.index = [ ...
                    data(1:7:end), ...
                    data(2:7:end), ...
                    data(3:7:end), ...
                    data(4:7:end), ...
                    data(5:7:end), ...
                    data(6:7:end), ...
                    data(7:7:end)] ;

                mesh.wedg6.index(:,1:6) = ...
                mesh.wedg6.index(:,1:6) + 1;
     
            case 'pyra5'

        %-- read "PYRA5" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 6;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,5),'%i'],numr);
                
                mesh.pyra5.index = [ ...
                    data(1:6:end), ...
                    data(2:6:end), ...
                    data(3:6:end), ...
                    data(4:6:end), ...
                    data(5:6:end), ...
                    data(6:6:end)] ;
            
                mesh.pyra5.index(:,1:5) = ...
                mesh.pyra5.index(:,1:5) + 1;
            
            case 'value'

        %-- read "VALUE" data

                stag = regexp(tstr{2},';','split');
                
                nnum = str2double(stag{1}) ;
                vnum = str2double(stag{2}) ;
                
                base = strtrim(stag{3}); 
                name = strtrim(stag{4});
                
                numr = nnum * (vnum+1) ;
                
                data = ...
            fscanf(ffid,[repmat(real,1,vnum),'%i'],numr);
                
                mesh.(base).(name).index = ...
                    data(vnum+1:vnum+1:end);
                
                mesh.(base).(name).value = ...
                    zeros(nnum, vnum);
                
                for vpos = +1 : vnum
                
                mesh.(base).(name).value(:,vpos) = ...
                    data(vpos+0:vnum+1:end);

                end
                
                mesh.(base).(name).index = ...
                mesh.(base).(name).index + 1 ;
                
            end
                       
        end
           
        else
    %-- if(~ischar(lstr)) //i.e. end-of-file
            break ;
        end
        
    end
    
    fclose(ffid) ;
    
    catch err

%-- ensure that we close the file regardless
    if (ffid>-1)
    fclose(ffid) ;
    end
    
    rethrow(err) ;
    
    end

end


