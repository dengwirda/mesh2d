function [vert,edge,tria] = triread(name)
%TRIREAD read two-dimensional triangulation data from file.
%   [VERT,EDGE,TRIA] = TRIREAD(NAME) returns the vertex, ed-
%   ge and triangle sets stored in the mesh-file NAME. Data
%   is returned non-empty if it is present in the file. 

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 17/01/2014

    addpath('mesh-file');

    vert = []; edge = []; tria = [];

%---------------------------------------------- basic checks
    if (~ischar(name))
        error('triread:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end

%----------------------------------- borrow JIGSAW I/O func!
   [mesh] = readmsh(name) ;

%----------------------------------- extract data if present
    if (meshhas(mesh,'point'))
    vert = mesh.point.coord(:,1:2) ;    
    end
    if (meshhas(mesh,'edge2'))
    edge = mesh.edge2.index(:,1:2) ;
    end
    if (meshhas(mesh,'tria3'))
    tria = mesh.tria3.index(:,1:3) ;
    end

end



