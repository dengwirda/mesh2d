function [vert,conn,tria,tnum] = ...
          deltri2(vert,conn,node,PSLG,part)
%DELTRI2 compute a constrained 2-simplex Delaunay triangula-
%tion in the two-dimensional plane.
%   [VERT,CONN,TRIA,TNUM]=DELTRI2(VERT,CONN,NODE,PSLG,PART)
%   computes the Delaunay trianguation {VERT,TRIA}, the con-
%   straints CONN, and the "inside" status vector TNUM. VERT
%   is an V-by-2 array of XY coordinates to be triangulated,
%   TRIA is a T-by-3 array of vertex indexing, where each 
%   row defines a triangle, such that VERT(TRIA(II,1),:), 
%   VERT(TRIA(II,2),:) and VERT(TRIA(II,3),:) are the coord-
%   inates of the II-TH triangle. CONN is a C-by-2 array of
%   constraining edges, where each row defines an edge, as
%   per TRIA. The additional arguments NODE,PSLG and PART 
%   define a (mutliply-connected) polygonal region, where 
%   NODE is an N-by-2 array of vertices and PSLG is a P-by-2 
%   array of edges (a piecewise-straight-line-graph), where 
%   each row defines an edge as a pair of indices into NODE.
%   PART is a cell-array of polygonal "parts", where each
%   element PART{KK} is an array of edge indices defining a
%   polygonal region. PSLG(PART{KK},:) is the set of edges
%   in the KK-TH part. TNUM is a T-by-1 array of part index-
%   ing, such that TNUM(II) is the index of the part in whi-
%   ch the II-TH triangle resides.
%
%   See also DELAUNAYTRIANGULATION, DELAUNAYTRI, DELAUNAYN

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 14/01/2017

%---------------------------------------------- basic checks    
    if (~isnumeric(vert) || ~isnumeric(conn) || ...
        ~isnumeric(node) || ~isnumeric(PSLG) || ...
            ~iscell(part) )
        error('deltri2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end

    nvrt = size(vert,+1) ; nnod = size(node,+1) ;
    nedg = size(PSLG,+1) ;

%---------------------------------------------- basic checks
    if (ndims(vert) ~= +2 || ndims(conn) ~= +2)
        error('deltri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(vert,2)~= +2 || size(conn,2)~= +2)
        error('deltri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    if (ndims(node) ~= +2 || ndims(PSLG) ~= +2)
        error('deltri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(node,2)~= +2 || size(PSLG,2)~= +2)
        error('deltri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
%---------------------------------------------- basic checks
    if (min([conn(:)])<+1 || max([conn(:)])>nvrt)
        error('deltri2:invalidInputs', ...
            'Invalid CONN input array.') ;
    end
    if (min([PSLG(:)])<+1 || max([PSLG(:)])>nnod)
        error('deltri2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
    if (min([part{:}])<+1 || max([part{:}])>nedg)
        error('deltri2:invalidInputs', ...
            'Invalid PART input array.') ;
    end

    if (exist( ...
    'delaunayTriangulation','class') )
%------------------------------------ use class if available    
    dtri = ...
    delaunayTriangulation(vert,conn) ;
    vert = dtri.Points;
    conn = dtri.Constraints;
    tria = dtri.ConnectivityList;
    else
    if (exist('delaunayTri','class') )
%------------------------------------ use class if available
    dtri = delaunayTri   (vert,conn) ;
    vert = dtri.X;
    conn = dtri.Constraints;
    tria = dtri.Triangulation;
    else
    error('deltri2:unsupportedMATLAB', ...
    'Delaunay triangulation is not supported');
    end
    end

%------------------------------------ calc. "inside" status!    
    tnum = zeros(size(tria,+1),+1) ;
    
    tmid = vert(tria(:,1),:) ...
         + vert(tria(:,2),:) ...
         + vert(tria(:,3),:) ;
    tmid = tmid / +3.0;
    
    for ppos = 1 : length(part)

       [stat] = inpoly2( ...
        tmid,node,PSLG(part{ppos},:)) ;
 
        tnum(stat) = ppos ;
        
    end    

%------------------------------------ keep "interior" tria's  
    tria = tria(tnum>+0,:);
    tnum = tnum(tnum>+0,:);

end



