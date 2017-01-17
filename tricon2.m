function [ee,tt] = tricon2(varargin)
%TRICON2 edge-centred connectivity for a conforming 2-simpl-
%ex triangulation in the two-dimensional plane.
%   [EE,TT] = TRICON2(TT,CC) returns the edge-based adjacen-
%   cy for a mesh of 2-simlexes (triangles). EE = [V1,V2,T1,
%   T2,CE] is the set of unique 1-simplexes (edges) in the 
%   mesh TT. Each row of {V1,V2} defines an edge, each row
%   of {T1,T2} defines the two triangles adjacent to an edge 
%   and CE is a "constraint" flag. TT = [V1,V2,V3,E1,E2,E3],
%   is the set of unique 2-simplexes in the mesh, where
%   {E1,E2,E3} define the tria-to-edge mapping. Each row of 
%   {E1,E2,E3} are the indicies of the three edges that make 
%   up each triangle.

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 14/01/2014

%---------------------------------------------- extract args
    tt = []; cc = [];

    if (nargin>=1), tt = varargin{1}; end
    if (nargin>=2), cc = varargin{2}; end

%---------------------------------------------- basic checks
    if (~isnumeric(tt))
        error('tricon2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    if (ndims(tt) ~= +2 || size(tt,2) ~= +3)
        error('tricon2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (min(min(tt)) < +1 )
        error('tricon2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end
    
    if (~isempty(cc))
    if (~isnumeric(cc))
        error('tricon2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    if (ndims(cc) ~= +2 || size(cc,2) ~= +2)
        error('tricon2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    end

    nt = size(tt,1);
    nc = size(cc,1);

%------------------------------------------ non-unique edges
    ee = zeros(nt*3,2);
    ee((1:nt)+nt*0,:) = tt(:,[1,2]);
    ee((1:nt)+nt*1,:) = tt(:,[2,3]);
    ee((1:nt)+nt*2,:) = tt(:,[3,1]);
    
%------------------------------ unique edges and re-indexing
   [ee, iv, jv] = ...
        unique(sort(ee, 2), 'rows');
    
%------------------- tria-to-edge indexing: 3 edges per tria 
    tt = [tt, zeros(nt*1,3)] ;
    tt(:,4) = jv((1:nt)+nt*0);
    tt(:,5) = jv((1:nt)+nt*1);
    tt(:,6) = jv((1:nt)+nt*2);
    
%------------------- edge-to-tria indexing: 2 trias per edge
    ne = size(ee,1);
    ee = [ee, zeros(ne*1,3)] ;
    ep = +3 * ones (ne*1,1)  ;
    for ti = +1 : nt
        ei = tt(ti,4) ; 
        ee(ei,ep(ei)) = ti;
        ej = tt(ti,5) ;
        ee(ej,ep(ej)) = ti;        
        ek = tt(ti,6) ;
        ee(ek,ep(ek)) = ti;
    
        ep(ei) = ep(ei)+1 ;
        ep(ej) = ep(ej)+1 ;
        ep(ek) = ep(ek)+1 ;
    end
       
    if (isempty(cc)), return; end
    
%------------------------------------ find constrained edges
    is = ismember( ...
       ee(:,1:2),sort(cc,2),'rows');
%------------------------------------ mark constrained edges
    ee(is,5) = +1 ;
    
end



