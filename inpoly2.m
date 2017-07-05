function [stat] = inpoly2(varargin)
%INPOLY2 compute "points-in-polygon" queries.  
%   [STAT] = INPOLY2(VERT,NODE,EDGE) returns the "inside/ou-
%   tside" status for a set of vertices VERT and a polygon 
%   {NODE,EDGE} embedded in a two-dimensional plane. General
%   non-convex and multiply-connected polygonal regions can 
%   be handled. VERT is an N-by-2 array of XY coordinates to 
%   be tested. STAT is an associated N-by-1 logical array,
%   with STAT(II) = TRUE if VERT(II,:) is an interior point.
%   The polygonal region is defined as a piecewise-straight-
%   line-graph, where NODE is an M-by-2 array of polygon ve-
%   rtices and EDGE is a P-by-2 array of edge indexing. Each
%   row in EDGE represents an edge of the polygon, such that
%   NODE(EDGE(KK,1),:) and NODE(EDGE(KK,2),:) are the coord-
%   inates of the endpoints of the KK-TH edge. If the argum-
%   ent EDGE is omitted it assumed that the vertices in NODE
%   are connected in ascending order.
%
%   See also INPOLYGON

%   This algorithm is based on a "crossing-number" test, co-
%   unting the number of times a line extending from each 
%   point past the right-most region of the polygon interse-
%   cts with the polygonal boundary. Points with odd counts 
%   are "inside". A simple implementation requires that each
%   edge intersection be checked for each point, leading to 
%   O(N*M) complexity...
%
%   This implementation seeks to improve these bounds:
%
% * Binning the query points by y-value and determining can-
%   didate edge intersection sets via bin overlap. This red-
%   uces overall algorithmic complexity.
%
% * Carefully checking points against the bounding-box asso-
%   ciated with each polygon edge. This minimises the number
%   of calls to the (relatively) expensive edge intersection 
%   test.

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 05/07/2017

%---------------------------------------------- extract args
    node = []; edge = []; vert = [];

    if (nargin>=+1), vert = varargin{1}; end
    if (nargin>=+2), node = varargin{2}; end
    if (nargin>=+3), edge = varargin{3}; end
    
%---------------------------------------------- default args
    nnod = size(node,1) ;
    nvrt = size(vert,1) ;
    
    if (isempty(edge))
        edge = [(1:nnod-1)',(2:nnod)'; nnod,1];
    end
    
%---------------------------------------------- basic checks    
    if ( ~isnumeric(node) || ...
         ~isnumeric(edge) || ...
         ~isnumeric(vert) )
        error('inpoly2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end

%---------------------------------------------- basic checks
    if (ndims(node) ~= +2 || ...
        ndims(edge) ~= +2 || ...
        ndims(vert) ~= +2 )
        error('inpoly2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(node,2)~= +2 || ... 
        size(edge,2)~= +2 || ...
        size(vert,2)~= +2 )
        error('inpoly2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
%---------------------------------------------- basic checks
    if (min([edge(:)]) < +1 || ...
            max([edge(:)]) > nnod)
        error('inpoly2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end

%-------------- flip to ensure the y-axis is the "long" axis
    ddxy = max(vert,[],1) - ...
           min(vert,[],1) ;
    
    if (ddxy(1) > ddxy(2))
        vert = vert(:,[2,1]) ;
        node = node(:,[2,1]) ;
    end
    
%----------------------------------- sort points via y-value
    swap = node(edge(:,2),2) ...
         < node(edge(:,1),2) ;
         
    edge(swap,[1,2]) = ...
        edge(swap,[2,1]) ;
     
%----------------------------------- setup "bins" on y-value
    ybin = ...
      unique(node (:,2)) ;
    
    nbin = length (ybin) ;

%----------------------------------- "bin" edges via y-value
   [junk,bmin] = histc( ...
        node(edge(:,1),2),ybin);
        
   [junk,bmax] = histc( ...
        node(edge(:,2),2),ybin);
    
%----------------------------------- "bin" points on y-value
   [junk,vmap] = ...
        histc(vert( :,+2),ybin);
        
   [vmap,vidx] = sort(vmap);
   
    vidx = vidx(vmap~=0) ;
    vmap = vmap(vmap~=0) ;
   
    iptr = find( diff(vmap)) ;
   
    nidx = length (vidx) ;
    
    ibin = vmap([iptr;nidx]) ;
   
    vptr = zeros(nbin,2) ;
    vptr(ibin,1) = [+1; iptr+1];
    vptr(ibin,2) = [iptr; nidx];
 
%----------------------------------- do crossing-number test
    isoctave = exist( ...
        'OCTAVE_VERSION','builtin') > +0 ;
        
    stat = false(nvrt,1);
        
    for epos = +1 : size(edge,1)
        
    %------------------------------- find intersecting boxes   
        imin = bmin(epos);
        imax = bmax(epos);
        
        while (vptr(imin,1)==0 && imin < imax)
            imin = imin+1;
        end
        
        while (vptr(imax,1)==0 && imax > imin)
            imax = imax-1;
        end
   
        head = vptr(imin,1) ;
        tail = vptr(imax,2) ;
        
        if (head==+0), continue ; end
   
        if (isoctave)
        
    %-- OCTAVE is *shockingly* bad at executing loops, so -- 
    %-- even though it involves far more operations! -- call
    %-- the vectorised version below.
            
        inod = edge(epos,1) ;
        jnod = edge(epos,2) ;

    %------------------------------- calc. edge bounding-box
        yone = node(inod,2) ;
        ytwo = node(jnod,2) ;
        xone = node(inod,1) ;
        xtwo = node(jnod,1) ;
        
        ydel = ytwo - yone;
        xdel = xtwo - xone;
     
    %------------------------------- calc. edge-intersection
        vset = ...
            vidx(head:tail) ;
        
        xpos = vert(vset,1) ;
        ypos = vert(vset,2) ;
    
        cros = ypos >=yone  & ...
               ypos < ytwo  & ...
        ydel* (xpos - xone)<= ...
        xdel* (ypos - yone) ;
   
        stat(vset(cros)) = ...
            ~stat(vset(cros)) ;
    
        else

    %-- MATLAB is actually pretty good at JIT-ing code these
    %-- days, so use the asymptotically faster version based
    %-- on the pre-computed ordering.
        
        inod = edge(epos,1) ;
        jnod = edge(epos,2) ;

    %------------------------------- calc. edge bounding-box
        yone = node(inod,2) ;
        ytwo = node(jnod,2) ;
        xone = node(inod,1) ;
        xtwo = node(jnod,1) ;
        
        ydel = ytwo - yone;
        xdel = xtwo - xone;

        xmin = min(xone,xtwo) ;        
          
    %------------------------------- calc. edge-intersection  
        for kpos = head:tail
        
            jpos = vidx(kpos,1) ;            
            ypos = vert(jpos,2) ;
            if (ypos >=yone && ...
                ypos < ytwo )
                xpos = vert(jpos,1) ;
                if (xpos >= xmin)
                    if ( ...
                    ydel* (xpos-xone)<= ...
                    xdel* (ypos-yone) )
                    
                    stat(jpos) = ...
                        ~stat(jpos) ;
                    end
                else
                    stat(jpos) = ...
                        ~stat(jpos) ;
                end
            end
            
        end
              
        end
            
    end
         
 end
 
 
 
 