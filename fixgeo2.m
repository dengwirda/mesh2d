function [node,PSLG,part] = fixgeo2(varargin)
%FIXGEO2 attempts to "fix" issues with geometry definitions.
%   [NNEW,ENEW,PNEW] = FIXGEO2(NODE,EDGE,PART) returns a new
%   "repaired" geometry definition. Currently, the following
%   operations are performed:
%
%   (1) redundant nodes are "zipped" together.
%   (2) redundant edges are deleted.
%   (3) edges are split about intersecting nodes.
%   (4) ...
%
%   See also REFINE2

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 12/06/2017
%-----------------------------------------------------------
    
    filename = mfilename('fullpath') ;
    filepath = fileparts( filename ) ;
    
    addpath([filepath,'/aabb-tree']) ;

%---------------------------------------------- extract ARGS
    node = []; PSLG = []; part = {};

    if (nargin>=+1), node = varargin{1}; end
    if (nargin>=+2), PSLG = varargin{2}; end
    if (nargin>=+3), part = varargin{3}; end

%---------------------------------------------- default EDGE
    nnum = size(node,1);
    
    if (isempty(PSLG))
        PSLG = [(1:nnum-1)',(2:nnum-0)';nnum,1];
    end
    
%---------------------------------------------- default PART    
    enum = size(PSLG,1);
    
    if (isempty(part)), part{1} = (1:enum)'; end
    
%---------------------------------------------- quick return    
    if (isempty(node)), return; end
    if (isempty(PSLG)), return; end
    if (isempty(part)), return; end
  
%---------------------------------------------- basic checks    
    if ( ~isnumeric(node) || ...
         ~isnumeric(PSLG) || ...
         ~iscell(part) )
        error('fixgeo2:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(node) ~= +2 || ...
        ndims(PSLG) ~= +2 )
        error('fixgeo2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;
    end
    if (size(node,2)~= +2 || ...
        size(PSLG,2)~= +2 )
        error('fixgeo2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;
    end
    
%---------------------------------------------- basic checks
    if (min([PSLG(:)])<+1 || ...
            max([PSLG(:)]) > nnum)
        error('fixgeo2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
    
    pmin = cellfun(@min,part);
    pmax = cellfun(@max,part);
    
    if (min([pmin(:)])<+1 || ...
            max([pmax(:)]) > enum)
        error('fixgeo2:invalidInputs', ...
            'Invalid PART input array.') ;
    end
    
%------------------------------------- prune redundant nodes
   [node,PSLG,part] = ...
        prunenode(node,PSLG,part) ;
   
%------------------------------------- prune redundant edges
   [node,PSLG,part] = ...
        pruneedge(node,PSLG,part) ;
   
%------------------------------------- node-edge & edge-edge
   [node,PSLG,part] = ...
        intersect(node,PSLG,part) ;

    if (size(node,1) ~= nnum || ...
        size(PSLG,1) ~= enum ) 
    
%------------------------------------- iterate if any change
   [node,PSLG,part] = ...
          fixgeo2(node,PSLG,part) ;
    
    end
        
end

function [node,PSLG,part] = prunenode(node,PSLG,part)
%PRUNENODE "prune" redundant nodes by "zipping" those within
%tolerance of each other.

%------------------------------------- calc. "zip" tolerance
    nmin = min(node,[],+1) ;
    nmax = max(node,[],+1) ;  
    ndel = nmax - nmin;
    
    ztol = eps ^ 0.80 ;
    
    zlen = ztol * max(ndel);
    
%------------------------------------- index clustered nodes 
    ball = zeros(size(node,1),3);
    ball(:,1:2) = node(:,1:2);
    ball(:,  3) = zlen * zlen;

   [vp,vi] = ...
      findball(ball,node(:,1:2));

%------------------------------------- "zip" clustered nodes
   [vt,iv] = ...
      sort(vp(:,2) - vp(:,1));
    
    izip = zeros(size(node,1),1);
    imap = zeros(size(node,1),1);
    
    for kk = size(vp,1):-1:+1
        ii = iv(kk);
        for ip = vp(ii,1) ...
               : vp(ii,2)
            jj = vi(ip) ;
            if (izip(ii) == 0 && ...
                izip(jj) == 0 && ...
                ii~= jj  )
            
                izip(jj) = ii ;
            
            end 
        end
    end

%------------------------------------- re-index nodes//edges
    next = +1 ;
    for kk = +1:+1:size(vp,1)
        if (izip(kk) == +0)
            imap(kk) = next ;
            next = next + 1 ;
        end
    end
    
    imap(izip ~= 0) = ...
       imap(izip(izip ~= 0));
    
    PSLG = imap(PSLG) ;

    node = node(izip == 0,:);

end

function [node,PSLG,part] = pruneedge(node,PSLG,part)
%PRUNEEDGE "prune" redundant edges.

%------------------------------------- prune redundant topo. 
   [ptmp,ivec,jvec] = ...
        unique(sort(PSLG,+2),'rows') ;

    PSLG = PSLG(ivec,:);
    
    for ppos = +1:length(part)
    
        part{ppos} = ...
            unique(jvec(part{ppos})) ;
    
    end
    
%------------------------------------- prune collapsed topo.
    keep = diff(PSLG,[],2) ~= +0 ;
    
    jvec = zeros(size(PSLG,1),1) ;
    jvec(keep) = +1;
    jvec = cumsum(jvec);

    PSLG = PSLG(keep,:);
    
    for ppos = +1:length(part)
    
        part{ppos} = ...
            unique(jvec(part{ppos})) ;
    
    end

end
    
function [node,PSLG,part] = intersect(node,PSLG,part)
%INTERSECT split PSLG about intersecting nodes and/or edges.

%------------------------------------- node//edge intersect!
   [lp,li,TR] = findline(...
    node(PSLG(:,1),1:2), ...
    node(PSLG(:,2),1:2),node(:,1:2)) ;
    
    seen = false(size(PSLG,1),1); 
    ediv = zeros(size(PSLG,1),1);

    next = size(PSLG,1) ;
    
    PSLG = [PSLG; zeros(size(PSLG))] ;
    
    for ii = +1:+1:size(lp,1)
        for ip = lp(ii,1) ...
               : lp(ii,2)
            jj = li(ip) ;
            ni = PSLG(jj,1) ;
            nj = PSLG(jj,2) ;
            if (ni ~= ii && ...
                nj ~= ii && ...
                ~seen(jj) )
                
            PSLG(  jj,1) = ni ;
            PSLG(  jj,2) = ii ;

            next = next + 1;

            PSLG(next,1) = ii ;
            PSLG(next,2) = nj ;

            ediv(jj) = next;             
            seen(jj) = true;
                
            end
        end
    end
    PSLG = PSLG(1:next,:);
    
%------------------------------------- re-index edge in part
    for ppos = +1:length(part)
    
        enew = ediv(part{ppos});
        enew = enew(enew ~= 0) ;
        
        part{ppos} = ...
            [part{ppos}; enew] ;
    
    end

%------------------------------------- edge//edge intersect!
    
    %%!!todo: ...
    
    
end


