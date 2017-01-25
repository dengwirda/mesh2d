function [tp,tj,tr] = findtria(pp,tt,pj,ty,varargin)
%FINDTRIA spatial queries for collections of d-simplexes.
%   [TP,TI] = FINDTRIA(PP,TT,PI,TY) finds the set of simple-
%   xes that intersect with a given spatial query. Simplexes 
%   are specified via the vertex array PP = [X1,X2,...,XN] 
%   and the indexing array TT = [T1,T2,...,TM], such that 
%   the vertex positions for the II-th simplex are the poin-
%   ts [PP(TT(II,1),:),PP(TT(II,2),:),...,PP(TT(II,M),:)].
%
%   Simplexes are NOT required to form a conforming triangu-
%   lation. Specifically, non-delaunay, non-convex and even 
%   overlapping configurations are supported. Multiple matc-
%   hes may be returned if the collection is overlapping.
%
%   A set of intersecting simplexes is returned for each 
%   query point in PI, such that the II-th point is associa-
%   ted with the simplexes TI(TP(II,1):TP(II,2)). Unenclosed 
%   points have TP(II,1)==+0 and TP(II,2)==+0. Currently, 
%   the following spatial queries are supported:
%
%   * Point-location: TY = 'PTS', PI = [X1,X2,...,XN] is an
%     array of coordinates, and the set [TP,TI] represents
%     the enclosing simplxes for each point in PI.
%   * ...
%
%   Additional query types may be added in the future.
%
%   In general, query points may be matched to multiple sim-
%   plexes, but in cases when single matches are guaranteed, 
%   or if only a single match is desired, the following ret-
%   urns a singly-matched indexing array that is consistent 
%   with MATLAB's existing point-location routines:
%
%      [tp,tj] = findtria(pp,tt,pi,'pts') ;
%       ti = nan(size(tp,1),1);
%       in = tp(:,1) > +0;
%       ti(in) = tj(tp(in,+1));
%
%   [TP,TI,TR] = FINDTRIA(PP,TT,PI) additionally returns the
%   supporting aabb-tree used internally to compute the que-
%   ry. If the underlying collection [PP,TT] is static, the 
%   tree TR may be passed to subsequent calls, via [...] = 
%   FINDTRIA(PP,TT,PI,TY,TR). This syntax may lead to impro-
%   ved performance, especially when the number of simplexes 
%   is large w.r.t. the number of query points. Note that in 
%   such cases the underlying simplexes are NOT permitted to 
%   change between calls, or erroneous results may be retur-
%   ned. Additional parameters used to control the creation 
%   of the underlying aabb-tree may also be passed via [...] 
%   = FINDTRIA(PP,TT,PI,TY,TR,OP). See MAKETREE for additio-
%   nal information.
%
% See also FINDTRIADEMO, MAKETREE

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 25/01/2017

    tp = []; tj = []; tr = []; op = [];
%------------------------------ quick return on empty inputs
    if (isempty(pj)), return; end
%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ...
        ~isnumeric(tt) || ...
        ~isnumeric(pj) || ...
        ~ischar(ty)  )
        error('Incorrect input class.') ;
    end
%---------------------------------------------- basic checks
    if (nargin < +3 || nargin > +6)
        error('Incorrect number of inputs.');
    end
    if (ndims(pp) ~= +2 || size(pp,2) < +2 || size(pp,2) > size(tt,2))
        error('Incorrect input dimensions.');
    end
    if (ndims(tt) ~= +2 || size(tt,2) < +3)
        error('Incorrect input dimensions.');
    end
    
    switch (lower(ty))
        case 'pts'
    %----------------------------- check query array for PTS
    if (ndims(pj) ~= +2 || size(pj,2) ~= size(pp,2)*+1 )
        error('Incorrect input dimensions.');
    end    
    
        case 'ray'
    %----------------------------- check query array for RAY
    if (ndims(pj) ~= +2 || size(pj,2) ~= size(pp,2)*+2 )
        error('Incorrect input dimensions.');
    end
    
        otherwise
    %----------------------------- query must be unsupported
        error('Unsupported query type.');
    end
%------------------------------- extract user-defined inputs
    if (nargin >= +5), tr = varargin{1}; end
    if (nargin >= +6), op = varargin{2}; end
%---------------------------------------------- basic checks
    if (~isempty(tr) && ~isstruct(tr) )
        error('Incorrect input class.') ;
    end
    if (~isempty(op) && ~isstruct(op) )
        error('Incorrect input class.') ;
    end
    
    if (isempty(tr))
%------------------------------ compute aabb's for triangles
    bi = pp(tt(:,1),:); bj = pp(tt(:,1),:);
    for ii = 2 : size(tt,2)    
        bi = min(bi,pp(tt(:,ii),:)) ;
        bj = max(bj,pp(tt(:,ii),:)) ;
    end
    bb = [bi,bj];
%------------------------------ compute aabb-tree for aabb's
    if (isempty(op))                % scale tree with |pj|
        op.nobj = ceil(size(tt,1)/...
                       size(pj,1)) * +4 ;                   
        op.nobj = max( +32,op.nobj);  
        op.nobj = min(+256,op.nobj);
    end
    tr = maketree(bb,op);           % compute aabb-tree
    else
%---------------------------------- check existing aabb-tree
    if (~isfield(tr,'xx') || ...
        ~isfield(tr,'ii') || ...
        ~isfield(tr,'ll') )
        error('Incorrect aabb-tree.') ;
    end
    end
    
%---------------------------------- create tree-item mapping
    switch (lower(ty))
        case 'pts'
    tm = pushvert(tr,pj);
        case 'ray'
%%!!todo:
%%!!tm = findray (tr,pj);
    end
    
%------- compute spatial query using tree to induce a tiling
    ic = cell(size(tm.ii,1),1);
    jc = cell(size(tm.ii,1),1);
    for ip = 1 : size(tm.ii,1)
    %-------------------------- extract trias/items per tile
        ni = tm.ii(ip,1) ;          % node (in tree)
    %-------------------------- do O(n*m) search within tile
        switch (lower(ty))
            case 'pts'
        %---------------------- calc. pts-tria intersections
        [pi,ti] = testpts(pp,tt,pj,tm.ll{ip,1},tr.ll{ni,1});
                
            case 'ray'
        %---------------------- calc. ray-tria intersections
            %%!! todo:   
        
        end     
        ic{ip} = pi;
        jc{ip} = ti;        
    end
  
%-------------------------------- concat matches into arrays
    ii = vertcat(ic{:});
    tj = vertcat(jc{:});
    nj = length(tj);
%-------------------------------- form sparse-style indexing
   [ii,ix] = sort (ii) ; 
    tj = tj(ix);
    ix = find(diff(ii) > +0);
%------------------------- IN = TRUE if we found any matches
    in = false(size(pj,1),1);
    in(ii) = true;
%------------------------- ptrs into TJ for each point in PJ
    tp = zeros(size(pj,1),2);
    tp(in,1) = [+1; ix+1];
    tp(in,2) = [ix; nj+0];
    
end

function [ip,it] = testpts(pp,tt,pi,pk,tk)
%TESTPTS compute the "point-tria" matches within a tile.

    mp = length(pk); 
    mt = length(tk);

    switch (size(tt,2))
        case 3
    %------------------------------------ pts in 2-simplexes
        pk = pk';
        pk = pk(ones(mt,1),:); pk = pk(:);
        tk = tk(:,ones(1,mp)); tk = tk(:);

        in = intria2(pp,tt(tk,:),pi(pk,:));
        ip = pk(in);
        it = tk(in);

        case 4
    %------------------------------------ pts in 3-simplexes
        pk = pk';
        pk = pk(ones(mt,1),:); pk = pk(:);
        tk = tk(:,ones(1,mp)); tk = tk(:);

        in = intria3(pp,tt(tk,:),pi(pk,:));
        ip = pk(in);
        it = tk(in);

        otherwise
    %------------------------------------ pts in d-simplexes
       [il,jl] = intrian(pp,tt(tk,:),pi(pk,:));
        ip = pk(il(:));
        it = tk(jl(:));
    end

end

function [in] = intria2(pp,tt,pi)
%INTRIA2 returns TRUE for points enclosed by 2-simplexes.

    vi = pp(tt(:,1),:) - pi;
    vj = pp(tt(:,2),:) - pi;
    vk = pp(tt(:,3),:) - pi;
%------------------------------- compute sub-volume about PI
    aa = zeros(size(tt,1),3) ;
    aa(:,1) = vi(:,1).*vj(:,2) - ...
              vj(:,1).*vi(:,2) ;
    aa(:,2) = vj(:,1).*vk(:,2) - ...
              vk(:,1).*vj(:,2) ;
    aa(:,3) = vk(:,1).*vi(:,2) - ...
              vi(:,1).*vk(:,2) ;         
%------------------------------- PI is internal if same sign
    rt = eps^.8 * vol2(pp,tt).^2 ;
    in = aa(:,1).*aa(:,2) >= -rt ...
       & aa(:,2).*aa(:,3) >= -rt ...
       & aa(:,3).*aa(:,1) >= -rt ;
        
end

function [in] = intria3(pp,tt,pi)
%INTRIA3 returns TRUE for points enclosed by 3-simplexes.

    v1 = pi - pp(tt(:,1),:);
    v2 = pi - pp(tt(:,2),:);
    v3 = pi - pp(tt(:,3),:);
    v4 = pi - pp(tt(:,4),:);
%------------------------------- compute sub-volume about PI
    aa = zeros(size(tt,1),4) ;
    aa(:,1) = ...
   +v1(:,1).*(v2(:,2).*v3(:,3) - ...
              v2(:,3).*v3(:,2) ) ...
   -v1(:,2).*(v2(:,1).*v3(:,3) - ...
              v2(:,3).*v3(:,1) ) ...
   +v1(:,3).*(v2(:,1).*v3(:,2) - ...
              v2(:,2).*v3(:,1) ) ;
    aa(:,2) = ...
   +v1(:,1).*(v4(:,2).*v2(:,3) - ...
              v4(:,3).*v2(:,2) ) ...
   -v1(:,2).*(v4(:,1).*v2(:,3) - ...
              v4(:,3).*v2(:,1) ) ...
   +v1(:,3).*(v4(:,1).*v2(:,2) - ...
              v4(:,2).*v2(:,1) ) ;
    aa(:,3) = ...
   +v2(:,1).*(v4(:,2).*v3(:,3) - ...
              v4(:,3).*v3(:,2) ) ...
   -v2(:,2).*(v4(:,1).*v3(:,3) - ...
              v4(:,3).*v3(:,1) ) ...
   +v2(:,3).*(v4(:,1).*v3(:,2) - ...
              v4(:,2).*v3(:,1) ) ;
    aa(:,4) = ...
   +v3(:,1).*(v4(:,2).*v1(:,3) - ...
              v4(:,3).*v1(:,2) ) ...
   -v3(:,2).*(v4(:,1).*v1(:,3) - ...
              v4(:,3).*v1(:,1) ) ...
   +v3(:,3).*(v4(:,1).*v1(:,2) - ...
              v4(:,2).*v1(:,1) ) ;
%------------------------------- PI is internal if same sign
    rt = eps^.8 * vol3(pp,tt).^2 ;
    in = aa(:,1).*aa(:,2) >= -rt ...
       & aa(:,1).*aa(:,3) >= -rt ...
       & aa(:,1).*aa(:,4) >= -rt ...
       & aa(:,2).*aa(:,3) >= -rt ...
       & aa(:,2).*aa(:,4) >= -rt ...
       & aa(:,3).*aa(:,4) >= -rt ;
       
end

function [vol2] = vol2(pp,t2)
%VOL2 returns (signed) "volumes" for 2-simplexes.

    vv12 = pp(t2(:,2),:)-pp(t2(:,1),:);
    vv13 = pp(t2(:,3),:)-pp(t2(:,1),:);

    switch (size(pp,2))
        case +2
           
        vol2 = vv12(:,1).*vv13(:,2) ...
             - vv12(:,2).*vv13(:,1) ;
        vol2 = 0.5 * vol2;
            
        case +3
            
        avec = cross(vv12,vv13);
        vol2 = sqrt(sum(avec.^2,2)) ;
        vol2 = 0.5 * vol2;    
        
    end
    
end

function [vol3] = vol3(pp,t3)
%VOL2 returns (signed) "volumes" for 3-simplexes.

    vv14 = pp(t3(:,4),:)-pp(t3(:,1),:);
    vv24 = pp(t3(:,4),:)-pp(t3(:,2),:);
    vv34 = pp(t3(:,4),:)-pp(t3(:,3),:);

    vdet = + vv14(:,1) .* ...
            (vv24(:,2).*vv34(:,3)  ...
           - vv24(:,3).*vv34(:,2)) ...
           - vv14(:,2) .* ...
            (vv24(:,1).*vv34(:,3)  ...
           - vv24(:,3).*vv34(:,1)) ...
           + vv14(:,3) .* ...
            (vv24(:,1).*vv34(:,2)  ...
           - vv24(:,2).*vv34(:,1)) ;

    vol3 = vdet / 6.0;

end

function [ii,jj] = intrian(pp,tt,pi)
%INTRIAN return a list of points and enclosing "n"-simplexes.

   [np,pd] = size(pi);
   [nt,td] = size(tt);
%---------------- coefficient matrices for barycentric coord.
    mm = zeros(pd,pd,nt);
    for id = +1 : pd
    for jd = +1 : pd
        mm(id,jd,:) = pp(tt(:,jd),id) - ...
                      pp(tt(:,td),id) ;
    end
    end
%---------------- solve linear systems for barycentric coord.
    xx = zeros(pd,np,nt);
    vp = zeros(pd,np,+1);
    for ti = +1 : nt
    %---------------------------------------- form rhs coeff.
        for id = +1 : pd
            vp(id,:) = ...
            pi(:,id) - pp(tt(ti,td),id) ;
        end
    %---------------------------- actually faster to call LU-
       [ll,uu] = lu(mm(:,:,ti));
    %---------------------------- and then forward/back solve
        xx(:,:,ti) = uu\(ll\vp);
    end    
%-------------------- PI is internal if coord. have same sign
    in = all(xx >= +.0-(eps^.9),1) & ...
         sum(xx,1) <= +1.+(eps^.9) ;
%-------------------- find lists of matching points/simplexes
   [ii,jj] = find(reshape(in,[np,nt])) ;
       
end



