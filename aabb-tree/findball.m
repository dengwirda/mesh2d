function [bp,bj,tr] = findball(varargin)
%FINDBALL spatial queries for collections of d-balls.
%   [BP,BI] = FINDBALL(BB,PI)
%
%   A set of intersecting balls is returned for each query 
%   point in PI, such that the II-th point is associated with 
%   the balls BI(BP(II,1):BP(II,2)). Unenclosed points have 
%   BP(II,1) == 0 and BP(II,2) == 0.
%
%   [BP,BI,TR] = FINDBALL(BB,PI) additionally returns the
%   supporting aabb-tree used internally to compute the que-
%   ry. If the underlying collection BB is static, the tree 
%   TR may be passed to subsequent calls, via [BP,BI,TR] = 
%   FINDBALL(BB,PI,TR). This syntax may lead to improved pe-
%   rformance, especially when the number of balls is large 
%   w.r.t. the number of query points. Note that in such ca-
%   ses the distribution of underlying balls is NOT permitt-
%   ed to change between calls, or erroneous results may be 
%   returned. Additional parameters used to control the cre-
%   ation of the underlying aabb-tree may also be passed via 
%   [...] = FINDBALL(BB,PI,TR,OP). See MAKETREE for additio-
%   nal information.
%
%   See also MAKETREE

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 05/01/2017

    bp = []; bj = []; tr = []; op = [];

%---------------------------------------------- basic checks
    if (nargin < +2 || nargin > +4)
        error('findball:incorrectNumInputs', ...
            'Incorrect number of inputs.');
    end

    if (nargin >= +1), bb = varargin{1}; end
    if (nargin >= +2), pp = varargin{2}; end
    if (nargin >= +3), tr = varargin{3}; end
    if (nargin >= +4), op = varargin{4}; end

%------------------------------ quick return on empty inputs
    if (isempty(bb)), return; end
%---------------------------------------------- basic checks    
    if (~isnumeric(bb) || ~isnumeric(pp))
        error('findball:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    
    if (ndims(bb) ~= +2 || size(bb,2) < +3 )
        error('findball:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end
    if (ndims(pp) ~= +2 || ...
            size(bb,2) ~= size(pp,2)+1)
        error('findball:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end

    if (~isempty(tr) && ~isstruct(tr) )
        error('findball:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    if (~isempty(op) && ~isstruct(op) )
        error('findball:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
  
    nd = size(pp,2);    
    if (isempty(tr))
%------------------------------ compute aabb-tree for d-ball
    rs = sqrt(bb(:,nd+1));
    rs = rs(:,ones(1,nd));   
    ab =[bb(:,1:nd)-rs, ...         % compute aabb-tree
         bb(:,1:nd)+rs]; 
    tr = maketree(ab,op) ;
    else
%---------------------------------- check existing aabb-tree
    if (~isfield(tr,'xx') || ...
        ~isfield(tr,'ii') || ...
        ~isfield(tr,'ll') )
        error('findball:incorrectAABBstruct', ...
            'Incorrect aabb-tree.') ;
    end
    end

    tm = pushvert(tr,pp) ;          % tree-vert mapping

%------------------------------ spatial query over tree-node
    ic = cell(size(tm.ii,1),1);
    jc = cell(size(tm.ii,1),1);
    for ip = 1 : size(tm.ii,1)
    %-------------------------- extract balls/verts per tile
        ni = tm.ii(ip,1) ;          % node (in tree)
        
    %-------------------------- do O(n*m) search within tile        
        pk = tm.ll{ip,1}';          % verts in tile
        bk = tr.ll{ni,1} ;          % balls in tile
        
    %-------------------------- push ball/vert onto n*m tile
        mp = length(pk); 
        mb = length(bk);       
        
        pk = pk(ones(mb,1),:); pk = pk(:); 
        bk = bk(:,ones(1,mp)); bk = bk(:);
  
    %-------------------------- compute O(n*m) loc. distance     
        dd = ...
    sum((pp(pk,+1:nd)-bb(bk,+1:nd)).^2,2);
        
        in = dd<=bb(bk,nd+1) ;
        
    %-------------------------- push loc. ball=>vert matches
        ic{ip} = pk(in);
        jc{ip} = bk(in);
    end

%-------------------------------- concat matches into arrays
    ii = vertcat(ic{:});
    bj = vertcat(jc{:});
    nj = length(bj);
%-------------------------------- form sparse-style indexing
   [ii,ix] = sort (ii) ; 
    bj = bj(ix);
    ix = find(diff(ii) > +0);
%------------------------- IN = TRUE if we found any matches
    in = false(size(pp,1),1);
    in(ii) = true;
%------------------------- ptrs into BJ for each point in PP
    bp = zeros(size(pp,1),2);
    bp( :,1) =  +1;
    bp(in,1) = [+1; ix+1];
    bp(in,2) = [ix; nj+0];

end



