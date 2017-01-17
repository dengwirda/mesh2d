function [xx] = blocksolve(aa,bb)
%BLOCKSOLVE solve a collection of linear equations A*x = b.
%   X = BLOCKSOLVE(A,B) solves the block-wise collection of 
%   linear equations A*X = B. Each block is solved independ-
%   ently, such that A(:,:,IB) * X(:,IB,:) = B(:,IB,:) for 
%   the IB-th block. All blocks must be the same size.
%
%   BLOCKSOLVE is often faster than repeated calls to SLASH
%   when the block-size is sufficiently small (i.e. BSIZ <= 
%   +15) and the number of blocks is sufficiently large.
%   BLOCKSOLVE converts to a sparse block-diagonal format
%   internally and relies on MATLAB's sparse banded solver.
%
%   Example:
%     bs = +5 ;               % block size
%     nb = +10000 ;           % num. block
%     aa = rand(bs,bs,nb);
%     bb = rand(bs,nb,+1);
%     
%     tic, x1 = blocksolve(aa,bb) ; toc
%     
%     tic
%     x2 = zeros(bs,nb,+1);
%     for ib = +1 : nb
%         x2(:,ib,:) = aa(:,:,ib) \ ...
%                      bb(:,ib,:) ;
%     end
%     toc
%
%   See also SLASH.

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 15/12/2014

%---------------------------------------------- basic checks    
    if (~isnumeric(aa) || ...
        ~isnumeric(bb) )
        error('blocksolve:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
%---------------------------------------------- basic checks
    if (ndims(aa) < +2 || ndims(aa) > +3 || ...
        ndims(bb) < +2 || ndims(bb) > +3 )
        error('blocksolve:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end
    if (size(aa,1) ~= size(bb,1))
        error('blocksolve:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end

%---------------------------------------- extract dimensions
   [br,bc,nb] = size(aa) ;

    nv = size(bb,+3);       %% No. RHS
    nc = min(br*nb,+4096);  %% process in chunks of NC eqn.'s
    mb = ceil(nc/br);       %% blocks per chunk(s)
    nf = floor(nb/mb);      %% No. "full" chunk(s)
    
    xx = zeros(bc,nb, nv); kb = +0 ;
    
    if (nf > +1)
%---------------------------------- indexing for full chunks
    ii = reshape(+1:br*mb,br,+1,mb);
    ii = ii(:,ones(+1,bc),:);
    jj = reshape(+1:bc*mb,+1,bc,mb);
    jj = jj(ones(+1,br),:,:);
   
%-------- assemble & solve "internal" block-diagonal systems
    for jb = +1 : nf       
%---------------------------------- coeff. for current chunk
    ak = aa(:,:,kb+1:kb+mb);
    bk = bb(:,kb+1:kb+mb,:);
    bk = reshape(bk,br*mb,[]);  
%---------------------------------- _soln. for current chunk
    xx(:,kb+1:kb+mb,:) = ...
        reshape(sparse(ii(:),jj(:),ak(:))\bk,[bc,mb,nv]);
%-------------------------------------- advance block offset
    kb = kb + mb;
    end
    end

%-------- assemble & solve "trailing" block-diagonal systems    
    mb = nb - kb;
%------------------------------ indexing for trailing blocks
    ii = reshape(+1:br*mb,br,+1,mb);
    ii = ii(:,ones(+1,bc),:);
    jj = reshape(+1:bc*mb,+1,bc,mb);
    jj = jj(ones(+1,br),:,:);
%---------------------------------- coeff. in trailing chunk
    ak = aa(:,:,kb+1:kb+mb);
    bk = bb(:,kb+1:kb+mb,:);
    bk = reshape(bk,br*mb,[]) ;    
%---------------------------------- soln. for trailing chunk
    xx(:,kb+1:kb+mb,:) = ...
        reshape(sparse(ii(:),jj(:),ak(:))\bk,[bc,mb,nv]);

end

