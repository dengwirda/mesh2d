function [vert,conn,tria,tnum] = refine2(varargin)
%REFINE2 (Frontal)-Delaunay-refinement for two-dimensional,
%polygonal geometries.
%   [VERT,EDGE,TRIA,TNUM] = REFINE2(NODE,EDGE) returns a co-
%   nstrained Delaunay triangulation of the polygonal region
%   {NODE,EDGE}. NODE is an N-by-2 array of polygonal verti-
%   ces and EDGE is an E-by-2 array of edge indexing. Each
%   row in EDGE represents an edge of the polygon, such that
%   NODE(EDGE(JJ,1),:) and NODE(EDGE(JJ,2),:) are the coord-
%   inates of the endpoints of the JJ-TH edge. If the argum-
%   ent EDGE is omitted it assumed that the vertices in NODE
%   are connected in ascending order.
%
%   [...] = REFINE2(NODE,EDGE,PART) computes a triangulation
%   for a multiply-connected geometry. PART is a cell-array 
%   of polygonal "parts", where each element PART{KK} is an 
%   array of edge indices defining a given polygonal region. 
%   EDGE(PART{KK}, :) is the set of edges in the KK-TH part.
%
%   VERT is a V-by-2 array of XY coordinates in the triangu-
%   lation, EDGE is an array of constrained edges, TRIA is a
%   T-by-3 array of triangles, and TNUM is a T-by-1 array of
%   part indices. Each row of TRIA and EDGE define an eleme-
%   nt. VERT(TRIA(II,1),:), VERT(TRIA(II,2),:) and VERT(TRIA
%   (II,3),:) are the coordinates of the II-TH triangle. The
%   edges in EDGE are defined in a similar manner. NUM is an
%   array of part indexing, such that TNUM(II) is the index 
%   of the part in which the II-TH triangle resides.
%
%   [...] = REFINE2(..., OPTS) passes an additional options 
%   structure OPTS, containing various user-defined paramet-
%   ers, including:
%
% - OPTS.KIND = {'DELFRONT'}, 'DELAUNAY' -- the type of ref-
%   inement employed. The 'DELFRONT' algorithm is typically
%   slower, but produces higher quality output.
%
% - OPTS.RHO2 = {1.025} -- the maximum allowable radius-edge 
%   ratio. Refinement proceeds until all interior triangles
%   satisfy the radius-edge threshold. Smaller radius-edge
%   ratios lead to improved triangle shape, with RHO2=1 req-
%   uiring that all angles exceed 30 degrees. Setting RHO2<1 
%   may lead to non-convergence.
%
% - OPTS.SIZ1 = {1.333} -- the normalised rel.-length th-
%   reshold for edge-elements. Each exterior edge is refined 
%   until LL/HH<SIZ1, where LL is the edge-length, HH is the
%   edge-centred mesh-size value.
% 
% - OPTS.SIZ2 = {1.300} -- the normalised rel.-length th-
%   reshold for tria-elements. Each interior tria is refined
%   until RE/HH<SIZ2, where RE is an effective tria length, 
%   based on the circumradius, HH is the tria-centred mesh-
%   size value.
%
% - OPTS.DISP = { +10 } -- refinement verbosity. Set to INF
%   for quiet execution.
%
%   [...] = REFINE2(..., HFUN,ARGS) also passes an optional
%   mesh-size function argument. Setting HFUN = HMAX, where 
%   HMAX is a scalar value, imposes a constant size constra-
%   int over the full domain. HFUN can also be defined as a 
%   general function handle [HH] = HFUN(PP), where PP is an
%   N-by-2 array of XY coordinates and HH is the associated
%   vector of mesh-size values. User-defined HFUN must be
%   fully vectorised. Additional arguments {A1,A2,...AN} for 
%   HFUN can be passed as trailing parameters to REFINE2. In
%   such cases, HFUN must adopt a signature [HH] = HFUN(PP,
%   A1,A2,...,AN). HFUN must return positive values.
%
%   See also SMOOTH2, DRAWSCR, TRIDEMO

%   This routine implements a "multi-refinement" variant of
%   Delaunay-refinement type mesh-generation. Both standard
%   Delaunay-refinement and Frontal-Delaunay type algorithms
%   are available. The Frontal-Delaunay approach is a simpl-
%   ified version of the algorithm described in: D. Engwirda
%   "Locally optimal Delaunay-refinement and optimisation-
%   based mesh generation", Ph.D. Thesis, Univ. of Sydney, 
%   2014. This work is an extension of the "off-centre" type
%   methodology introduced in: H. Erten and A. Ungor, "Qual-
%   ity triangulations with locally optimal Steiner points",
%   SIAM Journal on Sci. Comp. 31(3) 2009, pp: 2103--2130.

%   A more advanced, and fully three-dimensional implementa-
%   tion is available as part of the JIGSAW pacakge. For de-
%   tails, see github.com/dengwirda/jigsaw-matlab.

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 23/01/2017
    
    filename = mfilename('fullpath') ;
    filepath = fileparts( filename ) ;
    
    addpath([filepath,'/aabb-tree']) ;

%---------------------------------------------- extract args
    node = []; PSLG = []; part = {}; opts = [] ; 
    hfun = []; harg = {};

    if (nargin>=+1), node = varargin{1}; end
    if (nargin>=+2), PSLG = varargin{2}; end
    if (nargin>=+3), part = varargin{3}; end
    if (nargin>=+4), opts = varargin{4}; end
    if (nargin>=+5), hfun = varargin{5}; end
    if (nargin>=+6), harg = varargin(6:end); end

   [opts] = makeopt(opts);

%---------------------------------------------- default EDGE
    nnod = size(node,1);
    
    if (isempty(PSLG))
        PSLG = [(1:nnod-1)',(2:nnod)'; nnod,1] ;
    end
    
%---------------------------------------------- default PART    
    ncon = size(PSLG,1);
    
    if (isempty(part))
        part{1} = (1:ncon)';
    end
    
%---------------------------------------------- basic checks    
    if (~isnumeric(node) || ~isnumeric(PSLG) || ...
        ~iscell(part) || ~isstruct(opts))
        error('refine2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(node) ~= +2 || ndims(PSLG) ~= +2)
        error('refine2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(node,2)~= +2 || size(PSLG,2)~= +2)
        error('refine2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
%---------------------------------------------- basic checks
    if (min([PSLG(:)])<+1 || max([PSLG(:)])>nnod)
        error('refine2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
    
    pmin = cellfun(@min,part);
    pmax = cellfun(@max,part);
    
    if (min([pmin(:)])<+1 || max([pmax(:)])>ncon)
        error('refine2:invalidInputs', ...
            'Invalid PART input array.') ;
    end

%---------------------------------------------- output title
    if (~isinf(opts.disp))
        fprintf(1,'\n') ;
        fprintf(1,' Refine triangulation...\n') ;
        fprintf(1,'\n') ;
        fprintf(1,[...
' -------------------------------------------------------\n', ...
'      |ITER.|          |CDT1(X)|          |CDT2(X)|     \n', ...
' -------------------------------------------------------\n', ...
             ] ) ;
    end

%-------------------------------- PASS 0: inflate box bounds
    vert = node ; tria = []; tnum = []; 
    conn = PSLG ; iter = +0;

    vmin = min(vert,[],+1) ;    % inflate bbox for stability
    vmax = max(vert,[],+1) ;
    
    vdel = vmax - 1.*vmin;
    vmin = vmin - .5*vdel;
    vmax = vmax + .5*vdel;

    vbox = [vmin(1), vmin(2)
            vmax(1), vmin(2)
            vmax(1), vmax(2)
            vmin(1), vmax(2)
           ] ;
    vert = [vert; vbox] ;

%-------------------------------- PASS 1: refine 1-simplexes
   [vert,conn,tria,tnum,iter] = ...
        cdtref1(vert,conn,tria,tnum, ...
            node,PSLG,part,opts,hfun,harg,iter);

%-------------------------------- PASS 2: refine 2-simplexes
   [vert,conn,tria,tnum,iter] = ...
        cdtref2(vert,conn,tria,tnum, ...
            node,PSLG,part,opts,hfun,harg,iter);
    
%-------------------------------- trim extra adjacency info.  
    tria = tria(:,1:3) ;
    
    if (~isinf(opts.disp)), fprintf(1,'\n'); end
    
end

function [vert,conn,tria,tnum,iter] = ...
            cdtref1(vert,conn,tria,tnum, ...
                node,PSLG,part,opts,hfun,harg,iter)
%CDTREF1 constrained Delaunay-refinement for 1-simplex elem-
%nts embedded in R^2.
%   [...] = CDTREF1(...) refines the set of 1-simplex eleme-
%   nts embedded in the triangulation until all constraints 
%   are satisfied. Specifically, edges are refined until all
%   local mesh-spacing and encroachment conditions are met.
%   Refinement proceeds according to either a Delaunay-refi-
%   nement or Frontal-Delaunay type approach, depending on
%   user-settings. In either case, new steiner vertices are
%   introduced to split "bad" edges - those that violate the
%   set of prescribed constraints. In the "-DR" type process
%   edges are split about their circumballs (midpoints). In
%   the "-FD" approach, new vertices are positioned such th-
%   at mesh-spacing constraints are satisfied in a "locally-
%   optimal" fashion.

    ntol = +1.60;

    while  (true)
    
        iter = iter + 1 ;
    
        if (iter>=opts.iter),break; end
    
    %------------------------------------- calc. circumballs
        bal1 = cdtbal1(vert,conn) ;
        
    %------------------------------------- eval. length-fun.
        if (~isempty(hfun))
            if (isnumeric(hfun))
            fun0 = hfun * ...
              ones(size(vert,1),1);
            fun1 = hfun ;
            else
            fun0 = feval( ...
                hfun,vert,harg{:});
            fun1 = fun0(conn(:,1))...
                 + fun0(conn(:,2));
            fun1 = fun1 / +2. ;
            end
        else
            fun0 = +inf * ...
              ones(size(vert,1),1);
            fun1 = +inf ;
        end
    
        siz1 = ...
         +4. * bal1(:,3)./(fun1.*fun1) ;
  
    %------------------------------------- test encroachment
        bal1(:,3) = ...
            (1.-eps^.75) * bal1(:,3) ;
  
       [vp,vi] = ...
          findball(bal1,vert(:,1:2)) ;

    %------------------------------------- near=>[vert,edge]
        ebad = false(size(conn,1),1) ;
        near = [];
        for ii = +1 : size(vp,1)
            for ip = vp(ii,1):vp(ii,2)
                jj = vi(ip) ;
                if (ii ~= conn(jj,1) ...
                &&  ii ~= conn(jj,2) )
                near = [near; ii,jj] ;
                end
            end
        end
        
        if (~isempty(near))
    %-- mark edge "encroached" if there is a vert within its
    %-- dia.-ball that is not joined to either of its vert's
    %-- via an existing edge...         
            ivrt = conn(near(:,2),1);
            jvrt = conn(near(:,2),2);
 
            conn = sort(conn,2);
        
            pair = [near(:,1), ivrt];   
            ivec = ismember( ...
            sort(pair,2),conn,'rows') ;
            
            pair = [near(:,1), jvrt];
            jvec = ismember( ...
            sort(pair,2),conn,'rows') ;
           
            okay = ~ivec & ~jvec ;
            
            ebad(near(okay,2))=true ;
      
        end
        
    %------------------------------------- refinement queues
        ref1 = false(size(conn,1),1);      
        ref1(ebad)           = true ;   %- edge encroachment
        ref1(siz1>opts.siz1* ...        %- bad equiv. length
                  opts.siz1) = true ;
        
        num1 = find(ref1)  ;
      
    %------------------------------------- dump-out progess!
        if (mod(iter,opts.disp)==+0 || ...
                isempty(num1) )
            numc = size(conn,1) ;
            numt = size(tria,1) ;
            fprintf(+1, ...
            '%11i %18i %18i\n', ...
            [iter,numc,numt]) ;
        end
      
    %------------------------------------- nothing to refine
        if (isempty(num1)), break; end
        
    %------------------------------------- refine "bad" tria
        switch (lower(opts.kind))
        case 'delaunay'
    %------------------------------------- do circ-ball pt's
        new1 = bal1(ref1, 1:2) ;
        
        vnew = (1:size(new1,1))' ...
                + size(vert,1) ;
        
        cnew = [conn( ref1,1), vnew
                conn( ref1,2), vnew];
        conn = [conn(~ref1,:); cnew];
            
    %------------------------------------- update vertex set    
        vert = [vert; new1(:,1:2)];
  
        
        case 'delfront'
    %-- symmetric off-centre scheme:- refine edges from both
    %-- ends simultaneously, placing new vertices to satisfy
    %-- the worst of mesh-spacing and local voronoi constra-
    %-- ints.
  
        evec = vert(conn(ref1,2),:) ...
             - vert(conn(ref1,1),:) ;
        elen = sqrt(sum(evec.^2,2)) ;
        evec = evec ./ [elen, elen] ;
  
    %------------------------------------- "voro"-type dist.
        olen = sqrt(bal1(ref1,3));
        
    %------------------------------------- "size"-type dist.
        ihfn = fun0(conn(ref1,1));
        jhfn = fun0(conn(ref1,2));
        
    %------------------------------------- bind "safe" dist.
        ilen = min(olen,ihfn) ;
        jlen = min(olen,jhfn) ;
 
    %------------------------------------- locate offcentres       
        inew = vert(conn(ref1,1),:) ...
             + [ilen,ilen].*evec ;
        jnew = vert(conn(ref1,2),:) ...
             - [jlen,jlen].*evec ;
             
    %------------------------------------- iter. "size"-type
        for ioff = +1 : +4
    %------------------------------------- eval. length-fun.
        if (~isempty(hfun))
            if (isnumeric(hfun))
            iprj = hfun * ...
              ones(size(inew,1),1);
            jprj = hfun * ...
              ones(size(jnew,1),1);
            else
            iprj = feval( ...
                hfun,inew,harg{:});
            jprj = feval( ...
                hfun,jnew,harg{:});
            end
        else
            iprj = +inf * ...
              ones(size(inew,1),1);
            jprj = +inf * ...
              ones(size(jnew,1),1);
        end
        
        iprj = 0.5*ihfn + 0.5*iprj;
        jprj = 0.5*jhfn + 0.5*jprj;

    %------------------------------------- bind "safe" dist.
        ilen = min(olen,ihfn) ;
        jlen = min(olen,jhfn) ;
 
    %------------------------------------- locate offcentres       
        inew = vert(conn(ref1,1),:) ...
             + [ilen,ilen].*evec ;
        jnew = vert(conn(ref1,2),:) ...
             - [jlen,jlen].*evec ;
        
        end
        
    %------------------------------------- merge i,j if near        
        near = ...
            ilen+jlen>=olen*ntol ;
        
        znew = inew(near,:) * .5 ...
             + jnew(near,:) * .5 ;
        
        inew = inew(~near,1:2) ;
        jnew = jnew(~near,1:2) ;
 
    %------------------------------------- split constraints
        zset = (1:size(znew,1))' ...
                + size(vert,1) ;
                
        iset = (1:size(inew,1))' ...
                + size(znew,1) + ...
                + size(vert,1) ;
                
        jset = (1:size(jnew,1))' ...
                + size(znew,1) + ...
                + size(inew,1) + ...
                + size(vert,1) ;
        
        set1 = num1( near);
        set2 = num1(~near);
        
        cnew = [conn( set1,1), zset
                conn( set1,2), zset
                conn( set2,1), iset
                conn( set2,2), jset
                iset, jset ] ;
        conn = [conn(~ref1,:); cnew];
            
    %------------------------------------- update vertex set    
        vert = [vert; znew(:,1:2)];
        vert = [vert; inew(:,1:2)];
        vert = [vert; jnew(:,1:2)];
           
           
        end % switch(lower(opts.kind))
    
    end

end

function [vert,conn,tria,tnum,iter] = ...
            cdtref2(vert,conn,tria,tnum, ...
                node,PSLG,part,opts,hfun,harg,iter)
%CDTREF2 constrained Delaunay-refinement for 2-simplex elem-
%nts embedded in R^2.
%   [...] = CDTREF2(...) refines the set of 2-simplex eleme-
%   nts embedded in the triangulation until all constraints 
%   are satisfied. Specifically, triangles are refined until
%   all local mesh-spacing and element-shape conditions are
%   met. Refinement proceeds according to either a Delaunay-
%   refinement or Frontal-Delaunay type approach, depending 
%   on user-settings. In either case, new steiner points are
%   introduced to split "bad" triangles - those that violate 
%   the set of prescribed constraints. In the "-DR" type pr-
%   ocess triangles are split about their circumballs. In
%   the "-FD" approach, new vertices are positioned such th-
%   at mesh-spacing and element-shape constraints are satis-
%   fied in a "locally-optimal" fashion.

    bias = +.775;

    while  (true)
    
        iter = iter + 1 ;

    %------------------------------------- build current CDT
       [vert,conn, ...
        tria,tnum]= deltri2(vert,conn, ...
                            node,PSLG, ...
                            part) ;

       [edge,tria]= tricon2(tria,conn) ;
       
        if (iter>=opts.iter),break; end

    %------------------------------------- calc. circumballs
        bal1 = cdtbal1(vert,conn) ;
        bal2 = cdtbal2(vert,edge,tria) ;
        len2 = minlen2(vert,tria) ;

        rho2 = bal2(:,+3) ./ len2 ;

    %------------------------------------- refinement scores
        scr2 = rho2 .* bal2(:,+3) ;
  
    %------------------------------------- eval. length-fun.     
        if (~isempty(hfun))
            if (isnumeric(hfun))
            fun0 = hfun * ...
              ones(size(vert,1),1);
            fun2 = hfun ;
            else
            fun0 = feval( ...
                hfun,vert,harg{:});
            fun2 = fun0(tria(:,1))...
                 + fun0(tria(:,2))...
                 + fun0(tria(:,3));
            fun2 = fun2 / +3. ;
            end
        else
            fun0 = +inf * ...
              ones(size(vert,1),1);
            fun2 = +inf ;
        end
        
        siz2 = ...
         +3. * bal2(:,3)./(fun2.*fun2) ;

    %------------------------------------- refinement queues
        ref1 = false(size(conn,1),1);
        ref2 = false(size(tria,1),1);
        
        skip = isfeat2(vert,edge,tria) ;
        
        ref2(rho2>opts.rho2* ...        %- bad rad-edge len.
                  opts.rho2) = true ;
        ref2(skip) = false ;
        ref2(siz2>opts.siz2* ...        %- bad equiv. length
                  opts.siz2) = true ;
        
        num2 = find(ref2);
      
    %------------------------------------- dump-out progess!
        if (mod(iter,opts.disp)==+0 || ...
                isempty(num2) )
            numc = size(conn,1) ;
            numt = size(tria,1) ;
            fprintf(+1, ...
            '%11i %18i %18i\n', ...
            [iter,numc,numt]) ;
        end
      
    %------------------------------------- nothing to refine
        if (isempty(num2)), break; end 
        
       [scr2,idx2] = sort( ...
            scr2(num2),'descend');
        num2 = num2(idx2);

    %------------------------------------- refine "bad" tria
        switch (lower(opts.kind))
        case 'delaunay'
    %------------------------------------- do circ-ball pt's
        new2 = zeros(length(num2),3);
        new2(:,1:2) = bal2(num2,1:2);
        new2(:,  3) = ...
            bal2(num2,3) * bias^2;
        
        
        case 'delfront'
    %-- off-centre scheme -- refine triangles by positioning
    %-- new vertices along a local segment of the voronoi
    %-- diagram, bounded by assoc. circmballs. New points
    %-- are placed to satisfy the worst of local mesh-length 
    %-- and element-shape constraints.
             
        ftri = false(length(num2),1);
        tadj = zeros(length(num2),1);

    %------------------------------------- find frontal edge
       [lmin,emin] = ...
            minlen2(vert,tria(num2,:)) ;

        epos = zeros(length(num2),1);

        for ii = +1 : length(epos)
            epos(ii) = ...
              tria(num2(ii),emin(ii)+3);
        end

    %------------------------------------- find neighb. tria
        for ii = +1 : length(epos)
        if (num2(ii)~= edge(epos(ii),3))
            tadj(ii) = edge(epos(ii),3);
        else
            tadj(ii) = edge(epos(ii),4);
        end
        end
        
    %------------------------------------- find frontal tria
        for ii = +1 : length(tadj)
        if (edge(epos(ii),5) > +0)
            ftri(ii) = true ;       
        else
            ftri(ii) = ~ref2(tadj(ii)) ;
        end
        end
        
    %------------------------------------- do circ-ball pt's
        new2 = zeros(length(num2),3);
        new2(:,1:2) = bal2(num2,1:2);
        new2(:,  3) = ...
            bal2(num2,3)*bias^2 ;
       
    %------------------------------------- locate offcentres 
        fedg = epos(ftri);
        
        emid = vert(edge(fedg,+1),:) ...
             + vert(edge(fedg,+2),:) ;
        emid = emid * +0.50 ;
        
        elen = sqrt(lmin(ftri)) ;
    
    %------------------------------------- "voro"-type dist.    
        ovec = new2(ftri,1:2)-emid ;
        olen = sqrt(sum(ovec.^2,2));
        ovec = ovec ./ [olen,olen] ;
        
        hmid = fun0(edge(fedg,+1),:) ...
             + fun0(edge(fedg,+2),:) ;
        hmid = hmid * +0.50 ;
        
    %------------------------------------- "ball"-type dist.
        rtri = elen * opts.off2 ;
        rfac = elen * 0.50;
        dsqr = rtri.^2 - rfac.^2;
        doff = rtri + ...
            sqrt(max(+0.,dsqr)) ;
        
    %------------------------------------- "size"-type dist.
        dsiz = sqrt(3.) / 2.* hmid ;
   
    %------------------------------------- bind "safe" dist.
        dist = ...
          min([dsiz,doff,olen],[],2) ;

    %------------------------------------- locate offcentres
        off2 = ...
        emid + [dist,dist] .* ovec ;
  
    %------------------------------------- iter. "size"-type
        for ioff = +1 : +4
    %------------------------------------- eval. length-fun.           
        if (~isempty(hfun))
            if (isnumeric(hfun))
            hprj = hfun * ...
              ones(size(off2,1),1);
            else
            hprj = feval( ...
                hfun,off2,harg{:});
            end
        else
            hprj = +inf * ...
              ones(size(off2,1),1);
        end

    %------------------------------------- "size"-type dist.        
        hprj = 0.5*hmid + 0.5*hprj ;
        
        dsiz = sqrt(3.) / 2.* hprj ;
        
    %------------------------------------- bind "safe" dist. 
        dist = ...
          min([dsiz,doff,olen],[],2) ;

    %------------------------------------- locate offcentres
        off2 = ...
        emid + [dist,dist] .* ovec ;
            
        end
   
    %------------------------------------- keep frontal pt's
        orad = ...    
        sqrt((elen*.5).^2 + dist.^2) ;
        
        if (any(ftri))
            new2 = ... 
        [off2(:,1:2),(bias*orad).^2] ;
        end
        
        end % switch(lower(opts.kind))

    %------------------------------------- inter.-ball dist.
       [vp,vi] = ...
          findball(new2,new2(:,1:2));

    %------------------------------------- proximity filters
        keep = true (size(new2,1),1);
        for ii = size(vp,1):-1:+1
            for ip = vp(ii,1) ...
                   : vp(ii,2)
                jj = vi(ip);
                if (jj<ii && keep(jj))           
                keep(ii) = false ;
                end 
            end
        end

        new2 = new2(keep,:);

    %------------------------------------- find encroachment
       [vp,vi] = ...
          findball(bal1,new2(:,1:2));
        
        keep = true (size(new2,1),1);
        for ii = size(vp,1):-1:+1
            for ip = vp(ii,1) ...
                   : vp(ii,2)
                jj = vi(ip);
                ref1(jj) =  true ;
                keep(ii) = false ;
            end
        end
        
        new2 = new2(keep,:);
        new1 = bal1(ref1,:);
        
    %------------------------------------- split constraints
        vnew = (1:length( ...
          find(ref1)))'+size(vert,1);
        
        cnew = [conn( ref1,1), vnew
                conn( ref1,2), vnew];
        conn = [conn(~ref1,:); cnew];
 
    %------------------------------------- update vertex set              
        vert = [vert; new1(:,1:2)];
        vert = [vert; new2(:,1:2)];
        
    end

end

function [opts] = makeopt(opts)
%MAKEOPT setup the options structure for REFINE2.

    if (~isfield(opts,'kind'))
        opts.kind = 'delfront';
    else
    if (~strcmpi(opts.kind,'delfront') && ...
        ~strcmpi(opts.kind,'delaunay') )
        error( ...
        'refine2:invalidOption','Invalid refinement KIND.'); 
    end
    end
    
    if (~isfield(opts,'iter'))
        opts.iter = +inf;
    else
    if (~isnumeric(opts.iter))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.iter)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    end
    
    if (~isfield(opts,'disp'))
        opts.disp = +10 ;
    else
    if (~isnumeric(opts.disp))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.disp)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    end
    
    if (~isfield(opts,'rho2'))
        opts.rho2 = 1.025;
    else
    if (~isnumeric(opts.rho2))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.rho2)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    end
    
    if (~isfield(opts,'off2'))
        opts.off2 = 0.933;
    else
    if (~isnumeric(opts.off2))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.off2)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    end
    
    if (~isfield(opts,'siz1'))
        opts.siz1 = 1.333;
    else
    if (~isnumeric(opts.siz1))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.siz1)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    end
    
    if (~isfield(opts,'siz2'))
        opts.siz2 = 1.300;
    else
    if (~isnumeric(opts.siz2))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.siz2)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    end
    
end



