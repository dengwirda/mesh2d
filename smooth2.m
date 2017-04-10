function [vert,conn,tria,tnum] = smooth2(varargin)
%SMOOTH2 "hill-climbing" mesh-smoothing for two-dimensional,
%2-simplex triangulations.
%   [VERT,EDGE,TRIA,TNUM] = SMOOTH2(VERT,EDGE,TRIA,TNUM) re-
%   turns a "smoothed" triangulation {VERT,TRIA}, incorpora-
%   ting "optimised" vertex coordinates and mesh topology.
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
%   [VERT,EDGE,TRIA,TNUM] = SMOOTH2(... ,OPTS) passes an ad-
%   ditional options structure OPTS, containing user-defined 
%   parameters, including:
%
% - OPTS.VTOL = {+1.0E-02} -- relative vertex movement tole-
%   rance, smoothing is converged when (VNEW-VERT) <= VTOL *
%   VLEN, where VLEN is a local effective length-scale.
%
% - OPTS.ITER = {+32} -- max. number of smoothing iterations
%
% - OPTS.DISP = {+ 4} -- smoothing verbosity. Set to INF for 
%   quiet execution.
%
%   See also REFINE2, DRAWSCR, TRIDEMO

%   [VERT,EDGE,TRIA,TNUM] = SMOOTH2(... ,HFUN,HARG)

%   This routine is loosely based on the DISTMESH algorithm,
%   employing a "spring-based" analogy to redistribute mesh
%   vertices. Such an approach is described in: P.O. Persson 
%   and Gilbert Strang. "A simple mesh generator in MATLAB." 
%   SIAM review 46(2) 2004, pp: 329--345. Details of the al-
%   gorithm used here are somewhat different, with an alter-
%   ative spring-based update employed, in addition to hill-
%   climbing element quality guarantees, and vertex density
%   controls.

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 05/04/2017
    
    filename = mfilename('fullpath');
    filepath = fileparts( filename );
    
    addpath([filepath,'/aabb-tree']);

%---------------------------------------------- extract args
    vert = []; conn = []; tria = []; tnum = [] ; 
    opts = []; hfun = []; harg = {};

    if (nargin>=+1), vert = varargin{1}; end
    if (nargin>=+2), conn = varargin{2}; end
    if (nargin>=+3), tria = varargin{3}; end
    if (nargin>=+4), tnum = varargin{4}; end
    if (nargin>=+5), opts = varargin{5}; end
    if (nargin>=+6), hfun = varargin{6}; end
    if (nargin>=+7), harg = varargin(7:end); end

   [opts] = makeopt(opts) ;

%---------------------------------------------- default TNUM
    if (isempty(tnum)), 
        tnum = ones(size(tria, 1), 1) ; 
    end

%---------------------------------------------- basic checks    
    if ( ~isnumeric(vert) || ...
         ~isnumeric(conn) || ...
         ~isnumeric(tria) || ...
         ~isnumeric(tnum) || ...
         ~isstruct (opts) )
        error('smooth2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(vert) ~= +2 || ...
        ndims(conn) ~= +2 || ...
        ndims(tria) ~= +2 || ...
        ndims(tnum) ~= +2 )
        error('smooth2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    if (size(vert,2)~= +2 || ...
        size(conn,2)~= +2 || ...
        size(tria,2)~= +3 || ...
        size(tnum,2)~= +1 || ...
        size(tria,1)~= size(tnum,1) )
        error('smooth2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nvrt = size(vert,1) ;

%---------------------------------------------- basic checks
    if (min(min(conn(:,1:2))) < +1 || ...
            max(max(conn(:,1:2))) > nvrt )
        error('smooth2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
    
    if (min(min(tria(:,1:3))) < +1 || ...
            max(max(tria(:,1:3))) > nvrt )
        error('smooth2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end

%---------------------------------------------- output title
    if (~isinf(opts.disp))
        fprintf(1,'\n') ;
        fprintf(1,' Smooth triangulation...\n') ;
        fprintf(1,'\n') ;
        fprintf(1,[...
' -------------------------------------------------------\n', ...
'      |ITER.|          |MOVE(X)|          |DTRI(X)|     \n', ...
' -------------------------------------------------------\n', ...
             ] ) ;
    end

%---------------------------------------------- polygon bnds    
    node = vert; PSLG = conn; part = {};

    pmax = max(tnum(:));
    for ppos = +1 : pmax
    
        tsel = tnum == ppos ;
        tcur = tria(tsel,:) ;
        
       [ecur,tcur] = tricon2(tcur) ;
    
        ebnd = ecur(:,4)==0 ;
    
        same = setset2( ...
            PSLG,ecur(ebnd,1:2));
    
        part{ppos} = find(same) ;
    
    end

%---------------------------------------------- DO MESH ITER
    tcpu.full = +0. ;
    tcpu.dtri = +0. ;
    tcpu.tcon = +0. ;
    tcpu.hfun = +0. ;
    tcpu.iter = +0. ;
    tcpu.undo = +0. ;

    tnow =  tic ;
     
    for iter = +1 : opts.iter

    %------------------------------------------ inflate adj.
        ttic = tic ;
       
       [edge,tria] = tricon2(tria,conn);

        tcpu.tcon = ...
            tcpu.tcon + toc(ttic) ;

    %------------------------------------------ vertex |deg|
        vdeg = zeros(size(vert,1),1) ;
        for epos = +1 : size(edge,1)
            
            ivrt = edge(epos,1);
            jvrt = edge(epos,2);
 
            vdeg(ivrt) = ...
                vdeg(ivrt) + 1 ;
            vdeg(jvrt) = ...
                vdeg(jvrt) + 1 ;
    
        end
        free = (vdeg == 0) ; 
        
    %------------------------------------------ compute scr.
        oscr = triscr2(vert,tria) ;
        
    %------------------------------------------ vert. iter's
        ttic =  tic ;
        
        vold = vert ;
        for isub = +1 : max(+2,min(+8,iter))
    
        %-- compute HFUN at vert/midpoints
            hvrt = evalhfn( ...
                vert,edge,hfun,harg) ;
                
            hmid = hvrt(edge(:,1),:) ...
                 + hvrt(edge(:,2),:) ;
            hmid = hmid * +.5 ;
            
        %-- calc. relative edge extensions        
            evec = vert(edge(:,2),:) ...
                 - vert(edge(:,1),:) ;
            elen = sqrt(sum(evec.^2,2));
            
            scal = +1. - elen./hmid;
            scal = min (+1.0, scal);
            scal = max (-1.0, scal);
            
        %-- projected points from each end
            ipos = vert(edge(:,1),:) ...
                 -.67*[scal,scal].*evec;
            jpos = vert(edge(:,2),:) ...
                 +.67*[scal,scal].*evec;
            
            scal = abs(scal);
            
        %-- sum contributions edge-to-vert
            vsum = zeros(size(vert,1),1) ;
            vsum = vsum + eps^.8;
            vnew = vert * eps^.8;
            for epos = +1 : size(edge,1)
            
                ivrt = edge(epos,1);
                jvrt = edge(epos,2);

                wval = scal(epos,1);
               %wval = +1. ;

                vsum(ivrt,1) = ...
                vsum(ivrt,1) + wval;
                vsum(jvrt,1) = ...
                vsum(jvrt,1) + wval;

                vnew(ivrt,1) = ...
                vnew(ivrt,1) + ...
                    wval * ipos(epos,1);
                vnew(ivrt,2) = ...
                vnew(ivrt,2) + ...
                    wval * ipos(epos,2);
            
                vnew(jvrt,1) = ...
                vnew(jvrt,1) + ...
                    wval * jpos(epos,1);
                vnew(jvrt,2) = ...
                vnew(jvrt,2) + ...
                    wval * jpos(epos,2);
            
            end
            vnew = vnew ./ [vsum, vsum];
            
        %-- fixed points. edge projection?
            vnew(conn(:),1:2) = ...
            vert(conn(:),1:2) ;
            
        %-- reset for the next local iter.         
            vert = vnew ;
    
        end
 
        tcpu.iter = ...
            tcpu.iter + toc(ttic) ;
    
    %------------------------------------------ hill-climber
        ttic = tic ;
        
    %-- find worst tria adj. to each vert.
        lscr = ones (size(vert,1),1) ;
        ilow = zeros(size(vert,1),1) ;
        for tpos = +1 : size(tria,1)

        ivrt = tria(tpos,1) ;
        jvrt = tria(tpos,2) ;
        kvrt = tria(tpos,3) ;
        
        if (lscr(ivrt) > oscr(tpos))
            ilow(ivrt) = tpos;
            lscr(ivrt) = oscr(tpos);
        end
        if (lscr(jvrt) > oscr(tpos))
            ilow(jvrt) = tpos;
            lscr(jvrt) = oscr(tpos);
        end
        if (lscr(kvrt) > oscr(tpos))
            ilow(kvrt) = tpos;
            lscr(kvrt) = oscr(tpos);
        end
        
        end
        tlow = false(size(tria,1),1) ;
        tlow(ilow(ilow > +0)) = true ;
        
    %-- unwind vert. upadte if score lower
        nscr = ones(size(tria,1),1);
        btri = true(size(tria,1),1);
        
        for undo = +1 : size(vert,+1)
        
            nscr(btri) = triscr2( ...
                vert,tria(btri,:)) ;
      
        %-- TRUE if tria needs "unwinding" 
            smin = +.70 ;
            smax = +.90 ;
            sdel = .025 ;
        
            stol = smin+iter*sdel;
            stol = min (smax,stol) ;
            
            btri = nscr <= stol ...
                 & nscr <  oscr ...
                 & tlow ;
             
            if (~any(btri)), break; end
             
        %-- relax toward old vert. coord's
            ivrt = ...
              unique(tria(btri,1:3));
            
            bvrt = false(size(vert,1),1);
            bvrt(ivrt) = true;
            
            bnew =  +.80 ^ undo ;
            bold =  +1.0 - bnew ;
            
            vert(bvrt,:) = ...
                bold * vold(bvrt,:) ... 
              + bnew * vert(bvrt,:) ;
      
            btri = ...
                any(bvrt(tria(:,1:3)),2);
        
        end
    
        tcpu.undo = ...
            tcpu.undo + toc(ttic) ;
    
    %------------------------------------- test convergence!
        vdel = sqrt(sum((vert-vold).^2, 2)) ;
    
    
        % todo!!
        % delete small edges
        % refine large edges
   
        
    %------------------------------------- |deg|-based prune
        keep = false(size(vert,1),1);
        keep(vdeg>+4) = true;
        keep(conn(:)) = true;
        keep(free(:)) = true;
        
    %------------------------------------- reindex vert/conn 
        redo = zeros(size(vert,1),1);
        redo(keep) = +1;
        redo = cumsum(redo);
        conn = redo  (conn);
        
        vert = vert(keep,:);

  
    %------------------------------------- build current CDT
        ttic = tic ;
       
       [vert,conn,tria,tnum] = ...
          deltri2(vert,conn,node,PSLG,part) ;
          
        tcpu.dtri = ...
            tcpu.dtri + toc(ttic) ;
        
    %------------------------------------- dump-out progess!
        if (mod(iter,opts.disp)==+0)
            move = abs (vdel) > ...
                opts.vtol * hvrt;
            numm = sum (move,1) ;
            numt = size(tria,1) ;
            fprintf(+1, ...
            '%11i %18i %18i\n', ...
            [iter,numm,numt]) ;
        end
        
    %------------------------------------- loop convergence!
        if (all(abs(vdel)<=opts.vtol*hvrt))
            break ;
        end
        
    end

    tria = tria(:, 1:3) ;
    
    tcpu.full = ...
        tcpu.full + toc(tnow) ;

    if (opts.dbug)
    %------------------------------------- print debug timer 
        fprintf(1,'\n') ;
        fprintf(1,' Mesh smoothing timer...\n');
        fprintf(1,'\n') ;
        fprintf(1, ...
        ' FULL: %f \n', tcpu.full);
        fprintf(1, ...
        ' DTRI: %f \n', tcpu.dtri);
        fprintf(1, ...
        ' TCON: %f \n', tcpu.tcon);
        fprintf(1, ...
        ' HFUN: %f \n', tcpu.hfun);
        fprintf(1, ...
        ' ITER: %f \n', tcpu.iter);
        fprintf(1, ...
        ' UNDO: %f \n', tcpu.undo);   
        fprintf(1,'\n') ;
    end
    
    if (~isinf(opts.disp)), fprintf(1,'\n'); end

end

function [hvrt] = evalhfn(vert,edge,hfun,harg)
%EVALHFN eval. the spacing-fun. at mesh vertices.

    if (~isempty (hfun))
        if (isnumeric(hfun))
            hvrt = hfun * ...
              ones(size(vert,1),1) ;
        else
            hvrt = feval( ...
                hfun,vert,harg{:}) ;
        end
    else
    
%-- if no HFUN, HFUN is "weighted" edge-len. at vertices
    
        hvrt = zeros(size(vert,1),1) ;
        wval = zeros(size(vert,1),1) ;
    
        evec = vert(edge(:,2),:) - ...
               vert(edge(:,1),:) ;
        elen = sqrt(sum(evec.^2,2));
        
        for epos = +1 : size(edge,1)
            
            ivrt = edge(epos,1);
            jvrt = edge(epos,2);
 
           %wcur = +1./elen(epos) ;
            wcur = +1. ;
            
            hvrt(ivrt) = ...
            hvrt(ivrt) + ...
                wcur * elen(epos) ;
            hvrt(jvrt) = ...
            hvrt(jvrt) + ...
                wcur * elen(epos) ;
            
            wval(ivrt) = ...
            wval(ivrt) + wcur;
            wval(jvrt) = ...
            wval(jvrt) + wcur;
                
        end
        hvrt = hvrt ./ max(eps,wval) ;
        
        free = true (size(vert,1),1) ;
        free(edge(:,1)) = false;
        free(edge(:,2)) = false;
        
        hvrt(free)      =  +inf;
        
    end

end

function [opts] = makeopt(opts)
%MAKEOPT setup the options structure for SMOOTH2.
    
    if (~isfield(opts,'iter'))
        opts.iter = +32 ;
    else
    if (~isnumeric(opts.iter))
        error('smooth2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.iter)~= +1)
        error('smooth2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.iter <= +0)
        error('smooth2:invalidOptionValues', ...
            'Invalid OPT.ITER selection.') ;
    end
    end
    
    if (~isfield(opts,'disp'))
        opts.disp = + 4;
    else
    if (~isnumeric(opts.disp))
        error('smooth2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.disp)~= +1)
        error('smooth2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.disp <= +0)
        error('smooth2:invalidOptionValues', ...
            'Invalid OPT.DISP selection.') ;
    end
    end
    
    if (~isfield(opts,'vtol'))
        opts.vtol = +1.0E-02;
    else
    if (~isnumeric(opts.vtol))
        error('smooth2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.vtol)~= +1)
        error('smooth2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.vtol <= 0.)
        error('smooth2:invalidOptionValues', ...
            'Invalid OPT.VTOL selection.') ;
    end
    end
    
    if (~isfield(opts,'dbug'))
        opts.dbug = false;
    else
    if (~islogical(opts.dbug))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.dbug)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    end
    
end



