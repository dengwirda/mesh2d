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
%   Last updated    : 26/01/2017
    
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
    for iter = +1 : opts.iter

    %------------------------------------------ inflate adj.
       [edge,tria] = tricon2(tria,conn);

    %------------------------------------------ vertex |deg|
        vadj = zeros(size(vert,1),1) ;
        for epos = +1 : size(edge,1)
            
            ivrt = edge(epos,1);
            jvrt = edge(epos,2);
 
            vadj(ivrt) = ...
                vadj(ivrt) + 1 ;
            vadj(jvrt) = ...
                vadj(jvrt) + 1 ;
    
        end
        free = (vadj == 0) ; 
       
    %------------------------------------------ compute HFUN 
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
     
                wcur = +1./elen(epos) ;
               %wcur = +1. ;
                
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
            hvrt(free) =  +inf ;
            
        end
        
        hmid =(hvrt(edge(:,1),:) ...
             + hvrt(edge(:,2),:)) * 0.50 ;
        
    %------------------------------------------ compute scr.
        oscr = triscr2(vert,tria) ;
        
    %------------------------------------------ vert. iter's
        vold = vert ;
        for isub = +1 : max(+2,min(+8,iter))
            
            evec = vert(edge(:,2),:) - ...
                   vert(edge(:,1),:) ;
            elen = sqrt(sum(evec.^2,2));
    
        %-- calc. relative edge extensions        
            scal = +1. - elen./hmid;
            scal = min (+0.5, scal);
            scal = max (-0.5, scal);
            
        %-- projected points from each end
            ipos = vert(edge(:,1),:) ...
                 - 1.*[scal,scal].*evec;
            jpos = vert(edge(:,2),:) ...
                 + 1.*[scal,scal].*evec;
            
        %-- sum contributions edge-to-vert
            vnew = vert ;
            for epos = +1 : size(edge,1)
            
                ivrt = edge(epos,1);
                jvrt = edge(epos,2);
            
                vnew(ivrt,1) = ...
                vnew(ivrt,1)+ipos(epos,1);
                vnew(ivrt,2) = ...
                vnew(ivrt,2)+ipos(epos,2);
                
                vnew(jvrt,1) = ...
                vnew(jvrt,1)+jpos(epos,1);
                vnew(jvrt,2) = ...
                vnew(jvrt,2)+jpos(epos,2);
            
            end
            vnew = vnew./[vadj+1,vadj+1] ;
            
        %-- fixed points. edge projection?
            vnew(conn(:),1:2) = ...
            vert(conn(:),1:2) ;
            
        %-- reset for the next local iter.         
            vert = vnew ;
    
        end
    
    %------------------------------------------ hill-climber
        for undo = +1 : size(vert,+1)
        
            nscr = triscr2(vert,tria) ;
        
        %-- find worst tria for each vert.
            lscr = ones(size(vert,1),1) ;
            ilow = ones(size(vert,1),1) ;
            for tpos = 1 : size(tria,1)

            ivrt = tria(tpos,1) ;
            jvrt = tria(tpos,2) ;
            kvrt = tria(tpos,3) ;
            
            if (lscr(ivrt) > nscr(tpos))
                ilow(ivrt) = tpos;
                lscr(ivrt) = nscr(tpos);
            end
            if (lscr(jvrt) > nscr(tpos))
                ilow(jvrt) = tpos;
                lscr(jvrt) = nscr(tpos);
            end
            if (lscr(kvrt) > nscr(tpos))
                ilow(kvrt) = tpos;
                lscr(kvrt) = nscr(tpos);
            end
            
            end
            tlow = false(size(tria,1),1);
            tlow(ilow(ilow > 0)) = true ;
      
        %-- TRUE if tria needs "unwinding" 
            stol = +.75 ;   
            btri = nscr <= stol ...
                 & nscr <  oscr ...
                 & tlow ;
             
            if (~any(btri)), break; end
             
        %-- relax toward old vert. coord's
            bvrt = ...
              unique(tria(btri,1:3));
            
            bnew =  +.67 ^ undo ;
            bold =  +1.0 - bnew ;
            
            vert(bvrt,:) = ...
                bold * vold(bvrt,:) ... 
              + bnew * vert(bvrt,:) ;
        
        end
    
    %------------------------------------- test convergence!
        vdel = sqrt(sum((vert-vold).^2, 2)) ;
    
    
        % todo!!
        % delete small edges
        % refine large edges
  
  
    %------------------------------------- build current CDT
       [vert,conn,tria,tnum] = ...
          deltri2(vert,conn,node,PSLG,part) ;
        
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

    tria = tria(:,1:3) ;
    
    if (~isinf(opts.disp)), fprintf(1,'\n'); end

end

function [opts] = makeopt(opts)
%MAKEOPT setup the options structure for SMOOTH2.
    
    if (~isfield(opts,'iter'))
        opts.iter = +32;
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
    
end



