function tridemo(demo)
%TRIDEMO run various triangulation demos for MESH2D.
%   TRIDEMO(N) runs the N-TH demo problem. The following de-
%   mo problems are currently available:
%
% - DEMO-1: explore the impact of the "radius-edge" thresho-
%   ld (RHO2) on mesh density/quality.
%
% - DEMO-2: explore the impact of the "Frontal-Delaunay" vs.
%   "Delaunay-refinement " algorithms. 
%
% - DEMO-3 explore impact of user-defined mesh-size constra-
%   ints.
%
% - DEMO-4 explore impact of "hill-climbing" mesh optimisat-
%   ions.
%
% - DEMO-5 assemble triangulations for multi-part geometry 
%   definitions.
%
% - DEMO-6 explore impact of user-defined mesh-size constra-
%   ints.
%
% - DEMO-7 larger-scale problem, mesh refinement + optimisa-
%   tion. 
%
% - DEMO-8 medium-scale problem, mesh refinement + optimisa-
%   tion. 
%
%   See also REFINE2, SMOOTH2

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 23/01/2017

    close all;

    switch (demo)
        case 1, demo1();
        case 2, demo2();
        case 3, demo3();
        case 4, demo4();
        case 5, demo5();
        case 6, demo6();
        case 7, demo7();
        case 8, demo8();
            
        otherwise
        error('tridemo:invalidSelection', 'Invalid selection!');
    end

end

function demo1
%DEMO1 explore impact of RHO2 threshold on mesh density/qua-
%lity.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/lake.msh'];

   [node,edge] = triread( meshfile );
 
    fprintf(1, [ ...
' The REFINE2 routine can be used to build guaranteed-quality \n', ...
' Delaunay triangulations for general polygonal geometries in \n', ...
' the two-dimensional plane. The "quality" of elements in the \n', ...
' triangulation can be controlled using the "radius-edge" bo- \n', ...
' und RHO2. \n', ...
        ] ) ;
        
%---------------------------------------------- RHO2 = +1.50   
    fprintf(1, ' \n') ;
    fprintf(1, [ ...
' Setting large values for RHO2, (RHO2 = 1.50 here) generates \n', ...
' sparse triangulations with poor worst-case angle bounds. \n', ...
        ] ) ;
    
    opts.kind = 'delaunay';
    opts.rho2 = +1.50 ;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[]  ,opts) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    title(['TRIA-MESH: RHO2<=+1.50, |TRIA|=' , ...
        num2str(size(tria,1))]) ;
    
%---------------------------------------------- RHO2 = +1.00    
    fprintf(1, [ ...
' Setting small values for RHO2, (RHO2 = 1.00 here) generates \n', ...
' dense triangulations with good worst-case angle bounds. \n', ...
        ] ) ;
    
    opts.kind = 'delaunay';
    opts.rho2 = +1.00 ;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[]  ,opts) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    title(['TRIA-MESH: RHO2<=+1.00, |TRIA|=' , ...
        num2str(size(tria,1))]) ;
        
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
    
end

function demo2
%DEMO2 explore impact of refinement "KIND" on mesh quality/-
%density.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/lake.msh'];

   [node,edge] = triread( meshfile );
 
    fprintf(1, [ ...
' The REFINE2 routine supports two Delaunay-based refinement  \n', ...
' algorithms: a "standard" Delaunay-refinement type approach, \n', ...
' and a "Frontal-Delaunay" technique. For problems constrain- \n', ...
' ed by element "quality" alone, the Frontal-Delaunay approa- \n', ...
' ch typically produces sigificantly sparser meshes. in both  \n', ...
' cases, the same worst-case element quality bounds are sati- \n', ...
' fied in a guaranteed manner.'
        ] ) ;
 
%---------------------------------------------- = "DELAUNAY"
    opts.kind = 'delaunay';
    opts.rho2 = +1.00 ;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[]  ,opts) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert , ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    title(['TRIA-MESH: KIND=DELAUNAY, |TRIA|=', ...
        num2str(size(tria,1))]) ;
    
%---------------------------------------------- = "DELFRONT"
    opts.kind = 'delfront';
    opts.rho2 = +1.00 ;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[]  ,opts) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert , ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    title(['TRIA-MESH: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
    
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
    
end

function demo3
%DEMO3 explore impact of user-defined mesh-size constraints.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/airfoil.msh'];

   [node,edge] = triread( meshfile );
 
    fprintf(1, [ ...
' Additionally, the REFINE2 routine supports size-driven ref- \n', ...
' inement, producing meshes that satisfy constraints on elem- \n', ...
' ent edge-lengths. The LFSHFN2 routine can be used to create \n', ...
' mesh-size functions based on an estimate of the "local-fea- \n', ...
' ture-size" associated with a polygonal domain. The Frontal- \n', ...
' Delaunay refinement algorithm discussed in DEMO-2 is espec- \n', ...
' ially good at generating high-quality triangulations in the \n', ...
' presence of mesh-size constraints. \n', ...
        ] ) ;
 
%---------------------------------------------- do size-fun. 
    olfs.dhdx = +0.15;
 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge, ...
                    []  ,olfs) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
    figure;
    patch('faces',tlfs(:,1:3),'vertices',vlfs , ...
        'facevertexcdata' , hlfs, ...
        'facecolor','interp', ...
        'edgecolor','none') ;
    hold on; axis image off;
    title(['MESH-SIZE: KIND=DELAUNAY, |TRIA|=', ...
        num2str(size(tlfs,1))]) ;
 
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun , ...
                    vlfs,tlfs,slfs,hlfs);
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert , ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    title(['TRIA-MESH: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
    
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
    
end

function demo4
%DEMO4 explore impact of "hill-climbing" mesh optimisations.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/airfoil.msh'];

   [node,edge] = triread( meshfile );
 
    fprintf(1, [ ...
' The SMOOTH2 routine provides iterative mesh "smoothing" ca- \n', ...
' pabilities, seeking to improve triangulation quality by ad- \n', ...
' justing the vertex positions and mesh topology. Specifical- \n', ...
' ly, a "hill-climbing" type optimisation is implemented, gu- \n', ...
' aranteeing that mesh-quality is improved monotonically. The \n', ...
' DRAWSCR routine provides detailed analysis of triangulation \n', ...
' quality, plotting histograms of various quality metrics. \n', ...
        ] ) ;
 
%---------------------------------------------- do size-fun. 
    olfs.dhdx = +0.15;
 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge, ...
                    []  ,olfs) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun , ...
                    vlfs,tlfs,slfs,hlfs);
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert , ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    title(['MESH-REF.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
        
%---------------------------------------------- do mesh-opt.
   [vnew,enew, ...
    tnew,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tnew(:,1:3),'vertices',vnew , ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tnew,1))]) ;
           
    drawscr(vert,etri,tria,tnum);
    drawscr(vnew,enew,tnew,tnum);
           
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
        
    set(figure(3),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;
    set(figure(4),'units','normalized', ...
        'position',[.35,.05,.30,.35]) ;
    
end

function demo5
%DEMO5 assemble triangulations for multi-part geometry defi-
%nitions.

    fprintf(1, [ ...
' Both the REFINE2 and SMOOTH2 routines also support "multi-  \n', ...
' part" geometry definitions -- generating conforming triang- \n', ...
' ulations that conform to internal and external constraints. \n', ...
        ] ) ;

%---------------------------------------------- create geom.
    nod1 = [
        -1., -1.; +1., -1.
        +1., +1.; -1., +1.
        ] ;
    edg1 = [
         1 ,  2 ;  2 ,  3
         3 ,  4 ;  4 ,  1
        ] ;
    edg1(:,3) = +0;
    
    
    nod2 = [
        +.1, +0.; +.8, +0.
        +.8, +.8; +.1, +.8
        ] ;
    edg2 = [
         1 ,  2 ;  2 ,  3
         3 ,  4 ;  4 ,  1
        ] ;
    edg2(:,3) = +1;
        
        
    adel = 2.*pi / +64 ;
    amin = 0.*pi ;
    amax = 2.*pi - adel;
    
    xcir = +.33 * ...
        cos(amin:adel:amax)';
    ycir = +.33 * ...
        sin(amin:adel:amax)';
    xcir = xcir - .33;
    ycir = ycir - .25;
    ncir = [xcir,ycir] ;
    numc = size(ncir,1);
        
    ecir(:,1) = ...
        [(1:numc-1)'; numc] ;
    ecir(:,2) = ...
        [(2:numc-0)'; +1  ] ;
    ecir(:,3) = +2;
    
    edg2(:,1:2) = ...
    edg2(:,1:2)+size(nod1,1);
    edge = [edg1; edg2];
    node = [nod1; nod2];
        
    ecir(:,1:2) = ...
    ecir(:,1:2)+size(node,1);
    edge = [edge; ecir];
    node = [node; ncir];
    
    part{1} = [ ...
        find(edge(:,3) == 0) 
        find(edge(:,3) == 1)
        find(edge(:,3) == 2)
        ] ;
    part{2} = [ ...
        find(edge(:,3) == 1)
        ] ;
    part{3} = [ ...
        find(edge(:,3) == 2)
        ] ;
        
    edge = edge(:,1:2) ;
    
%---------------------------------------------- do size-fun.
    hmax = +0.075;
 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge, ...
                    part) ;
    
    hlfs = min(hmax,hlfs) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;

   [vert,etri, ...
    tria,tnum] = refine2(node,edge,part,  [], ...
                         hfun, ...
                         vlfs,tlfs,slfs,hlfs) ;
                         
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tria(tnum==1,1:3), ...
        'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    patch('faces',tria(tnum==2,1:3), ...
        'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor','b') ;
    patch('faces',tria(tnum==3,1:3), ...
        'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor','m') ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
    
    figure;
    patch('faces',tlfs(:,1:3),'vertices',vlfs , ...
        'facevertexcdata' , hlfs, ...
        'facecolor','interp', ...
        'edgecolor','none') ;
    hold on; axis image off;
    title(['MESH-SIZE: KIND=DELAUNAY, |TRIA|=', ...
        num2str(size(tlfs,1))]) ;
        
    drawscr(vert,etri,tria,tnum);
           
    drawnow;
        
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ; 
    set(figure(3),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;

end

function demo6
%DEMO6 explore impact of "hill-climbing" mesh optimisations.

%---------------------------------------------- create geom.
    node = [
        -1., -1.; +3., -1.
        +3., +1.; -1., +1.
        ] ;
    edge = [
         1 ,  2 ;  2 ,  3
         3 ,  4 ;  4 ,  1
        ] ;
        
    adel = 2.*pi / +64 ;
    amin = 0.*pi ;
    amax = 2.*pi - adel;
    
    xcir = +.20 * ...
        cos(amin:adel:amax)';
    ycir = +.20 * ...
        sin(amin:adel:amax)';
    ncir = [xcir,ycir] ;    
    numc = size(ncir,1);
        
    ecir(:,1) = ...
        [(1:numc-1)'; numc] ;
    ecir(:,2) = ...
        [(2:numc-0)'; +1  ] ;
        
    ecir = ecir+size(node,1);
    edge = [edge; ecir];
    node = [node; ncir];
    
%---------------------------------------------- do mesh-gen.
    hfun = @hfun6 ;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun);
    
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert , ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
        
    figure;
    patch('faces',tria(:,1:3),'vertices',vert , ...
        'facevertexcdata' , hfun6(vert), ...
        'facecolor','interp', ...
        'edgecolor','none') ;
    hold on; axis image off;
    title('MESH-SIZE function.');
   
    drawscr(vert,etri,tria,tnum);
           
    drawnow;
        
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ; 
    set(figure(3),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;  
  
end

function [hfun] = hfun6(test)
%HFUN5 user-defined mesh-size function for DEMO5.

    hmax = +.05 ;
    hmin = +.01 ;

    xmid = +0.0 ;
    ymid = +0.0 ;
    
    hcir = exp( -1.*(test(:,1)-xmid).^2 ...
                -2.*(test(:,2)-ymid).^2 ) ;

    hfun = hmax - (hmax-hmin) * hcir  ;

end

function demo7
%DEMO7 larger-scale problem, mesh refinement + optimisation.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/islands.msh'];

   [node,edge] = triread( meshfile );
 
%---------------------------------------------- do size-fun. 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun , ...
                    vlfs,tlfs,slfs,hlfs);
        
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert , ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
           
    drawscr(vert,etri,tria,tnum);
           
    drawnow;
    
end

function demo8
%DEMO8 medium-scale problem, mesh refinement + optimisation.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/river.msh'];

   [node,edge] = triread( meshfile );
 
%---------------------------------------------- do size-fun. 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun , ...
                    vlfs,tlfs,slfs,hlfs);
        
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert , ...
        'facecolor','w', ...
        'edgecolor','k') ;
    hold on; axis image off;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
           
    drawscr(vert,etri,tria,tnum);
           
    drawnow;
    
end



