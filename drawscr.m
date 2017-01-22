function drawscr(vert,conn,tria,tnum)
%DRAWSCR draw quality-metrics for a 2-simplex triangulation
%embedded in the two-dimensional plane.
%   DRAWSCR(VERT,EDGE,TRIA,TNUM) draws histograms of quality
%   metrics for the triangulation.
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
%   See also DRAWTRI, REFINE2, SMOOTH2

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 17/01/2017

%---------------------------------------------- basic checks    
    if ( ~isnumeric(vert) || ...
         ~isnumeric(conn) || ...
         ~isnumeric(tria) || ...
         ~isnumeric(tnum) )
        error('drawscr:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(vert) ~= +2 || ...
        ndims(conn) ~= +2 || ...
        ndims(tria) ~= +2 || ...
        ndims(tnum) ~= +2 )
        error('drawscr:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(vert,2)~= +2 || ...
        size(conn,2) < +2 || ...
        size(tria,2) < +3 )
        error('drawscr:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nvrt = size(vert,1) ;

%---------------------------------------------- basic checks
    if (min(min(conn(:,1:2))) < +1 || ...
            max(max(conn(:,1:2))) > nvrt )
        error('drawscr:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
 
    if (min(min(tria(:,1:3))) < +1 || ...
            max(max(tria(:,1:3))) > nvrt )
        error('drawscr:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end

%-- borrowed from the JIGSAW library!

    dolabel = true ;

%-- draw sub-axes directly -- sub-plot gives
%-- silly inconsistent spacing...!
    
    axpos21 = [.125,.60,.80,.30] ;
    axpos22 = [.125,.15,.80,.30] ;
    
    axpos31 = [.125,.75,.80,.15] ;
    axpos32 = [.125,.45,.80,.15] ;
    axpos33 = [.125,.15,.80,.15] ;
    
%-- draw cost histograms for 2-tria elements
    figure;
    set(gcf,'color','w','position',[128,128,640,300]);
   %if (isfield(cost.tria3,'hfunc') )
    if (false)
    
%-- have size-func data
    axes('position',axpos31); hold on;
    scrhist(triscr2(vert,tria),'tria3');
    if (dolabel)
    title('Quality metrics (TRIA-2)');
    end
    axes('position',axpos32); hold on;
    anghist(triang2(vert,tria),'tria3');
    axes('position',axpos33); hold on;
   %hfnhist(trihfn2(vert,tria),'tria3');
    
    else
    
%-- null size-func data
    axes('position',axpos21); hold on;
    scrhist(triscr2(vert,tria),'tria3');
    if (dolabel)
    title('Quality metrics (TRIA-2)');
    end
    axes('position',axpos22); hold on;
    anghist(triang2(vert,tria),'tria3');
    
    end
      
end

function anghist(ad,ty)
%ANGHIST draw histogram for "angle" quality-metric.

    ad = ad(:);
    be = linspace(0.,180.,91);
    bm =(be(1:end-1)+be(2:end))/2.;
    hc = histc(ad,be);
    
    switch (ty)
    case 'tria4'
        poor = bm <  10.  | bm >= 160. ;
        okay =(bm >= 10.  & bm <  20. )| ...
              (bm >= 140. & bm <  160.);
        good =(bm >= 20.  & bm <  30. )| ...
              (bm >= 120. & bm <  140.);
        best = bm >= 30.  & bm <  120. ;
        
    case 'tria3'
        poor = bm <  15.  | bm >= 150. ;
        okay =(bm >= 15.  & bm <  30. )| ...
              (bm >= 120. & bm <  150.);
        good =(bm >= 30.  & bm <  45. )| ...
              (bm >= 90.  & bm <  120.);
        best = bm >= 45.  & bm <  90.  ;
    
    end
    
    r = [.85,.00,.00] ; y = [1.0,.95,.00] ;
    g = [.00,.90,.00] ; k = [.60,.60,.60] ;
    
    bar(bm(poor),hc(poor),1.05,...
        'facecolor',r,'edgecolor',r) ;
    bar(bm(okay),hc(okay),1.05,...
        'facecolor',y,'edgecolor',y) ;
    bar(bm(good),hc(good),1.05,...
        'facecolor',g,'edgecolor',g) ;
    bar(bm(best),hc(best),1.05,...
        'facecolor',k,'edgecolor',k) ;
    
    axis tight;
    set(gca,'ycolor', get(gca,'color'),'ytick',[],...
        'xtick',0:30:180,'layer','top','fontsize',...
            18,'linewidth',2.,'ticklength',[.025,.025],...
                'box','off','xlim',[0.,180.]) ;
    
    mina = max(1.000,min(ad)); %%!! so that axes don't obscure!
    maxa = min(179.0,max(ad));
    
    line([ mina, mina],...
        [0,max(hc)],'color','r','linewidth',1.5);
    line([ maxa, maxa],...
        [0,max(hc)],'color','r','linewidth',1.5);
    
    if ( mina > 25.0)
        text(mina-1.8,.9*max(hc),num2str(min(ad),'%16.1f'),...
            'horizontalalignment',...
                'right','fontsize',22) ;
    else
        text(mina+1.8,.9*max(hc),num2str(min(ad),'%16.1f'),...
            'horizontalalignment',...
                'left' ,'fontsize',22) ;
    end
    
    if ( maxa < 140.)
        text(maxa+1.8,.9*max(hc),num2str(max(ad),'%16.1f'),...
            'horizontalalignment',...
                'left' ,'fontsize',22) ;    
    else
        text(maxa-1.8,.9*max(hc),num2str(max(ad),'%16.1f'),...
            'horizontalalignment',...
                'right','fontsize',22) ;     
    end
   
    switch (ty)
    case 'tria4'
        text(-9.0,0.0,'$\theta_{\tau}$',...
            'horizontalalignment','right',...
                'fontsize',28,'interpreter','latex') ;
        
    case 'tria3'
        text(-9.0,0.0,'$\theta_{f}$',...
            'horizontalalignment','right',...
                'fontsize',28,'interpreter','latex') ;
        
    end
    
end

function scrhist(sc,ty)
%SCRHIST draw histogram for "score" quality-metric.

    sc = sc(:);
    be = linspace(0.,1.,101);
    bm = (be(1:end-1)+be(2:end)) / 2.;
    hc = histc(sc,be);

    switch (ty)   
    case 'tria4'
        poor = bm <  .25 ;
        okay = bm >= .25 & bm <  .50 ;
        good = bm >= .50 & bm <  .75 ;
        best = bm >= .75 ;
    
    case 'tria3'
        poor = bm <  .30 ;
        okay = bm >= .30 & bm <  .60 ;
        good = bm >= .60 & bm <  .90 ;
        best = bm >= .90 ;
        
    end
 
    r = [.85,.00,.00] ; y = [1.0,.95,.00] ;
    g = [.00,.90,.00] ; k = [.60,.60,.60] ;
    
    bar(bm(poor),hc(poor),1.05,...
        'facecolor',r,'edgecolor',r) ;
    bar(bm(okay),hc(okay),1.05,...
        'facecolor',y,'edgecolor',y) ;
    bar(bm(good),hc(good),1.05,...
        'facecolor',g,'edgecolor',g) ;
    bar(bm(best),hc(best),1.05,...
        'facecolor',k,'edgecolor',k) ;
    
    axis tight;    
    set(gca,'ycolor', get(gca,'color'),'ytick',[],...
        'xtick',.0:.2:1.,'layer','top','fontsize',...
            18,'linewidth',2.,'ticklength',[.025,.025],...
                'box','off','xlim',[0.,1.]) ;
    
    mins = max(0.010,min(sc)); %%!! so that axes don't obscure!
    maxs = min(0.990,max(sc));
    
    line([ mins, mins],...
        [0,max(hc)],'color','r','linewidth',1.5);
    line([mean(sc),mean(sc)],...
        [0,max(hc)],'color','r','linewidth',1.5);
    
    if ( mins > .4)
        text(mins-.01,.9*max(hc),num2str(min(sc),'%16.3f'),...
            'horizontalalignment',...
                'right','fontsize',22) ;
    else
        text(mins+.01,.9*max(hc),num2str(min(sc),'%16.3f'),...
            'horizontalalignment',...
                'left' ,'fontsize',22) ;
    end
    
    text(mean(sc)-.01,.9*max(hc),num2str(mean(sc),'%16.3f'),...
        'horizontalalignment','right','fontsize',22) ;
    
    switch (ty)
    case 'tria4'
        text(-.05,0.0,'$v_{\tau}$',...
            'horizontalalignment','right',...
                'fontsize',28,'interpreter','latex') ;
        
    case 'tria3'
        text(-.05,0.0,'$a_{f}$',...
            'horizontalalignment','right',...
                'fontsize',28,'interpreter','latex') ;
        
    end
    
end

function hfnhist(hf,ty)
%HFNHIST draw histogram for "hfunc" quality-metric.

    be = linspace(0.,2.,101);
    bm = (be(1:end-1)+be(2:end)) / 2.;
    hc = histc(hf,be);

    poor = bm <  .40 | bm >= 1.6  ;
    okay =(bm >= .40 & bm <  .60 )| ...
          (bm >= 1.4 & bm <  1.6 );
    good =(bm >= .60 & bm <  .80 )| ...
          (bm >= 1.2 & bm <  1.4 );
    best = bm >= .80 & bm <  1.2 ;
 
    r = [.85,.00,.00] ; y = [1.0,.95,.00] ;
    g = [.00,.90,.00] ; k = [.60,.60,.60] ;
    
    bar(bm(poor),hc(poor),1.05,...
        'facecolor',r,'edgecolor',r) ;
    bar(bm(okay),hc(okay),1.05,...
        'facecolor',y,'edgecolor',y) ;
    bar(bm(good),hc(good),1.05,...
        'facecolor',g,'edgecolor',g) ;
    bar(bm(best),hc(best),1.05,...
        'facecolor',k,'edgecolor',k) ;
    
    axis tight; 
    set(gca,'ycolor', get(gca,'color'),'ytick',[],...
        'xtick',.0:.5:2.,'layer','top','fontsize',...
            18,'linewidth',2.,'ticklength',[.025,.025],...
                'box','off','xlim',[0.,2.]);
    
    line([mean(hf),mean(hf)],...
        [0,max(hc)],'color','r','linewidth',1.5);
    
    text(mean(hf)+.02,.9*max(hc),num2str(mean(hf),'%16.2f'),...
        'horizontalalignment','left','fontsize',22);
    
    text(-0.100,0.0,'$h_{r}$','horizontalalignment','right',...
        'fontsize',28,'interpreter','latex');
    
end



