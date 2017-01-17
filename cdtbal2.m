function [cc] = cdtbal2(pp,ee,tt)
%CDTBAL2 compute the modified circumballs associated with a 
%constrained 2-simplex Delaunay triangulation in R^2.
%   [CC] = CDTBAL2(PP,EE,TT) returns the smallest enclosing 
%   balls associated with the triangles in [PP,TT], such th-
%   at CC = [XC,YC,RC.^2]. Such balls never lie outside the 
%   boundaries of the associated CDT. See TRICON2 for info-
%   mation regarding the edge array EE.

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 14/01/2017

%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ~isnumeric(ee) || ...
        ~isnumeric(tt) )
        error('cdtbal2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ndims(ee) ~= +2 || ...
        ndims(tt) ~= +2 )
        error('cdtbal2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(pp,2)~= +2 || size(ee,2) < +5 || ...
        size(tt,2) < +6 )
        error('cdtbal2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

%----------------------------------------- calc. circumballs
    cc = tribal2(pp,tt);
    
%------------------------ replace with face-balls if smaller
    cc = minfac2(cc,pp,ee,tt,1,2,3) ;
    cc = minfac2(cc,pp,ee,tt,2,3,1) ;
    cc = minfac2(cc,pp,ee,tt,3,1,2) ; 
        
end

function [cc] = minfac2(cc,pp,ee,tt,ni,nj,nk)
%MINFAC2 modify the set of circumballs to constrain centres
%to the boundaries of the CDT.
%   [CM] = MINFAC2(CC,PP,EE,TT,NI,NJ,NK) returns the set of
%   modified circmballs CM, where any ball CC lying outside
%   the boundaries of the CDT [PP,EE,TT] is replaced by the
%   edge-centred diametric ball. [NI,NJ] are the local inde-
%   xes associated with an edge to test. NK is the local in-
%   dex of the opposite vertex.

%------------------------------------------------ edge balls
    bc = (pp(tt(:,ni),:)+pp(tt(:,nj),:))*.50 ;
    
%------------------------------------------------ edge radii
    br = sum((bc(:,1:2)-pp(tt(:,ni),:)).^2,2)...
       + sum((bc(:,1:2)-pp(tt(:,nj),:)).^2,2);
    br = br * +0.5 ;
    
%------------------------------------------------ outer edge          
    ei = tt(:,ni+3);
    ok = ee(ei,5)>0;
              
%------------------------------------------- enclosing radii
    ll = sum((bc(:,1:2)-pp(tt(:,nk),:)).^2,2);
    
%------------------------------------------- replace if min.
    ki = br >= ll & br <= cc(:,3) & ok;
    
%------------------------------------------- replace is min.
    cc(ki,1:2) = bc(ki,:);
    cc(ki,  3) = br(ki,:);

end



