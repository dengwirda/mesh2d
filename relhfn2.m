function [hrel] = relhfn2(pp,tt,hv)
%RELHFN2 calc. relative edge-length for a 2-simplex triangu-
%lation embedded in the two-dimensional plane.
%   [HREL] = RELHFN2(VERT,TRIA) returns the rel. edge-len.,
%   indicating conformance to the imposed mesh-spacing cons-
%   traints, where HREL is a E-by-1 vector, VERT is a V-by-2 
%   array of XY coordinates, and TRIA is a T-by-3 array of
%   vertex indexing, where each row defines a triangle, such 
%   that VERT(TRIA(II,1),:), VERT(TRIA(II,2),:) and VERT(
%   TRIA(II,3),:) are the coordinates of the II-TH triangle.
%
%   See also TRISCR2, TRIAREA, TRIANG2, TRIBAL2

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 11/07/2017

%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ~isnumeric(tt) ||   ...
        ~isnumeric(hv) )
        error('relhfn2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ndims(tt) ~= +2 )
        error('relhfn2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(pp,2)~= +2 || size(tt,2) < +3 )
        error('relhfn2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(hv,2)~= +1 || size(hv,1) ~= size(pp,1) )
        error('relhfn2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nnod = size(pp,1) ;

%---------------------------------------------- basic checks
    if (min(min(tt(:,1:3))) < +1 || ...
            max(max(tt(:,1:3))) > nnod )
        error('relhfn2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end

%----------------------------------- compute rel. mesh-sizes
   [ee,tt] = tricon2(tt);
   
    evec = pp(ee(:,2),:)-pp(ee(:,1),:) ;
         
    elen = sqrt(sum(evec.^2,+2));

    hmid = hv(ee(:,2),:)+hv(ee(:,1),:) ;
    hmid = hmid * +0.50 ;
    hrel = elen ./ hmid ;

end


