function [bb] = tribal2(pp,tt)
%TRIBAL2 compute the circumballs associated with a 2-simplex
%triangulation embedded in R^2.
%   [BB] = TRIBAL2(PP,TT) returns the circumscribing balls
%   associated with the triangles in [PP,TT], such that BB = 
%   [XC,YC,RC.^2].

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 29/04/2017

%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ...
        ~isnumeric(tt) )
        error('tribal2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ndims(tt) ~= +2)
        error('tribal2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(pp,2)~= +2 || size(tt,2) < +3)
        error('tribal2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

%------------------------------------------------ lhs matrix     
    ab = pp(tt(:,2),:)-pp(tt(:,1),:) ;
    ac = pp(tt(:,3),:)-pp(tt(:,1),:) ;
   
%------------------------------------------------ rhs vector    
    rv11 = sum(ab.*ab,2) ;
    rv22 = sum(ac.*ac,2) ;
    
%------------------------------------------------ solve sys.
    dd      = ab(:,1) .* ac(:,2) - ...
              ab(:,2) .* ac(:,1) ;
              
    bb = zeros(size(tt,1),3) ;
    bb(:,1) = (ac(:,2) .* rv11 - ...
               ab(:,2) .* rv22 ) ...
            ./ dd * +.5 ;
            
    bb(:,2) = (ab(:,1) .* rv22 - ...
               ac(:,1) .* rv11 ) ...
            ./ dd * +.5 ;
      
    bb(:,3) = sum(bb(:,1:2).^2,2) ;
    bb(:,1:2) = ...
        pp(tt(:,1),:) + bb(:,1:2) ;

end



