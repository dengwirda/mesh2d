function [bb] = tribal2(pp,tt)
%TRIBAL2 compute the circumballs associated with a 2-simplex
%triangulation embedded in R^2.
%   [BB] = TRIBAL2(PP,TT) returns the circumscribing balls
%   associated with the triangles in [PP,TT], such that BB = 
%   [XC,YC,RC.^2].

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 14/01/2017

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
    if (size(pp,2)~= +2 || size(tt,2)<= +3)
        error('tribal2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    bb = zeros(size(tt,1),3);
    rv = zeros(size(tt,1),2);
    
%------------------------------------------------ tria edges
    ab = pp(tt(:,2),:)-pp(tt(:,1),:);
    ac = pp(tt(:,3),:)-pp(tt(:,1),:);
    
%------------------------------------------------ lhs matrix     
    aa = zeros(2,2,size(tt,1));
    aa(1,:,:) = +2.*ab';
    aa(2,:,:) = +2.*ac';
    
%------------------------------------------------ rhs vector
    rv(:,1) = sum(ab.*ab,2);
    rv(:,2) = sum(ac.*ac,2);
    rv = rv';
    
%------------------------------------------------ solve sys.
    bb(:,1:2) = blocksolve(aa,rv)';
    bb(:,3) = sum(bb(:,1:2).^2,2) ;
    bb(:,1:2) = ...
        pp(tt(:,1),:) + bb(:,1:2) ;

end



