function [ipos] = lsearch(test,tval)
%LSEARCH a binary search for the index TEST(IPOS) < TVAL. 
%   [IPOS] = LSEARCH(TEST,TVAL) returns the upper-most index
%   in the vector TEST, such that TEST(IPOS) < TVAL. A bina-
%   ry-search type method is employed, requiring approximat-
%   ely O(LOG(N)) per query. The vector TEST must be sorted
%   in ascending order on input.

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 16/01/2017

    ilow = 1; iupp = numel(test) ;

%------------------------------------------- check endpoints
    if (test(ilow) >= tval)
        ipos = ilow; return ;
    end
    if (test(iupp) <  tval)
        ipos = iupp; return ;
    end
    
%------------------------------------------- inner recursion
    while (ilow < iupp-1)
        imid = ilow ...
        + floor((iupp-ilow)/2) ;
        if (test(imid) < tval)
            ilow = imid ;
        else
            iupp = imid ;
        end
    end
    
    ipos =  ilow ;

end



