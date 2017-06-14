function [ok,pp,qq,tp,tq] = linenear(pa,pb,pc,pd)
%LINENEAR calc. the nearest points on line segments embedded 
%in d-dimensions.

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 13/06/2017
%-----------------------------------------------------------

    m1 = (pa+pb) * +.5 ;
    D1 = (pb-pa) * +.5 ;
    
    m2 = (pc+pd) * +.5 ;
    D2 = (pd-pc) * +.5 ;

    r1 = sum(m2.*D1,2) ...
       - sum(m1.*D1,2) ;
    r2 = sum(m1.*D2,2) ...
       - sum(m2.*D2,2) ;

    A1 = sum(D1.*D1,2) ;
    A2 =-sum(D1.*D2,2) ;
    A3 =-sum(D1.*D2,2) ;
    A4 = sum(D2.*D2,2) ;

    dd = A1.*A4 - A2.*A3 ;

    tp = A4.*r1 - A2.*r2 ;
    tq =-A3.*r1 + A1.*r2 ;

    rt = max(abs(tp),abs(tq)) ;
    rt = rt * eps ^ .8 ;
    
    ok = abs(dd) > +rt ;
    
    tp(~ok) = +0. ; 
    tq(~ok) = +0. ;
    
    tp(ok) = tp(ok) ./ dd(ok) ;
    tq(ok) = tq(ok) ./ dd(ok) ;
    
    tp = max(min(tp,+1.),-1.) ;
    tq = max(min(tq,+1.),-1.) ;
    
    pp = m1 + tp .* D1 ;
    qq = m2 + tq .* D2 ;
    
end


