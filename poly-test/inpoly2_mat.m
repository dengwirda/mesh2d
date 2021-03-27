function [stat,bnds] = ...
    inpoly2_mat(vert,node,edge,fTOL,lbar)
%INPOLY2_MAT the local m-code version of the crossing-number
%test. Loop over edges; do a binary-search for the first ve-
%rtex that intersects with the edge y-range; do crossing-nu-
%mber comparisons; break when the local y-range is exceeded.

%   Darren Engwirda : 2017 --
%   Email           : d.engwirda@gmail.com
%   Last updated    : 19/12/2020

    feps = fTOL * lbar ^ +1 ;
    veps = fTOL * lbar ^ +1 ;

    nvrt = size (vert,1) ;
    nnod = size (node,1) ;
    nedg = size (edge,1) ;

    stat = false(nvrt,1) ;
    bnds = false(nvrt,1) ;
    
%----------------------------------- loop over polygon edges
    for epos = +1 : size(edge,1)
    
        inod = edge(epos,1) ;
        jnod = edge(epos,2) ;

    %------------------------------- calc. edge bounding-box
        yone = node(inod,2) ;
        ytwo = node(jnod,2) ;
        xone = node(inod,1) ;
        xtwo = node(jnod,1) ;
        
        xmin = min(xone,xtwo) ;
        xmax = max(xone,xtwo) ;
        
        xmin = xmin - veps;
        xmax = xmax + veps;
        ymin = yone - veps;
        ymax = ytwo + veps;
        
        ydel = ytwo - yone;
        xdel = xtwo - xone;

        edel = abs(xdel) + ydel ;
        
    %------------------------------- find top VERT(:,2)<YONE
        ilow = +1; iupp = nvrt;
        
        while (ilow < iupp - 1)    % binary search    
            imid = ilow ...
            + floor((iupp-ilow) / 2);
            
            if (vert(imid,2) < ymin)
                ilow = imid ;
            else
                iupp = imid ;
            end           
        end
        
        if (vert(ilow,2) >= ymin)
            ilow = ilow - 1 ;
        end

    %------------------------------- calc. edge-intersection
        for jpos = ilow+1 : nvrt
       
            if (bnds(jpos)), continue ; end
        
            xpos = vert(jpos,1) ;
            ypos = vert(jpos,2) ;
            
            if (ypos <= ymax)
                if (xpos >= xmin)
                    if (xpos <= xmax)
            
                %------------------- compute crossing number    
                    mul1 = ...
                    ydel * (xpos - xone) ;
                    mul2 = ...
                    xdel * (ypos - yone) ;
                    
                    if ((feps * edel) >= ...
                        abs(mul2 - mul1) )
                    
                %------------------- BNDS -- approx. on edge
                        bnds(jpos)= true ;
                        stat(jpos)= true ;
                    
                    elseif (ypos == yone ...
                        &&  xpos == xone )
                        
                %------------------- BNDS -- match about ONE
                        bnds(jpos)= true ;
                        stat(jpos)= true ;
                        
                    elseif (ypos == ytwo ...
                        &&  xpos == xtwo )
                    
                %------------------- BNDS -- match about TWO
                        bnds(jpos)= true ;
                        stat(jpos)= true ;
                    
                    elseif (mul1 < mul2)
                    
                    if (ypos >= yone ...
                        && ypos <  ytwo)
                        
                %------------------- advance crossing number
                        stat(jpos) = ...
                            ~stat (jpos) ;
                    
                    end
                    
                    end
                
                    end
                else
                
                    if (ypos >= yone ...
                        && ypos <  ytwo)
                        
                %------------------- advance crossing number
                        stat(jpos) = ...
                            ~stat (jpos) ;
                    
                    end
                
                end
            else
            
                break ;            % done -- due to the sort
            
            end
                    
        end

    end

end



