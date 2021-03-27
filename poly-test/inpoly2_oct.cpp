
// a fast pre-sorted variant of the crossing-number test for 
// INPOLY2.m 

//----------------------------------------------------------
//  Darren Engwirda : 2017 --
//  Email           : d.engwirda@gmail.com
//  Last updated    : 19/12/2020
//----------------------------------------------------------

#include <octave/oct.h>

DEFUN_DLD (inpoly2_oct, args, ,
  "-- INPOLY2_OCT(TEST,NODE,EDGE,FEPS); \n"
  "-- \n"
  "-- INPOLY2-OCT: low-level routine called by INPOLY2 to \n"
  "-- compute point-in-polygon queries. \n" )
{    
    octave_value_list rval;
    
    int const nargin = args.length () ;
    if (nargin != +5)
    {
        print_usage () ;
        return rval ;
    }
    
    Matrix const vert (
        args(0).matrix_value ()) ;
    Matrix const node (
        args(1).matrix_value ()) ;
    Matrix const edge (
        args(2).matrix_value ()) ;
    double const fTOL (
        args(3).double_value ()) ;
    double const lbar (
        args(4).double_value ()) ;
    
    if (error_state) return rval ;

    octave_idx_type const nvrt 
        = vert.rows () ;
    octave_idx_type const nnod 
        = node.rows () ;
    octave_idx_type const nedg 
        = edge.rows () ;

//---------------------------------- init. crossing no. bool
    boolMatrix stat(nvrt, 1) ;
    boolMatrix bnds(nvrt, 1) ;
    
    octave_idx_type vpos ;
    for (vpos = +0; vpos != nvrt; ++vpos)
    {
        stat(vpos) = false ;
        bnds(vpos) = false ;
    }

//---------------------------------- loop over polygon edges
    double const veps = fTOL * lbar ;
    double const feps = fTOL * lbar ;
    
    octave_idx_type epos ;
    for (epos = +0; epos != nedg; ++epos)
    {
        octave_idx_type const inod 
            = edge(epos,0) - 1 ;
        octave_idx_type const jnod 
            = edge(epos,1) - 1 ;

    //------------------------------ calc. edge bounding-box
        double yone = node (inod,1) ;
        double ytwo = node (jnod,1) ;
        double xone = node (inod,0) ;
        double xtwo = node (jnod,0) ;
        
        double xmin = xone < xtwo 
                    ? xone : xtwo ;
        double xmax = xone < xtwo 
                    ? xtwo : xone ;               
     
        xmin-= veps ;
        xmax+= veps ;
        
        double ymin = yone - veps ;
        double ymax = ytwo + veps ;
        
        double ydel = ytwo - yone ;
        double xdel = xtwo - xone ;
            
        double edel = 
            std::abs(xdel) + ydel ;

    //------------------------------ find top VERT(:,2)<YONE
        octave_idx_type ilow = +0 ; 
        octave_idx_type iupp = 
            nvrt - 1;
        octave_idx_type imid = +0 ;
        
        while (ilow < iupp - 1)
        {
            imid = ilow 
                + (iupp - ilow) / 2;
            
            if (vert(imid,1) < ymin)
            {
                ilow = imid ;
            }
            else
            {
                iupp = imid ;
            }
        }
        
        {
            if (vert(ilow,1) >=ymin)
            {
                ilow -= +1 ; 
            }
        }
        
    //------------------------------ calc. edge-intersection
        octave_idx_type vpos = ilow+1 ;
        for ( ; vpos != nvrt; ++vpos)
        {
            if (bnds(vpos)) continue;
        
            double xpos = vert(vpos,0);
            double ypos = vert(vpos,1);
            
            if (ypos <= ymax)
            {
                if (xpos >= xmin)
                {
                    if (xpos <= xmax)
                    {             // compute crossing number
                    double mul1 =
                    ydel * (xpos - xone) ;
                    double mul2 =
                    xdel * (ypos - yone) ;
                    
                    if ((feps * edel) >=
                    std::abs(mul2 - mul1) )
                    {             // BNDS -- approx. on edge
                        bnds(vpos)= true ;
                        stat(vpos)= true ;
                    }
                    else
                    if (ypos == yone &&
                        xpos == xone )
                    {             // BNDS -- match about ONE
                        bnds(vpos)= true ;
                        stat(vpos)= true ;
                    }   
                    else
                    if (ypos == ytwo &&
                        xpos == xtwo )
                    {             // BNDS -- match about TWO
                        bnds(vpos)= true ;
                        stat(vpos)= true ;
                    }
                    else
                    if (mul1 <  mul2 )
                    {
                    if (ypos >= yone
                        && ypos <  ytwo)
                    {             // advance crossing number
                        stat(vpos) =
                           ! stat (vpos) ;
                    }
                    }
                
                    }
                }
                else
                {
                    if (ypos >= yone
                        && ypos <  ytwo)
                    {             // advance crossing number
                        stat(vpos) =
                           ! stat (vpos) ;
                    }
                }
            }
            else
            {
                break ;           // done -- due to the sort
            }
        }

    }

    rval(0) = stat;
    rval(1) = bnds;
    
    return ( rval ) ;
}



