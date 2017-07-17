function compile
%COMPILE compile any MESH2D dependencies into MEX//OCT files
%for MATLAB//OCTAVE.
%
%   See also TRIDEMO, REFINE2, SMOOTH2

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 17/07/2017
%-----------------------------------------------------------

    if (exist('OCTAVE_VERSION','builtin')  > 0)

    fprintf(1,'\n');
  
    fprintf(1,'OCTAVE detected: Compiling OCT files.\n');
    
    fprintf(1,'\n');
    
    fprintf(1,'Compiling INPOLY2...\n');
    mkoctfile('inpoly2_oct.cpp');

    %%!! add additional oct-file builds here... 


    fprintf(1,'\n');

    fprintf(1,'Compiling MESH2D for OCTAVE: DONE.\n');

    fprintf(1,'\n');

    else

    fprintf(1,'\n');
    
    fprintf(1,'MATLAB detected: nothing to compile!!\n');
    
    fprintf(1,'\n');
    
    fprintf(1,'Compiling MESH2D for MATLAB: DONE.\n');

    fprintf(1,'\n');

    end

end



