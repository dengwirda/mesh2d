function [okay] = meshhas(mesh,varargin)
%MESHHAS helper routine to safely query a MESH structure.

%---------------------------------------------------------------------
%   Darren Engwirda
%   github.com/dengwirda/jigsaw-matlab
%   22-Mar-2016
%   d_engwirda@outlook.com
%---------------------------------------------------------------------
%

    okay = false; base = ''; item = '';

    if (nargin >= +2), base = varargin{1}; end
    if (nargin >= +3), item = varargin{2}; end
    
    if (nargin <= +1 || nargin >= +4)
        error('Incorrect number of arguments!') ;
    end
   
    if (~isstruct(mesh))
        error('MESH must be a valid struct.');
    end
    if (~isempty(base) && ~ischar(base))
        error('BASE must be a valid string!'); 
    end
    if (~isempty(item) && ~ischar(item))
        error('ITEM must be a valid string!');
    end
    
    if (isempty(item))
%-- default ITEM kinds given BASE types
    switch (base)
        case 'point', item = 'coord';
        case 'edge2', item = 'index';
        case 'tria3', item = 'index';
        case 'quad4', item = 'index';
        case 'tria4', item = 'index';
        case 'hexa8', item = 'index';
        case 'wedg6', item = 'index';
        case 'pyra5', item = 'index';
    end
    end
 
%-- check whether MESH.BASE.ITEM exists
    okay = isfield(mesh,base) && ...
        isfield(mesh.(base),item) && ...
            ~isempty(mesh.(base).(item)) ;
    
end


