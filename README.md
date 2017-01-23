# `MESH2D: A MATLAB-based mesh generator`

`MESH2D` is a simple `MATLAB`-based Delaunay mesh generator for two-dimensional geometries. It is designed to generate high-quality constrained Delaunay triangulations for general polygonal regions in the plane. `MESH2D` provides simple and yet effective implementations of "Delaunay-refinement" and "Frontal-Delaunay" triangulation techniques, in additional to "hill-climbing" type mesh-optimisation. 

# `Starting Out`

After downloading and unzipping the current <a href="https://github.com/dengwirda/mesh2d/archive/master.zip">repository</a>, navigate to the installation directory within <a href="http://www.mathworks.com">`MATLAB`</a> / <a href="https://www.gnu.org/software/octave">`OCTAVE`</a> and run the set of examples contained in `tridemo.m`:
````
tridemo(1); % investigate the impact of the "radius-edge" threshold.
tridemo(2); % Frontal-Delaunay vs. Delaunay-refinement refinement.
tridemo(3); % explore impact of user-defined mesh-size constraints.
tridemo(4); % explore impact of "hill-climbing" mesh optimisations.
tridemo(5); % assemble triangulations for "multi-part" geometries.
tridemo(6); % larger-scale problem, mesh refinement + optimisation. 
tridemo(7); % medium-scale problem, mesh refinement + optimisation. 
````
# `Attribution!`

If you make use of `MESH2D` please reference appropriately! `MESH2D` is designed to provide a simple and easy-to-understand implementation of Delaunay-based mesh-generation techniques. For a much more advanced, and fully three-dimensional mesh-generation library, see the <a href="https://github.com/dengwirda/jigsaw-matlab/">`JIGSAW`</a> package.

`[1]` - Darren Engwirda, <a href="http://hdl.handle.net/2123/13148">Locally-optimal Delaunay-refinement and optimisation-based mesh generation</a>, Ph.D. Thesis, School of Mathematics and Statistics, The University of Sydney, September 2014.
