## `INPOLY: Fast point(s)-in-polygon queries in MATLAB`

A fast 'point(s)-in-polygon' routine for <a href="http://www.mathworks.com">`MATLAB`</a> / <a href="https://www.gnu.org/software/octave">`OCTAVE`</a>.

`INPOLY` returns the "inside/outside" status for a set of vertices `VERT` and a general polygon (`PSLG`) embedded in the two-dimensional plane. General non-convex and multiply-connected polygonal regions can be handled. `INPOLY` is intended as a (very) fast replacement for `MATLAB`'s default `INPOLYGON` routine.

<p align="center">
  <img src = "../master/test-data/query.png">
</p>

`INPOLY` is based on a 'crossing-number' test, counting the number of times a line extending from each point past the right-most region of the polygon intersects with the polygonal boundary. Points with odd counts are 'inside'. A simple implementation requires that each edge intersection be checked for each point, leading to (slow) `O(N*M)` overall complexity.

This implementation seeks to improve these bounds. Query points are sorted by `y-value` and candidate intersection sets are determined via binary-search. Given a configuration with `N` test points, `M` edges and an average point-edge 'overlap' of `H`, the overall complexity scales like `O(M*H + M*LOG(N) + N*LOG(N))`, where `O(N*LOG(N))` operations are required for the initial sorting, `O(M*LOG(N))` operations are required for the set of binary-searches, and `O(M*H)` operations are required for the actual intersection tests. `H` is typically small on average, such that `H << N`. Overall, this leads to fast `O((N+M)*LOG(N))` complexity for average cases.

### `Quickstart`

After downloading and unzipping the current <a href="https://github.com/dengwirda/inpoly/archive/master.zip">repository</a>, navigate to the installation directory within <a href="http://www.mathworks.com">`MATLAB`</a> / <a href="https://www.gnu.org/software/octave">`OCTAVE`</a> and run the set of examples contained in `polydemo.m`:
````
polydemo(1); % a simple example
polydemo(2); % multiply-connected domains
polydemo(3); % speed comparison
````
For good performance in `OCTAVE`, `INPOLY` can be compiled from `C++` as an `*.oct` file. See documentation on `MKOCTFILE` for additional information. Typically, `MATLAB`'s in-built `JIT`-acceleration leads to good performance by default.

### `License Terms`

This program may be freely redistributed under the condition that the copyright notices (including this entire header) are not removed, and no compensation is received through use of the software.  Private, research, and institutional use is free.  You may distribute modified versions of this code `UNDER THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE MODIFICATIONS`. Distribution of this code as part of a commercial system is permissible `ONLY BY DIRECT ARRANGEMENT WITH THE AUTHOR`. (If you are not directly supplying this code to a customer, and you are instead telling them how they can obtain it for free, then you are not required to make any arrangement with me.) 

`DISCLAIMER`:  Neither I nor the University of Sydney warrant this code in any way whatsoever. This code is provided "as-is" to be used at your own risk.

### `References`

`[1]` - J. Kepner, D. Engwirda, V. Gadepally, C. Hill, T. Kraska, M. Jones, A. Kipf, L. Milechin, N. Vembar: <a href="https://arxiv.org/abs/2005.03156">Fast Mapping onto Census Blocks</a>, IEEE HPEC, 2020.
