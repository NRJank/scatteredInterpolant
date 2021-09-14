########################################################################
##
## Copyright (C) 2007-2021 The Octave Project Developers
##
## See the file COPYRIGHT.md in the top-level directory of this
## distribution or <https://octave.org/copyright/>.
##
## This file is part of Octave.
##
## Octave is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <https://www.gnu.org/licenses/>.
##
########################################################################

## -*- texinfo -*-
## @deftypefn  {} {@var{T} =} delaunayn (@var{pts})
## @deftypefnx {} {@var{T} =} delaunayn (@var{pts}, @var{options})
## Compute the Delaunay triangulation for an N-dimensional set of points.
##
## The Delaunay triangulation is a tessellation of the convex hull of a set of
## points such that no N-sphere defined by the N-triangles contains any other
## points from the set.
##
## The input matrix @var{pts} of size [n, dim] contains n points in a space of
## dimension dim.  The return matrix @var{T} has size [m, dim+1].  Each row of
## @var{T} contains a set of indices back into the original set of points
## @var{pts} which describes a simplex of dimension dim.  For example, a 2-D
## simplex is a triangle and 3-D simplex is a tetrahedron.
##
## An optional second argument, which must be a string or cell array of
## strings, contains options passed to the underlying qhull command.  See the
## documentation for the Qhull library for details
## @url{http://www.qhull.org/html/qh-quick.htm#options}.
## The default options depend on the dimension of the input:
##
## @itemize
## @item 2-D and 3-D: @var{options} = @code{@{"Qt", "Qbb", "Qc"@}}
##
## @item 4-D and higher: @var{options} = @code{@{"Qt", "Qbb", "Qc", "Qx"@}}
## @end itemize
##
## If Qhull fails for 2-D input the triangulation is attempted again with
## the options @code{@{"Qt", "Qbb", "Qc", "Qz"@}} which may result in
## reduced accuracy.
##
## If @var{options} is not present or @code{[]} then the default arguments are
## used.  Otherwise, @var{options} replaces the default argument list.
## To append user options to the defaults it is necessary to repeat the
## default arguments in @var{options}.  Use a null string to pass no arguments.
##
## @seealso{delaunay, convhulln, voronoin, trimesh, tetramesh}
## @end deftypefn

function T = delaunayn (pts, varargin)

  if (nargin < 1)
    print_usage ();
  endif
%keyboard
  if (isempty (varargin) || isempty (varargin{1}))
    try
      T = __delaunayn__ (pts);
    catch
      if (columns (pts) <= 2)
        T = __delaunayn__ (pts, "Qt Qbb Qc Qz");
      endif
    end_try_catch
  else
    T = __delaunayn__ (pts, varargin{:});
  endif

  ##Begin check for and removal of trivial simplices

  if (! isequal (T, 0)) ## skip trivial simplex check if no simplexes

    if (isa (pts, "single"))
      tol = 1e3 * eps ("single");
    else
      tol = 1e3 * eps;
    endif

    ## Try to remove the zero volume simplices.  The volume of the i-th simplex is
    ## given by abs(det(pts(T(i,1:end-1),:)-pts(T(i,2:end),:)))/factorial(ndim+1)
    ## (reference http://en.wikipedia.org/wiki/Simplex).  Any simplex with a
    ## relative volume less than some arbitrary criteria is rejected.  The
    ## criteria we use is the volume of the simplex corresponding to an
    ## orthogonal simplex with (ndim-1) edge lengths equal to the edge lengths of
    ## the original simplex.  If the relative volume is smaller than 1e3*eps, the
    ## simplex is rejected.  Note division of the two volumes means that the
    ## factor factorial(ndim+1) is dropped from the volume calculations.

    [nt, nd] = size (T); ## nt = simplex count, nd = # of simplex points
    dim = nd - 1;

    ## calculate common origin edge vectors for each simplex (p2-p1, p3-p1, ...)
    ## stored in 3D array such that
    ## rows = nt simplexes, cols = coordinates, pages = simplex edges
    edg_vec =  permute(reshape(pts(T(:,2:nd),:).',[dim, nt, dim]),[2 1 3]) ...
                 - pts(T(:,1), :, ones (1,1,dim));

    ## Calculate simplex volumes according to dimensionality of problem
    if  (dim <= 6)
      ## faster vector-product code paths for up to 3D cases

      if (dim == 2)
        ## 2D Use simple component cross product to calculate 2D triangle area
        vol = abs (edg_vec(:,1,1) .* edg_vec(:,2,2)...
                 - edg_vec(:,1,2) .* edg_vec(:,2,1));
      else
        ## >=3D: Use scalar triple product to calculate 3D tetrahedron volume, and
        ## explicit laplace expansion for higher dimension determinants 
        vol = abs (laplace_det (dim, edg_vec));
      endif

    else
      ## higher dimensions

      ## Place edge vectors into a block diagonal matrix
      eqs = sparse (dim * nt, dim * nt);
      eqs(logical (kron (speye (nt, nt), true (dim)))) = ...
            permute (edg_vec, [2,3,1])(:);
      ## Extract diagonal of LU factorization of that block diagonal matrix
      [~, eqs, p, ~] = lu (eqs, "vector");
      R = abs (diag (eqs));

      ## extract simplex volumes as product of diagonal elements of u
      ## preserving order relative to simplexes from delaunay triangulation
      [~, rev_sort]= sort (kron (1:nt, ones (1, dim))(p));
      vol = prod (reshape (R(rev_sort), dim, nt), 1).';

    endif

    ## Check for small volumes: compare simplex to volume with orthogonal edges
    idx = (vol ./ prod (sqrt (sumsq (edg_vec, 2)), 3)) < tol;

    ##Remove trivially small simplexes 
    T(idx,:) = [];
  endif
endfunction

function retval = laplace_det (dim, v)
  ## using scalar triple product for 3D determinant, and for >3D use laplace
  ## expansion recursively down to 3D.  
  ## input vector shape: rows - simplex, column - coordinate, page - edgevectors
  switch dim
    case 3
##      retval = dot (v(:,:,1), cross (v(:,:,2), v(:,:,3), 2), 2);
      retval = v(:,1,1) .* (v(:,2,2) .* v(:,3,3) - v(:,3,2) .* v(:,2,3))...
        - v(:,2,1) .* (v(:,1,2) .* v(:,3,3) - v(:,3,2) .* v(:,1,3))...
        + v(:,3,1) .* (v(:,1,2) .* v(:,2,3) - v(:,2,2) .* v(:,1,3));

    case 4
      retval = v(:,1,1) .* laplace_det (3, v(:, 2:4, 2:4))...
        - v(:,2,1) .* laplace_det (3, v(:, [1,3,4], 2:4))...
        + v(:,3,1) .* laplace_det (3, v(:, [1,2,4], 2:4))...
        - v(:,4,1) .* laplace_det (3, v(:, 1:3, 2:4));

    case 5
      retval = v(:,1,1) .* laplace_det (4, v(:, 2:5, 2:5))...
        - v(:,2,1) .* laplace_det (4, v(:, [1,3:5], 2:5))...
        + v(:,3,1) .* laplace_det (4, v(:, [1,2,4,5], 2:5))...
        - v(:,4,1) .* laplace_det (4, v(:, [1:3,5], 2:5))...
        + v(:,5,1) .* laplace_det (4, v(:, 1:4, 2:5));

    case 6
      retval = v(:,1,1) .* laplace_det (5, v(:,2:6, 2:6))...
        - v(:,2,1) .* laplace_det (5, v(:, [1,3:6], 2:6))...
        + v(:,3,1) .* laplace_det (5, v(:, [1,2,4:6], 2:6))...
        - v(:,4,1) .* laplace_det (5, v(:, [1:3,5,6], 2:6))...
        + v(:,5,1) .* laplace_det (5, v(:, [1:4,6], 2:6))...
        - v(:,6,1) .* laplace_det (5, v(:, 1:5, 2:6));

    case 7
      retval = v(:,1,1) .* laplace_det (6, v(:, 2:7, 2:7))...
        - v(:,2,1) .* laplace_det (6, v(:, [1,3:7], 2:7))...
        + v(:,3,1) .* laplace_det (6, v(:, [1,2,4:7], 2:7))...
        - v(:,4,1) .* laplace_det (6, v(:, [1:3,5:7], 2:7))...
        + v(:,5,1) .* laplace_det (6, v(:, [1:4,6,7], 2:7))...
        - v(:,6,1) .* laplace_det (6, v(:, [1:5,7], 2:7))...
        + v(:,7,1) .* laplace_det (6, v(:, 1:6, 2:7));

    otherwise
      error ("delaunayn: DIM too high for laplace expansion deternimant");
  endswitch

endfunction


%!testif HAVE_QHULL
%! x = [-1, 0; 0, 1; 1, 0; 0, -1; 0, 0];
%! assert (sortrows (sort (delaunayn (x), 2)), [1,2,5;1,4,5;2,3,5;3,4,5]);

## Test 3-D input
%!testif HAVE_QHULL
%! x = [-1, -1, 1, 0, -1]; y = [-1, 1, 1, 0, -1]; z = [0, 0, 0, 1, 1];
%! assert (sortrows (sort (delaunayn ([x(:) y(:) z(:)]), 2)), [1,2,3,4;1,2,4,5]);

## FIXME: Need tests for delaunayn

## Input validation tests
%!error <Invalid call> delaunayn ()
