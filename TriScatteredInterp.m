## Copyright (C) 2013 Ben Kurtz
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {} {@var{f} =} TriScatteredInterp
## @deftypefnx {} {@var{f} =} TriScatteredInterp (@var{x}, @var{y}, @var{q})
## @deftypefnx {} {@var{f} =} TriScatteredInterp (@var{x}, @var{y}, @var{z}, @var{q})
## @deftypefnx {} {@var{f} =} TriScatteredInterp (@var{P}, @var{q})
## @deftypefnx {} {@var{f} =} TriScatteredInterp (@dots{}, @var{method})
##
## Generate a 2D or 3D interpolating function defined by an input set of
## scattered datapoints @var{x}, @var{y}, @var{z}, and function values, @var{q}
## at those points, all specified as equal length vectors. The input points may
## be specified as a single matrix, @var{P} where each row of @var{P} contains
## the coordinates of a single point.
##
## Interpolation is based on Delaunay triangulation and can be controlled by
## supplying options for @var{method}, which can take the values of:
##
## @table @code
## @item nearest
## Discontinuous interpolation assigning value of nearest sample data point.
##
## @item linear
## (Default) Linear interpolation between nearest grid points according to local
## triangulation. 
##
## @item natural
## Natural neighbor interpoltaion 
## @end table
##
## @code{TriScatteredInterp} provides no extrapolation and will retrun NaN for
## any query points outside the convex hull defined by the sample points used to
## define @var{f}.
##
## Matlab Compatibility Note:  @code{TriScatteredInterp} is impleneted as a
## wrapper for @code{scatteredInterpolant}. As such it will return a
## @code{scatteredInterpolant} class object instead of a
## @code{TriScatteredInterp} object and and does not currently accept
## DelaunayTri objects as input.
##
## @seealso{scatteredInterpolant, griddata}
## @end deftypefn

function F = TriScatteredInterp (varargin)
  if ((nargin == 1) || (nargin > 5))
    print_usage ();
  endif
  if (nargin == 0)
    ## empty input is valid, but need to modiy extrap method to none.
    F = scatteredInterpolant ();
    F.ExtrapolationMethod = "none";
  else
    ## need to specify none for extrap method, so must pass a method. therefore
    ## need to know if one was specified, and pass the default if it was not.
    ## Only char inputs should be a single method input as last input.
    method_inputs = cellfun(@ischar, varargin);
    if any (method_inputs)
      if (sum (method_inputs) > 1)
        error ("TriScatteredInterp: only a single method input may be used.")
      elseif (method_inputs(end)~= 1)
        error ("TriScatteredInterp: method input must be last.")
      endif
      method = varargin{end};
      if ~any (strcmp (method, {"linear", "natural", "nearest"}))
        error ("TriScatteredInterp: unknown method %s", method);
      endif
    else
      method = "linear";
    endif

    ## A Matlab compatible TriScatteredInterpolant requires column vector
    ## inputs, but as Octave's statteredInterpolant allows arbitrary orientation
    ## there's no need to check and restrict it here.
    F = scatteredInterpolant (varargin{~method_inputs}, method, "none");
  endif
endfunction
