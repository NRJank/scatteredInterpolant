########################################################################
##
## Copyright (C) 2011-2021 The Octave Project Developers
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

classdef scatteredInterpolant

  ## -*- texinfo -*-
  ## @deftypefn {} {@var{f} =} scatteredInterpolant
  ## @deftypefnx {} {@var{f} =} scatteredInterpolant (@var{x}, @var{y}, @var{q})
  ## @deftypefnx {} {@var{f} =} scatteredInterpolant (@var{x}, @var{y}, @var{z}, @var{q})
  ## @deftypefnx {} {@var{f} =} scatteredInterpolant (@var{P}, @var{q})
  ## @deftypefnx {} {@var{f} =} scatteredInterpolant (@dots{}, @var{method})
  ## @deftypefnx {} {@var{f} =} scatteredInterpolant (@dots{}, @var{method}, @var{extrapolationmethod})
  ##
  ## Generate a 2D or 3D interpolating function @var{f} defined by an input set
  ## of scattered datapoints @var{x}, @var{y}, @var{z}, and function values,
  ## @var{q} at those points.
  ##
  ## @var{x}, @var{y}, @var{z} must be vectors of the same length as @var{q}.
  ## The input points may also be specified as a single matrix, @var{P} where
  ## each row of @var{P} contains the coordinates of a single point
  ## corresponding to function values in @var{q}.
  ##
  ## The returned function can then be evaluated as any of:
  ##
  ## @example
  ## @group
  ## @var{Si} = @var{f} (@var{xi}, @var{yi})
  ## @var{Si} = @var{f} (@var{xi}, @var{yi}, @var{zi})
  ## @var{Si} = @var{f} (@var{Pi})
  ## @var{Si} = @var{f} (@{@var{Xg}, @var{Yg}@})
  ## @var{Si} = @var{f} (@{@var{Xg}, @var{Yg}, @var{Zg}@})
  ## @end group
  ## @end example
  ##
  ## where @var{xi}, @var{yi}, (and @var{zi} for 3D points) are equal sized
  ## arrays containing coordinates of every query point, and @var{Pi} refers to
  ## a single matrix of query points, with each row consisting of the @var{x},
  ## @var{y} (and @var{z}) coordinates of a single point.  @{@var{Xg}, @var{Yg},
  ## @var{Zg}@} refers to a cell array of grid vectors defining each axis of
  ## a grid space to be interpolated.
  ##
  ## The output Si will contain the interpolated values with the same array
  ## shape as the inputs, except for grid vector inputs, where the output will
  ## be the shape of the defined grid.  Note that the grid vector ordering uses
  ## @code{ndgrid}, where the first vector is expanded along dim-1, and not
  ## @code{meshgrid}, where the first vector is expanded along the 'horizontal'
  ## or dim-2 dimension.  @xref{Three-Dimensional Plots} for more detail.
  ##
  ## Interpolation is based on Delaunay triangulation and can be controlled by
  ## supplying options for @var{method}, which can take the values of:
  ##
  ## @table @code
  ## @item nearest
  ## Discontinuous interpolation assigning value of nearest sample data point.
  ##
  ## @item linear
  ## (Default) Linear interpolation between nearest grid points according to
  ## local triangulation.
  ##
  ## @item natural
  ## Natural neighbor interpoltaion
  ## @end table
  ##
  ## Extrapolation behavior can be defined by providing by supplying an optional
  ## @var{ExtrapolationMethod} option.  Valid extrapolation methods are:
  ##
  ## @table @code
  ## @item none
  ## No extrapolation is performed and any points outside of the convex hull of
  ## the points used to define @var{f} will return a value of NaN.
  ##
  ## @item linear
  ## Linear extrapolation will be performed based on the gradient at nearby
  ## triangulation boundaries.  This is the default behavior when interpolation
  ## method is either @option{linear} or @option{natural}.
  ##
  ## @item nearest
  ## Nearest neighbor extrapolation will return the value of the nearest
  ## neigboring point on the boundary.  This is the default when interpolation
  ## method is set to @option{nearest}.
  ## @end table
  ##
  ## @sc{Matlab} compatibility note:  Octave's @code{scatteredInterpolant} class
  ## uses @code{delaunayn} to calculate the triangulation for interpolation,
  ## whereas @sc{Matlab} likely uses the newer delaunayTriangulation objects.
  ## Certain point inputs fail for delaunayn in both programs that can be
  ## handled by delaunayTriangulation (for example the 8 corner points of a
  ## unit cube). Additionally, the  existance of multiple valid Delaunay
  ## triangulations for a given point set can produce valid, but differing
  ## differing interpolation results between the two programs.  This will be an
  ## unavoidable scatteredInterpolant output compatibilty limitation at least
  ## until a compatible delaunayTriangulation is implemented.
  ##
  ## @seealso{TriScatteredInterp, griddata, delaunay, delaunayn}
  ## @end deftypefn

  properties (Access = public)
    Points = []; # 2D or 3D array of sampling Points used to define interpolation
    Values = []; # function values at sampling Points
    Method = "linear"; # interpolation method: linear, nearest, or natural
    ExtrapolationMethod = "linear"; # extrapolation method: linear, nearest, or none
%%  endproperties
%%
%%  properties (Access = private, Hidden = true)
    dimension = 0 #values - 0 (empty), 2 (2D), or 3 (3D), from columns (Points)
    tri = [] #stored delaunay triangulation

    ## state variables for warnings/errors. start all as 'valid' for empty
    ## constructor
    valid_tri = true; #is current triangulation valid
    enough_points = true; # are there enough points (>3 for 2D, >4 for 3D)
    valid_points_vals = true; # do points and values match in number/dimensions
  endproperties

  methods (Access = public)

    function this = scatteredInterpolant (varargin) ## constructor

      if (nargin == 0)
        ## do nothing, defaults already set, avoid other error checks.

      elseif (nargin == 1 || nargin > 6)
        ## if not empty, varargin must contain 2-6 elements
        print_usage ("scatteredInterpolant");

      elseif (isempty (varargin{1}) || isempty (varargin{2}))
        error ("scatteredInterpolant: input points and values cannot be empty");

      else
        ## must have at least two array inputs at the front. can have one or two
        ## char inputs, but must be last.
        ## array inputs must either be 3 or 4 equal length vectors, or one
        ## 2 or 3 column array and one vector, with number of rows in the array
        ## equal to vector length. The last numeric input is the value vector q.
        ## anything that follows must be 1 or 2 char inputs.
        ##
        ## Compatibility note - Matlab 2021a errors if not column vector inputs, but
        ## docs just say vectors.  Allowing vectors as a superset of compatible
        ## function and in case later behavior changes to allow more general inputs.

        numer_input_loc = cellfun (@isnumeric, varargin);
        numer_vector_loc = cellfun (@isvector, varargin) & numer_input_loc;
        char_input_loc = cellfun (@ischar, varargin);

        ## make sure no input that isn't numeric or a char
        if (! all (numer_input_loc + char_input_loc))
          error ("scatteredInterpolant: inputs must be numeric arrays or option strings ");
        endif

        ## make sure all numeric input comes before char input
        ## no char inputs produce empty comparison result which if counts as false
        if (find (numer_input_loc, 1, "last") >= find (char_input_loc, 1, "first"))
          print_usage ("scatteredInterpolant");
        endif

        num_input_count = sum (numer_input_loc);
        char_input_count = sum (char_input_loc);

        ## check that number of data and method inputs don't fall out of bounds
        if ((num_input_count < 2 || num_input_count > 4) || (char_input_count > 2))
          print_usage ("scatteredInterpolant");
        endif

        ## make sure last numeric input, should be q, is a vector
        if (!(numer_vector_loc(num_input_count) == 1))  #######################can this be !=
          error ("scatteredInterpolant: Value input must be a numeric vector.");
        endif

        ## make all vectors are column vectors for easier processing later.
        varargin(numer_vector_loc) = cellfun (@(x) x(:), ...
          varargin(numer_vector_loc), "UniformOutput", false);


        ## make sure all numeric inputs have equal number of rows as q
        numpoints = numel (varargin{num_input_count});
        if (any (cellfun (@rows, varargin(numer_input_loc)) != numpoints))
          error ("scatteredInterpolant: Point and Value inputs must have the same number of rows.");
        endif

        ## process numeric inputs, assign to this.Points and this.Values.
        ## check that numeric inputs are either 3-4 columns or
        ## first one is 2-3 wide array and second is a column

        switch num_input_count
          case 2
            if ((! ismatrix (varargin{1})) || (! any (size (varargin{1}, 2) == [2,3])))
              error ("scatteredInterpolant: Point input must be a 2 or 3 column array.");
            endif

            ## already verified second input is column vector numeric input.

            this.Points = varargin{1};
            this.Values = varargin{2};

            this = check_points (this);
            this = setTriangulation (this);

          case {3, 4}
            ## already verified numeric inputs are column vectors.
            this.Points = [varargin{1:num_input_count-1}];
            this.Values = varargin{num_input_count};

            this = check_points (this);
            this = setTriangulation (this);

          otherwise
            print_usage ("scatteredInterpolant");
        endswitch

        switch char_input_count
          case 0
            # Do nothing, Method defaults already set when this created
            # this.Method = "linear"; this.ExtrapolationMethod = "linear";

          case 1
            this.Method = tolower (varargin{num_input_count + 1});
            switch this.Method
              case {"linear", "natural"}
                # Do nothing. extrap method default 'linear' set when this created.
              case {"nearest"}
                this.ExtrapolationMethod = "nearest";
            otherwise
                error ("scatteredInterpolant: Invalid METHOD '%s'", this.Method);
            endswitch

          case 2
            this.Method = tolower (varargin{num_input_count + 1});
            this.ExtrapolationMethod = tolower (varargin{num_input_count + 2});

            #verify both methods are valid.
            if (! any (strcmp (this.Method, {"linear", "nearest", "natural"})))
              error ("scatteredInterpolant: Invalid METHOD '%s'", this.Method);
            endif

            if (strcmp (this.ExtrapolationMethod, "boundary"))
              error ("scatteredInterpolant: EXTRAPOLATIONMETHOD 'boundary' not yet implemented");

            elseif (! any (strcmp (this.ExtrapolationMethod, {"linear", "nearest", "none"})))
              error ("scatteredInterpolant: Invalid EXTRAPOLATIONMETHOD '%s'", ...
                       this.ExtrapolationMethod);
            endif

          otherwise
            print_usage ("scatteredInterpolant");
        endswitch

      endif
    endfunction

    function v = subsref (this, S)
      ## subsref either returns a property value, or if () does the actual interpolation
      ## performing interpolation should return errors or warnings if in a bad state


      #issue warning if Point status flags set
      this.dimension = columns (this.Points); ## 0 (empty), 2 (2D), 3(3D)

      if (! this.valid_tri)
        if (! this.enough_points)
          warning("scatteredInterpolant: not enough points to create a triangulation\n");
        else
          warning("scatteredInterpolant: cannot calculate triangulation from given points\n");
        endif
      endif

      for S_idx = 1: numel(S)
        if (S_idx == 1)

          switch S(1).type(1)
            case "("

              if (isempty (this.Points) || (! this.valid_tri)) ...
                   || (any (cellfun (@isempty, S(1).subs))) ...
                     || isempty (S(1).subs)
                v = [];

              elseif (! this.valid_points_vals)
                ## unequal number of Points and values, throw error
                error ("scatteredInterpolant: unequal number of points and values, cannot interpolate")

              else

                ## Query points input validation

                ## interp output will be a column vector.  If input has any other
                ## size, flag must be changed to true and sz_output must be set
                ## for final reshape
                reshape_flag = false;

                num_query_elements = numel (S(1).subs);
                switch num_query_elements
                  case 1
                    if (isnumeric (S(1).subs{1}))
                      ## single numeric input array. must be 2D array.
                      ## columns must match dim.
                      ## can allow a col vector to be accepted, handled as a
                      ## single point by switching to row vector
                      qpts = S(1).subs{1};

                      if (ndims (qpts) > 2)
                        error ("scatteredInterpolant: query points must be 2D vectors or arrays");
                      endif

                      ## if vector, ensure row vector
                      if (isvector (qpts))
                        qpts = qpts(:).';
                      endif

                      #check for correct dimensionality
                      if (columns (qpts) != this.dimension)
                        error ("scatteredInterpolant: query points dimension must match interpolant");
                      endif

                    elseif (iscell (S(1).subs{1}))
                      ## cell inputs must the same number of vectors as
                      ## dimension.  process as grid vectors to build query point
                      ## array.
##                      keyboard
                      if (! all (cellfun (@isnumeric, S(1).subs{1})))
                        error ("scatteredInterpolant: query grid vectors must be numeric");

                      elseif (numel (S(1).subs{1}) != this.dimension)
                        error ("scatteredInterpolant: query grid vector count must match interpolant dimension");

                      elseif (! all (cellfun (@isvector, S(1).subs{1})))
                        error ("scatteredInterpolant: query grid vectors must be row or column vectors");
                      endif


                      ## extract grid vectors and produce point array
                      ## vector orientation doesn't matter for ndgrid
                      switch this.dimension
                        case 2
                          [qpts_x, qpts_y] = ndgrid (S(1).subs{:}{:});
                          qpts = [qpts_x(:), qpts_y(:)];

                        case 3
                          [qpts_x, qpts_y, qpts_z] = ndgrid (S(1).subs{:}{:});
                          qpts = [qpts_x(:), qpts_y(:), qpts_z(:)];
                      endswitch
                      sz_output = size (qpts_x);
                      reshape_flag = true;

                    else
                      print_query_points_usage (this);
                    endif

                  case {2,3}

                    ## all query inputs need to be numeric vectors or arrays
                    if (! all (cellfun (@isnumeric, S(1).subs)))
                      print_query_points_usage (this);

                    elseif (num_query_elements != this.dimension)
                      error ("scatteredInterpolant: query points dimension must match interpolant");
                    endif

                    ## check for vectors to be equal length
                    if ((all (cellfun (@isvector, S(1).subs))) ...
                        && (! all (isequal (cellfun (@numel, S(1).subs, ...
                              "UniformOutput", false){:}))))
                      error ("scatteredInterpolant: query point vectors must have equal length");

                    ## and nd arrays to be equal size
                    elseif (! isequal (cellfun (@size, S(1).subs, "UniformOutput", false){:}))
                      error ("scatteredInterpolant: query point inputs must have equal size");
                    endif

                    ## set output size based on first input element
                    if (! iscolumn (S(1).subs{1}))
                      sz_output = size(S(1).subs{1});
                      reshape_flag = true;
                    endif

                    switch this.dimension
                      case 2
                        qpts = [S(1).subs{1}(:), S(1).subs{2}(:)];
                      case 3
                        qpts = [S(1).subs{1}(:), S(1).subs{2}(:), S(1).subs{3}(:)];
                    endswitch

                  otherwise
                    ## must be 1,2, or 3 input elements. call query usage error
                    print_query_points_usage (this);
                endswitch

                ## tsearchn outputs vector of containing simplex, or NaN for Outside point
                ## plus barycentric coordinates of point within the simplex
                ## vi = [barycentric_coord].[local.Values]

                v = NaN (rows(qpts), 1);

                ##Perform interpolation and extrapolation according to method
                ## note - Matlab provides little indication of the underlying
                ## linear extrapolation method, compatibility unsure.
                switch this.Method
                  case "nearest"
                    ## add an IF to even check extrap
                    switch this.ExtrapolationMethod
                      case "none"
                        nearest_pt_idx = dsearchn (this.Points, this.tri, qpts, NaN);
                        v(! isnan (nearest_pt_idx)) = ...
                          this.Values(nearest_pt_idx(! isnan (nearest_pt_idx)));

                      case "nearest"
                        nearest_pt_idx = dsearchn (this.Points, this.tri, qpts);
                        v = this.Values(nearest_pt_idx);

                      case "linear"
                        outside_qpt_idx = isnan (tsearchn (this.Points, this.tri, qpts));
                        nearest_pt_idx = dsearchn (this.Points, this.tri, qpts);
                        v(! outside_qpt_idx) = ...
                          this.Values(nearest_pt_idx(! outside_qpt_idx));

                        if (any (outside_qpt_idx))
                          error ('not finished')
                        endif

                      otherwise
                        error ("scatteredInterpolant: invalid EXTRAPOLATIONMETHOD %s", this.ExtrapolationMethod);
                    endswitch


                  case "linear"

                    [tri_with_qpts, qpt_bary_coords] = tsearchn (this.Points, this.tri, qpts);
                    inside_qpt_idx = ! isnan (tri_with_qpts);

                    qpt_vertices = this.tri(tri_with_qpts(inside_qpt_idx), :);

                    if (isvector (qpt_vertices))
                      #ensure vertex values in row vector
                      qpt_vertex_values =  this.Values(qpt_vertices)(:)';
                    else
                      qpt_vertex_values =  this.Values(qpt_vertices);
                    endif

                    v(inside_qpt_idx) = dot (qpt_vertex_values, qpt_bary_coords(inside_qpt_idx,:), 2);


                    switch this.ExtrapolationMethod
                      case "none"
                        #nothing to do. NaN already pre-filled for outside points.

                      case "nearest"
                        outside_qpt_idx = ! inside_qpt_idx;
                        nearest_pt_idx = dsearchn (this.Points, this.tri, qpts(outside_qpt_idx,:));
                        v (outside_qpt_idx) = this.Values(nearest_pt_idx);

                      case "linear"
                        outside_qpt_idx = isnan (tsearchn (this.Points, this.tri, qpts));

                        ##check for external pts
                        if (any (outside_qpt_idx))

                          ## find external triangulation points (repeats last pt)
                          boundary_pts = convhull (this.Points);

                          switch this.dimension
                            case 2
                              ## find visible boundary points for each qpt
                              ## return cell array of col vectors for each qpt

                              boundary_pts_vis = find_convex_vis_pts (this, this.Points(boundary_pts(1:end-1),:), qpts(outside_qpt_idx,:));
                              boundary_pts_count = cellfun ('numel', boundary_pts_vis);

                              all_boundary_vis_pts = sort (boundary_pts(unique (cell2mat (boundary_pts_vis'))'));

                              [boundary_pts_vis_idx_to_all, ~] = arrayfun(@(A) find(all_boundary_pts_vis==boundary_pts_vis{A}),[1:numel(boundary_pts_vis)]', "UniformOutput",false);

                              ## compute gradients at visible boundary pts
                              ## find points connecting to each boundary pt
                              [boundary_tris ~] = arrayfun (@(A) find (this.tri==A), all_boundary_vis_pts, "UniformOutput", false);
                              boundary_connections = cellfun (@(A) unique (this.tri(A,:)), boundary_tris, "UniformOutput", false);

                              ## for each boundary point, make gradient system of equations
                              boundary_gradeqns = arrayfun(@(A) [this.Points(boundary_connections{A}(boundary_connections{A}~=A),:) -  this.Points(all_boundary_vis_pts(A),:), ...
                                  this.Values(boundary_connections{A}(boundary_connections{A}~=A)) - this.Values(all_boundary_vis_pts(A))], ...
                                  [1:numel(all_boundary_vis_pts)]', "UniformOutput", false);

                              ## solve overdefined eqn system for [dQdx, dQdy] at each boundary pt (each row n is for boundary pt n)
                              boundary_pt_grads = cell2mat(cellfun (@(A) mldivide(A(:,1:end-1), A(:,end)), boundary_gradeqns', "UniformOutput",false))';

                             ## use visible boundary point gradients to extrapolate value at each point
                             ##do 2-point visibility
                             v(outside_qpt_idx)(boundary_pts_count==2)



                             #do > 2-pt visibility


                              ## split boundary into 'visibility regions' projected from boundary edges.


                              ## perform linear extrapolation based on that region and visible gradients

                            case 3  ##3D

                          endswitch
                        endif
                      otherwise
                        error ("scatteredInterpolant: invalid EXTRAPOLATIONMETHOD %s", this.ExtrapolationMethod);
                    endswitch

                  case "natural"
##                    error ("scatteredInterpolant: 'natural' interpolation method not yet implemented")

                    ##currently requires matgeom package
                     if (! pkg('list','matgeom'){1}.loaded)
                       error ("scatteredInterpolant: 'NATURAL' method requires 'matgeom' package be installed and loaded");
                     endif

                    #identify point/tri associations and index of points inside convex hull
                    [tri_with_qpts, qpt_bary_coords] = tsearchn (this.Points, this.tri, qpts);
                    inside_qpt_idx = !isnan (tri_with_qpts);



                    switch this.dimension
                      case 2
                        ## find circumcircles of delaunay tri (center & rad)
                        ## which are also the voronoi facet intersection points.
##                        [v_facet_pts, v_facet_list] = voronoin (this.Points);
##                        tri_circ_radii = distancePoints (v_facet_pts(2:end,:), this.Points(this.tri(:,1),:), "diag");
##                        ##cirs/spheres containing qp are it's nn points

                        #find centers
                        ## use circumcenters instead of voronoi to keep data
                        ## matched to delaunay tri's.
                        tri_circ_centers = circumCenter (this.Points(this.tri(:,1),:),this.Points(this.tri(:,2),:),this.Points(this.tri(:,3),:));
                        tri_circ_radii = distancePoints (tri_circ_centers, this.Points(this.tri(:,1),:), "diag");

                        ## find which circles cover which qpts. those qpts
                        ## natural neighbors will be all of the points of
                        ## the triangles in those circles.
                        in_which_circle = distancePoints (tri_circ_centers, qpts(inside_qpt_idx,:)) <= tri_circ_radii; #rows are tri_circs, cols are for each qpt

                        ## Get the interpolation points that are nat neighbors
                        ## for each qpt. likely unequal lengths, so store in
                        ## cells.
                        qpt_nat_neighbors_pts =  arrayfun(...
                          @(x) unique(this.tri(in_which_circle(:,x),:))(:),...
                            [1:columns(in_which_circle)],"UniformOutput",false);






                        switch this.ExtrapolationMethod
                          case "none"
                          case "nearest"
                          case "linear"
                          otherwise
                            error ("scatteredInterpolant: invalid EXTRAPOLATIONMETHOD %s", this.ExtrapolationMethod);
                        endswitch

                      case 3
                        switch this.ExtrapolationMethod
                          case "none"
                          case "nearest"
                          case "linear"
                          otherwise
                            error ("scatteredInterpolant: invalid EXTRAPOLATIONMETHOD %s", this.ExtrapolationMethod);
                        endswitch
                    endswitch
                endswitch

                if (reshape_flag)
                  v = reshape (v, sz_output);
                endif

              endif

            case {".", "{"}
              v = builtin ("subsref", this, S(1));

            otherwise
              error ("scatteredInterpolant: invalid scatteredInterpolant index type '%s'", S(1).type);
          endswitch

        else
          v = builtin ("subsref", v, S(S_idx));
        endif
      endfor
    endfunction

    function this = subsasgn (this, S, val)

      # there are only four public properties (even if empty) at object
      # creation. They can be updated via subasgn, but not created, and new
      # properties cannot be added. Points and values may be individually updated
      # to have a count mismatch without error. Error will be produced on attempted
      # interpolation.
      #
      # subsasgn should error on attempt to add new property, add invalid method, or
      # add poorly formed points/values. Points/values can be updated individually.
      # so numrows doesn't have to match on assignment. But the object should
      # track whether or not the points/values are in a matched state.
      #
      # Points - can be empty, must be a 2 x n (2D) or 3 x n (3D) double array
      # Values - can be empty, must be a 1 x m double vector. if not column, store as column.
      # Method - must be char/string: 'nearest', 'linear', or 'natural' (can't be empty)
      # ExtrapolationMethod - must be char/string: 'none', 'linear', 'nearest' (can't be empty)
      #
      # Updating either points or values should recompute the triangulation.
      update_triangulation = false;
      update_point_value_state = false;
      switch S(1).type(1)
        case "."
          if (numel(S) == 1)
            #if only 1 level to S, validate input values/forms before assignment
            switch (S(1).subs)
              case "Points"
                if (!isnumeric (val))
                  error ("scatteredInterpolant: Points must be numeric");
                elseif ((!isempty(val)) && ((!ismatrix (val)) || (!any (size (val, 2) == [2,3]))))
                  error ("scatteredInterpolant: Points input must be a 2 or 3 column array.");
                endif
                update_triangulation = true;
                update_point_value_state = true;

              case "Values"
                if (!isnumeric (val))
                  error ("scatteredInterpolant: Values must be numeric");
                elseif (!isempty (val))
                  if (!isvector (val))
                    error ("scatteredInterpolant: Values must be in vector form");
                  elseif (!iscolumn (val))
                    val = val(:);
                  endif
                endif
                update_point_value_state = true;

              case "Method"
                if (!ischar (val))
                  error ("scatteredInterpolant: METHOD input must be a string");
                elseif (!any (strcmp (val, {"linear", "nearest", "natural"})))
                  error ("scatteredInterpolant: invalid METHOD '%s'", ...
                          val);
                endif

              case "ExtrapolationMethod"
                if (!ischar (val))
                  error ("scatteredInterpolant: EXTRAPOLATIONMETHOD input must be a string");
                elseif (!any (strcmp (val, {"linear","nearest","none"})))
                  error ("scatteredInterpolant: invalid EXTRAPOLATIONMETHOD '%s'", ...
                          val);
                endif

              otherwise
                error ("scatteredInterpolant: invalid property '%s'", ...
                        S(1).subs);
            endswitch
            this = builtin ("subsasgn", this, S, val);

          elseif (numel(S) == 2)
          ## if 2nd index level, must be (), check value before calling builtin
            switch S(2).type(1)
              case "("
                switch (S(1).subs)
                  case "Points"
                    if (!isnumeric (val))
                      error ("scatteredInterpolant: Points must be numeric");
                    elseif (isempty (val) || isempty (this.Points))
                      # empty values can cause invalid array shape
                      testval = builtin ("subsasgn", this.Points, S(2), val);
                      if (!any (size (testval, 2) == [2,3]))
                        error ("scatteredInterpolant: invalid Points shape");
                      endif
                    endif
                    update_triangulation = true;
                    update_point_value_state = true;

                  case "Values"
                    if (!isnumeric (val))
                      error ("scatteredInterpolant: Values must be numeric");
                    elseif ((isempty (val) || isempty (this.Values)))
                      # empty values can cause invalid array shape
                      testval = builtin ("subsasgn", this.Values, S(2), val);
                      if (!isvector (testval))
                        error ("scatteredInterpolant: invalid Values shape");
                      elseif (isempty (this.Values) && !isempty (val) && !iscolumn (val))
                        val = val(:);
                      endif
                    endif
                    update_point_value_state = true;

                  case "Method"
                    if (!ischar (val))
                      error ("scatteredInterpolant: METHOD input must be a string");
                    endif

                    #test substitution, error if it results in a bad value
                    testval = builtin ("subsasgn", this.Method, S(2), val);
                    if (!any (strcmp (testval, {"linear", "nearest", "natural"})))
                      error ("scatteredInterpolant: invalid METHOD '%s'", ...
                              testval);
                    endif

                  case "ExtrapolationMethod"
                    if (!ischar (val))
                      error ("scatteredInterpolant: EXTRAPOLATIONMETHOD input must be a string");
                    endif

                    #test substitution, error if it results in a bad value
                    testval = builtin ("subsasgn", this.ExtrapolationMethod, S(2), val);
                    if (!any (strcmp (val, {"linear","nearest","none"})))
                      error ("scatteredInterpolant: invalid EXTRAPOLATIONMETHOD '%s'", ...
                              testval);
                    endif

                  otherwise
                    error ("scatteredInterpolant: invalid property '%s'", ...
                            S(1).subs);
                endswitch
                this = builtin ("subsasgn", this, S, val);

              case "."
                error ("scatteredInterpolant: scatteredInterpolant properties cannot be subindexed by '.'");
              case "{"
                error ("scatteredInterpolant: scatteredInterpolant cannot be indexed by {}");
              otherwise
            endswitch

            else
              ## no valid function for levels greater than 2
              error ("scatteredInterpolant: scatteredInterpolant assignment depth > 2 undefined");
            endif

            if (update_triangulation);
              ## update triangulation only if Points changed
              this = setTriangulation (this);
            endif
            if (update_point_value_state)
              ##update p/v state for points or values change
              this = check_points (this);
            endif

        case "("
          # () indexing reserved for interpolation. No object array definition,
          #  cannot be used for array assignment.
          error ("scatteredInterpolant: scatteredInterpolant array value assignment undefined");
        case "{"
          error ("scatteredInterpolant: scatteredInterpolant cannot be indexed by {}");
        otherwise
          error ("scatteredInterpolant: invalid scatteredInterpolant index type '%s'", S(1).type);
      endswitch
    endfunction

##    function disp (this)
##      if (nargout > 0)
##        error ("scatteredInterpolant: output assignment not defined");
####      elseif (nargin != 1)
####        error ("scatteredInterpolant: only defined for one input.");
####      elseif (! strcmp (class (this), "scatteredInterpolant"))
####        error ("scatteredInterpolant: only defined for scatteredInterpolant objects.")
##      endif
##      printf('%s object with properties:\n\n', class (this));
##      printf('    Points = %dx%d matrix\n', rows(this.Points), columns(this.Points));
##      printf('    Values = %dx%d matrix\n', rows(this.Values), columns(this.Values));
##      printf('    Method = %s\n', this.Method);
##      printf('    ExtrapolationMethod = %s\n', this.ExtrapolationMethod);
##    endfunction

  endmethods

  methods (Access = public, Hidden = false)

    function this = setTriangulation (this)
      ##called by constructor and subsasgn whenever tri needs a recalc

      if (isempty (this.Points))
        this.tri = [];
      else
        #compute delaunay triangulation. if invalid, catch error but continue
        try
          ##disable qhull warnings??
          this.tri = delaunayn (this.Points);
        catch
          this.tri = [];
        end_try_catch
      endif

      if (isempty (this.tri) || isequal(this.tri, 0))
        this.valid_tri = false;
      else
        this.valid_tri = true;
      endif
    endfunction

    function this = check_points (this)
      this.dimension = columns (this.Points); ## 0 (empty) or 2/3 for 2D/3D
      this.valid_points_vals = isequal (rows (this.Points), rows (this.Values));
      this.enough_points = this.dimension && (rows (this.Points) > this.dimension);
    endfunction

    function visible_points = find_convex_vis_pts (this, pts, viewpts)
      %% pts = pointset given as 2d array, 1 pt per row. [x1, y1, (z1); x2, y2 (z2); etc]
      %% pts already assumed to be a convex hull, in CCW order, without repetition.
      %% viewpt = viewpoint origin (campera location) as row vector [x y (z)]

      #recast multiple viewpts along dim3
      viewpts = permute (viewpts, [3 2 1]);

      % shift all points so that viewpt is at the origin (simplifies vector math)
      pts = pts - viewpts;
      [num_pts, dim, num_vpts] = size(pts);

      % find point center of mass to set reference vector
      pts_com = mean (pts, 1);

      if (dim == 2)
      % 2D calculate sin(theta) for each angle relative to C.O.M. vector
      % as defined, in 2D, neg angles are CCW w.r.t. com vector, pos angle are CW
      % ( sin(theta) = A x B / |A||B| )
        pt_sin_thetas = (pts(:,1,:).*pts_com(:,2,:) - pts(:,2,:).*pts_com(:,1,:)) ...
                          ./ (sqrt (sumsq (pts, 2) .* sumsq (pts_com, 2)));

       [minangle, min_theta_idx] = min (pt_sin_thetas, [], 1);
       [maxangle, max_theta_idx] = max (pt_sin_thetas, [], 1);

       ## check for numerically colinear points. for a convex hull, max/min
       ## angle points can only have same angle as another point if
       ## viewpoint-point vector is colinear with pt-pt edge. in that case, need
       ## to select closer point, with shorter viewpoint-point distance.

        [same_min_idx, same_min_vpt] = find (pt_sin_thetas == minangle);
        [same_max_idx, same_max_vpt] = find (pt_sin_thetas == maxangle);
        for idx = 1 : num_vpts
          repeated = same_min_idx(same_min_vpt == idx);
          if (numel (repeated) > 1)
            [~, closer_min] = min (sumsq (pts(repeated,:,idx), 2));
            min_theta_idx(:,:,idx) = repeated(closer_min);
            minangle(:,:,idx) = pt_sin_thetas(min_theta_idx(:,:,idx));
          endif

          repeated = same_max_idx(same_max_vpt == idx);
          if (numel (repeated) > 1)
            [~, closer_max] = max (sumsq (pts(repeated,:,idx), 2));
            max_theta_idx(:,:,idx) = repeated(closer_max);
            maxangle(:,:,idx) = pt_sin_thetas(max_theta_idx(:,:,idx));
          endif
        endfor

        ## visible points are all those between min&max index. points from
        ## convex hull are in CCW order (cyclical), min-max angle always CW
        ## around viewpoint, which for external viewpoint is always CCW around
        ## convex hull. so points go from [min_theta_idx : max_theta_idx],
        ## including around the limits.

        if (any (min_theta_idx == max_theta_idx))
          error ("only 1 visible point, shouldn't be possible for a 2D convex hull");
        endif

        simple_order = max_theta_idx > min_theta_idx;

        visible_points = cell (num_vpts, 1);
        visible_points (simple_order) = arrayfun (@colon, min_theta_idx(simple_order), max_theta_idx(simple_order), "UniformOutput", false);
        visible_points (!simple_order) = arrayfun(@(A,B) sort([A:num_pts, 1:B]), min_theta_idx(!simple_order), max_theta_idx(!simple_order),"UniformOutput",false);

      else
        # 3D
        error ("no 3D yet");
      endif
    endfunction

##    function print_usage (this)
##
##      ## FIXME: overloading print_usage until bug ___ is fixed for classdefs.
##      ## Can be removed and errors reverted back to 'print_usage ()' once that
##      ## is fixed.
##      ident = "Octave:invalid-fun-call";
##      msg = sprintf(["Invalid call to scatteredInterpolant.  Correct usage is:\n\n", ...
##        "    -- F = scatteredInterpolant\n"...
##        "    -- F = scatteredInterpolant (X, Y, Q)\n"...
##        "    -- F = scatteredInterpolant (X, Y, Z, Q)\n"...
##        "    -- F = scatteredInterpolant (P, Q)\n"...
##        "    -- F = scatteredInterpolant (..., METHOD)\n"...
##        "    -- F = scatteredInterpolant (..., METHOD, EXTRAPOLATIONMETHOD)\n\n"...
##        "Additional help for built-in functions and operators is\n"...
##        "available in the online version of the manual.  Use the command\n"...
##        "'doc <topic>' to search the manual index.\n\n"...
##        "Help and information about Octave is also available on the WWW\n"...
##        "at https://www.octave.org and via the help@octave.org\n"...
##        "mailing list.\n"]);
##
##      error (struct ("message", msg, "identifier", ident, "stack", dbstack (1)));
##    endfunction

    function print_query_points_usage (this)
      msg = sprintf(["scatteredInterpolant: invalid query points form. Correct usage is:\n\n", ...
    "    -- Si = f(xi, yi)\n", ...
    "    -- Si = f(xi, yi, zi)\n", ...
    "    -- Si = f(Pi)\n", ...
    "    -- Si = f({Xg, Yg})\n", ...
    "    -- Si = f({Xg, Yg, Zg})\n\n", ...
    "    See 'help scatteredInterpolant' for more information."]);
      error (struct ("message", msg, "identifier", "", "stack", dbstack (1)));
    endfunction

  endmethods

endclassdef

##TEST CONSTRUCTOR
%!test # empty object with defaults
%! A = scatteredInterpolant ();
%! assert (class (A), "scatteredInterpolant");
%! assert ({A.Points, A.Values, A.Method, A.ExtrapolationMethod}, {[], [], 'linear', 'linear'});

%!test # simple object, class, verify value orientation independence
%! pts = [magic(3); 2*magic(3)];
%! A = scatteredInterpolant (pts, [1:6]');
%! B = scatteredInterpolant (pts, [1:6]);
%! assert ({class(A), class(B)}, {"scatteredInterpolant", "scatteredInterpolant"});
%! assert ({A.Points, A.Values, A.Method, A.ExtrapolationMethod}, {pts, [1:6]', 'linear', 'linear'});
%! assert ({B.Points, B.Values, B.Method, B.ExtrapolationMethod}, {pts, [1:6]', 'linear', 'linear'});

%!test # check that objects unable to produce triangulations still produce object
%! pts = magic (3);
%! A = scatteredInterpolant (pts, [1:3]');
%! warning ("off");
%! assert (class(A), "scatteredInterpolant");
%! assert ({A.Points, A.Values, A.Method, A.ExtrapolationMethod}, {pts, [1:3]', 'linear', 'linear'});
%! warning ("on");

%!test # 2D vector orientation checks
%! A = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4]', [1:5]');
%! B = scatteredInterpolant ([1 2 3 4 5], [5 1 3 2 4]', [1:5]');
%! C = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4], [1:5]');
%! D = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4]', [1:5]);
%! E = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4], [1:5]);
%! F = scatteredInterpolant ([1 2 3 4 5], [5 1 3 2 4]', [1:5]);
%! G = scatteredInterpolant ([1 5; 2 1; 3 3; 4 2; 5 4], [1:5]');
%! H = scatteredInterpolant ([1 5; 2 1; 3 3; 4 2; 5 4], [1:5]);
%! assert (class (A), "scatteredInterpolant");
%! assert ({A.Points, A.Values, A.Method, A.ExtrapolationMethod}, {[1 2 3 4 5; 5 1 3  2 4]', [1:5]', 'linear', 'linear'});
%! assert (isequal (A.Points, B.Points, C.Points, D.Points, E.Points, F.Points, G.Points, H.Points));
%! assert (isequal (A.Values, B.Values, C.Values, D.Values, E.Values, F.Values, G.Values, H.Values));

%!test # 3D vector orientation checks
%! A = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4]', [2 1 5 4 3]', [1:5]');
%! B = scatteredInterpolant ([1 2 3 4 5], [5 1 3 2 4]', [2 1 5 4 3]', [1:5]');
%! C = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4], [2 1 5 4 3]', [1:5]');
%! D = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4]', [2 1 5 4 3], [1:5]');
%! E = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4]', [2 1 5 4 3]', [1:5]);
%! F = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4], [2 1 5 4 3], [1:5]);
%! G = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4]', [2 1 5 4 3], [1:5]);
%! H = scatteredInterpolant ([1 5 2; 2 1 1; 3 3 5; 4 2 4; 5 4 3], [1:5]');
%! I = scatteredInterpolant ([1 5 2; 2 1 1; 3 3 5; 4 2 4; 5 4 3], [1:5]);
%! assert (class (A), "scatteredInterpolant");
%! assert ({A.Points, A.Values, A.Method, A.ExtrapolationMethod}, {[1 5 2; 2 1 1; 3 3 5; 4 2 4; 5 4 3], [1:5]', 'linear', 'linear'});
%! assert (isequal (A.Points, B.Points, C.Points, D.Points, E.Points, F.Points, G.Points, H.Points, I.Points));
%! assert (isequal (A.Values, B.Values, C.Values, D.Values, E.Values, F.Values, G.Values, H.Values, I.Values));

%!test # method input checks
%! A = scatteredInterpolant ([1 5; 2 1; 3 3; 4 2; 5 4], [1:5]', "nearest");
%! B = scatteredInterpolant ([1 5; 2 1; 3 3; 4 2; 5 4], [1:5]', "nearest", "none");
%! C = scatteredInterpolant ([1 5; 2 1; 3 3; 4 2; 5 4], [1:5]', "linear", "none");
%! D = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4]', [2 1 5 4 3], [1:5]', "nearest");
%! E = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4]', [2 1 5 4 3], [1:5]', "nearest", "none");
%! F = scatteredInterpolant ([1 2 3 4 5]', [5 1 3 2 4]', [2 1 5 4 3], [1:5]', "linear", "none");
%! assert ({A.Method, A.ExtrapolationMethod}, {"nearest", "nearest"});
%! assert ({B.Method, B.ExtrapolationMethod}, {"nearest", "none"});
%! assert ({C.Method, C.ExtrapolationMethod}, {"linear", "none"});
%! assert ({D.Method, D.ExtrapolationMethod}, {"nearest", "nearest"});
%! assert ({E.Method, E.ExtrapolationMethod}, {"nearest", "none"});
%! assert ({F.Method, F.ExtrapolationMethod}, {"linear", "none"});

## Test input validation
%!error scatteredInterpolant (1)
%!error scatteredInterpolant (1, 2, 3, 4, 5, 6, 7)
%!error scatteredInterpolant (1, 2, 3, 4, 5, 6)
%!error scatteredInterpolant (1, 2, "abc", "def", "ghi")
%!error <input points and values cannot be empty> scatteredInterpolant ([], [])
%!error <input points and values cannot be empty> scatteredInterpolant ([], [1 2 3]')
%!error <input points and values cannot be empty> scatteredInterpolant (magic(3), [])
%!error <inputs must be numeric arrays or> scatteredInterpolant (1, 1, false, 1, "abc" , "def")
%!error scatteredInterpolant (1, 1, 1, "abc", 1, "def")
%!error <Value input must be a numeric vector> scatteredInterpolant (magic (3), magic (3))
%!error <must have the same number of rows> scatteredInterpolant (magic (3), [1:4])
%!error <must have the same number of rows> scatteredInterpolant (magic (3), [1:2])
%!error <must have the same number of rows> scatteredInterpolant ([1:8]', [1:10]', [1:10]', [1:10]')
%!error <must have the same number of rows> scatteredInterpolant ([1:10]', [1:8]', [1:10]', [1:10]')
%!error <must have the same number of rows> scatteredInterpolant ([1:10]', [1:10]', [1:8]', [1:10]')
%!error <must have the same number of rows> scatteredInterpolant ([1:10]', [1:10]', [1:10]', [1:8]')
%!error <Point input must be a 2 or 3> scatteredInterpolant ([1:3]', [1:3]')
%!error <Point input must be a 2 or 3> scatteredInterpolant ([[1:4];[3:6];[4:7]], [1:3]')
%!error <Point input must be a 2 or 3> scatteredInterpolant (cat(3,magic(3),magic(3)), [1:3]')
%!error <Invalid METHOD 'foo'> scatteredInterpolant (magic (3), [1:3], "foo")
%!error <Invalid METHOD 'foo'> scatteredInterpolant (magic (3), [1:3], "foo", "none")
%!error <Invalid METHOD 'foo'> scatteredInterpolant (magic (3), [1:3], "foo", "bar")
%!error <Invalid EXTRAPOLATIONMETHOD 'foo'> scatteredInterpolant (magic (3), [1:3], "linear", "foo")


##TEST SUBSREF AND INTERPOLATION FUNCTION
%!assert (subsref (scatteredInterpolant (), substruct (".", "Points")), [])
%!assert (subsref (scatteredInterpolant (), substruct (".", "Values")), [])
%!assert (subsref (scatteredInterpolant (), substruct (".", "Method")), "linear")
%!assert (subsref (scatteredInterpolant (), substruct (".", "ExtrapolationMethod")), "linear")
%!assert (subsref (scatteredInterpolant (), substruct (".", "Method", "()", {2})), "i")

%!error subsref (scatteredInterpolant (), substruct (".", "blah"))
%!error subsref (scatteredInterpolant (), substruct ("{}", {1, 1}))
%!error <invalid scatteredInterpolant index type> subsref (scatteredInterpolant (), struct ("type", "/", "subs", 1))

%!test ## check point & triangulation errors
%! A = scatteredInterpolant ([magic(3); 2*magic(3)],[1:6]');
%! fail ("A(1,2)", "query points dimension");
%! A.Values = [1:5]';
%! fail ("A(1,2,3)", "unequal number");

%!test ##check triangulation warnings
%! A = scatteredInterpolant ([[1:5]',[1:5]'], [1:5]');
%! fail ("A(1,2)", "warning", "cannot calculate triangulation");

%!test ##check point sufficiency warnings
%! A = scatteredInterpolant (magic(3), [1:3]');
%! fail ("A(1,2,3)", "warning", "not enough points");
%! A = scatteredInterpolant (magic(2), [1:2]');
%! fail ("A(1,2)", "warning", "not enough points");

##interpolation input handling
%!test
%! A = scatteredInterpolant ([magic(3); 2*magic(3)], [1:6]);
%! qpts = cat (3, [1 2 3], [4 5 6]);
%! fail ("A(qpts)", "query points must be 2D vectors or arrays");
%! fail ("A([1,2])", "query points dimension must match");
%! fail ("A({'abc'})", "must be numeric");
%! fail ("A({[1,2]})", "query grid vector count must match");
%! fail ("A({[1,2,3]})", "query grid vector count must match");
%! fail ("A({[1,2],[3,4]})", "query grid vector count must match");
%! fail ("A({[1,2;3,4],1,1})", "query grid vectors must be row or column");
%! fail ("A(1,2)", "query points dimension must match interpolant");
%! fail ("A(1,2,[3 4])", "query point vectors must have equal length");
%! fail ("A([1 2 3],[1 2],[3 4 5])", "query point vectors must have equal length");
%! fail ("A([1 2 3]',[1 2]',[3 4 5]')", "query point vectors must have equal length");
%! fail ("A([1 2 3],[1 2],[3 4 5]')", "query point vectors must have equal length");
%! fail ("A([1,2;3,4],[1,2;3,4],[1,2])", "query point inputs must have equal size");
%! fail ("A(qpts,qpts,[1 2 3])", "query point inputs must have equal size");
%! fail ("A('foo')", "invalid query points form");
%! fail ("A({1 2 3},{1 2 3})", "invalid query points form");
%! fail ("A(1,2,3,4)", "invalid query points form");

##interpolation output tests
%!test
%! pts = [0 0 0; 0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 1 1 1; 0.5 0.5 0.5];
%! A = scatteredInterpolant (pts, [1:9]', 'nearest', 'none');
%! qpts = magic(3)/10 -0.2;
%! assert (A(qpts), [NaN, 9, 3]')
%! A.ExtrapolationMethod = "nearest";
%! assert (A(qpts), [4, 9, 3]')


####TEST DISP
##%!error <only defined for one input> disp (scatteredInterpolant (), 1)
##%!error <output assignment not defined> B = disp (scatteredInterpolant)



##TEST SUBSASGN
#Test first level assignment
%!test
%! A = scatteredInterpolant ();
%! A.Points = [magic(3); 2*magic(3)];
%! A.Values = [1:6]';
%! A.Method = "linear";
%! A.ExtrapolationMethod = "none";
%! assert ({A.Points, A.Values, A.Method, A.ExtrapolationMethod}, {[magic(3); 2*magic(3)], [1:6]', "linear", "none"});
%! A.Values = [1:6];
%! assert ({A.Points, A.Values, A.Method, A.ExtrapolationMethod}, {[magic(3); 2*magic(3)], [1:6]', "linear", "none"});

#Test second level assignments
%!test
%! A = scatteredInterpolant(magic (2), [1:2]', "nearest", "none");
%! A.Points(3:4,:) = [3 4; 5 6];
%! assert (A.Points, [4 3; 1 2; 3 4; 5 6]);
%! A.Points(:,3) = [0; 7; 2; 6];
%! assert (A.Points, [4 3 0; 1 2 7; 3 4 2; 5 6 6]);
%! A.Points(1:2,1:2) = [0 0; 0 0];
%! assert (A.Points, [0 0 0; 0 0 7; 3 4 2; 5 6 6]);
%! A.Points(6) = 10;
%! assert (A.Points, [0 0 0; 0 10 7; 3 4 2; 5 6 6]);
%! A.Points(6:8) = [9 8 7];
%! assert (A.Points, [0 0 0; 0 9 7; 3 8 2; 5 7 6]);
%! A.Values(2) = 4;
%! assert (A.Values, [1; 4]);
%! A.Values(3:4) = [3, 4];
%! assert (A.Values, [1; 4; 3; 4]);
%! A.Values(8) = 8;
%! assert (A.Values, [1; 4; 3; 4; 0; 0; 0; 8]);
%! A.Method(2:end) = "atural";
%! assert (A.Method, "natural");


#Test valid empty assignments for Points and Values
%!test
%! warning("off");
%! A = scatteredInterpolant ([magic(3); 2*magic(3)], [1:6]);
%! B = scatteredInterpolant ([magic(3); 2*magic(3)], [1:6]);
%! A.Points = [];
%! assert (A.Points, []);
%! A.Values = [];
%! assert (A.Values, []);
%! B.Points(3,:) = [];
%! assert (B.Points, [8 1 6; 3 5 7; 2*magic(3)]);
%! B.Points (:,2) = [];
%! assert (B.Points, [8 6; 3 7; 16 12; 6 14; 8 4]);
%! B.Values(2) = [];
%! assert (B.Values, [1 3 4 5 6]');
%! warning("on");

## Test input validation
%!error subsasgn (scatteredInterpolant ([magic(3); 2*magic(3)], [1:6]), substruct (".", "Points", "()", {4}), [1 2])
%!error <invalid Points shape> subsasgn (scatteredInterpolant ([magic(3); 2*magic(3)], [1:6]), substruct (".", "Points", "()", {4}), [])
%!error <invalid Points shape> subsasgn (scatteredInterpolant (), substruct (".", "Points", "()", {1:3,1:4}), rand (3,4))
%!error <Points must be numeric> subsasgn (scatteredInterpolant (), substruct (".", "Points"), "abc")
%!error <Points must be numeric> subsasgn (scatteredInterpolant (), substruct (".", "Points"), true)
%!error <Points must be numeric> subsasgn (scatteredInterpolant (), substruct (".", "Points", "()", {1}), "a")
%!error <Points must be numeric> subsasgn (scatteredInterpolant (), substruct (".", "Points", "()", {1}), true)
%!error <Points must be numeric> subsasgn (scatteredInterpolant (), substruct (".", "Points"), "")
%!error <Points input must be a 2 or 3> subsasgn (scatteredInterpolant (), substruct (".", "Points"), [1:4]')
%!error <Points input must be a 2 or 3> subsasgn (scatteredInterpolant (), substruct (".", "Points"), magic (4))
%!error <Points input must be a 2 or 3> subsasgn (scatteredInterpolant (), substruct (".", "Points"), cat (3, magic (3), magic(3)))
%!error <invalid Values shape> subsasgn (scatteredInterpolant (), substruct (".", "Values", "()", {1:2,1:2}), magic (2))
%!error <Values must be numeric> subsasgn (scatteredInterpolant (), substruct (".", "Values"), "abc")
%!error <Values must be numeric> subsasgn (scatteredInterpolant (), substruct (".", "Values"), true)
%!error <Values must be numeric> subsasgn (scatteredInterpolant (), substruct (".", "Values", "()", {1}), "a")
%!error <Values must be numeric> subsasgn (scatteredInterpolant (), substruct (".", "Values", "()", {1}), true)
%!error <Values must be numeric> subsasgn (scatteredInterpolant (), substruct (".", "Values"), "")
%!error <Values must be in vector form> subsasgn (scatteredInterpolant (), substruct (".", "Values"), [1 1; 2 2])
%!error <Values must be in vector form> subsasgn (scatteredInterpolant (), substruct (".", "Values"), cat(3,[1 2 3],[1 2 3]))
%!error <METHOD input must be a string> subsasgn (scatteredInterpolant (), substruct (".", "Method"), 1)
%!error <METHOD input must be a string> subsasgn (scatteredInterpolant (), substruct (".", "Method"), [])
%!error <METHOD input must be a string> subsasgn (scatteredInterpolant (), substruct (".", "Method", "()", {1}), 1)
%!error <invalid METHOD> subsasgn (scatteredInterpolant (), substruct (".", "Method"), "foo")
%!error <invalid METHOD> subsasgn (scatteredInterpolant (), substruct (".", "Method"), "")
%!error <invalid METHOD> subsasgn (scatteredInterpolant (), substruct (".", "Method", "()", {1}), "f")
%!error <EXTRAPOLATIONMETHOD input must be a string> subsasgn (scatteredInterpolant (), substruct (".", "ExtrapolationMethod"), 1)
%!error <EXTRAPOLATIONMETHOD input must be a string> subsasgn (scatteredInterpolant (), substruct (".", "ExtrapolationMethod"), [])
%!error <EXTRAPOLATIONMETHOD input must be a string> subsasgn (scatteredInterpolant (), substruct (".", "ExtrapolationMethod", "()", {1}), 1)
%!error <invalid EXTRAPOLATIONMETHOD> subsasgn (scatteredInterpolant (), substruct (".", "ExtrapolationMethod"), "foo")
%!error <invalid EXTRAPOLATIONMETHOD> subsasgn (scatteredInterpolant (), substruct (".", "ExtrapolationMethod"), "")
%!error <invalid EXTRAPOLATIONMETHOD> subsasgn (scatteredInterpolant (), substruct (".", "ExtrapolationMethod", "()", {1}), "f")
%!error <invalid property> subsasgn (scatteredInterpolant (), substruct (".", "foo"), 1)
%!error <invalid property> subsasgn (scatteredInterpolant (), substruct (".", "foo", "()", {1}), 1)
%!error <scatteredInterpolant properties cannot be sub> subsasgn (scatteredInterpolant (), substruct (".", "Method" , ".", "foo"), 1)
%!error <scatteredInterpolant cannot be indexed by> subsasgn (scatteredInterpolant (), substruct (".", "Method" , "{}", {1, 1}), 1)
%!error <scatteredInterpolant assignment depth> subsasgn(scatteredInterpolant (), substruct (".", "Method", "()", {1}, "()", {1}))
%!error <scatteredInterpolant array value assignment undefined> subsasgn (scatteredInterpolant (), substruct ("()", {2}), 1)
%!error <scatteredInterpolant cannot be indexed by> subsasgn (scatteredInterpolant (), substruct ("{}", {1, 1}), 1)
%!error <invalid scatteredInterpolant index type> subsasgn (scatteredInterpolant (), struct ("type", "/", "subs", 1), 1)
