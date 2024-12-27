    function visible_points = find_convex_vis_pts (pts, viewpts)
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
        pt_sin_thetas = (pts(:,1,:).* pts_com(:,2,:) - pts(:,2,:) .* pts_com(:,1,:))./...
                          (sqrt (sumsq (pts, 2) .* sumsq (pts_com, 2)));

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
