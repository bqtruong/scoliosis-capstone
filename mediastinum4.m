% Use with pulmonarytoolkit
function [ms] = mediastinum4(lungs)
	% Load lung segmentation from ptk
	% S = load('pedslungs.mat');
 %    lungs = S.image;
    ms = zeros(size(lungs));

    llung = lungs; rlung = lungs;
    llung(llung == 2) = 0;
    rlung(rlung == 1) = 0; rlung(rlung == 2) = 1;

 %    % Artificially shifting rlung down to test special case
 %    replace_rlung = zeros(size(rlung));
	% for z = size(rlung,3)-50:-1:1
	% 	for y = 1:size(rlung,2)
	% 		for x = 1:size(rlung,1)
	% 			if rlung(x,y,z) == 1
	% 				replace_rlung(x,y,z+50) = 1;
	% 			end
	% 		end
	% 	end
	% end
	% rlung = replace_rlung;

    first_detected = 'none';
    llung_start_idx = NaN;
    rlung_start_idx = NaN;

    % Slice to get starting indices for the left and right lung
    for z = 1:size(lungs,3)
    	% Slice left lung
    	llung_slice = llung(:,:,z);

    	% Slice right lung
    	rlung_slice = rlung(:,:,z);

    	if strcmp(first_detected, 'none')
    		if ~isempty(find(llung_slice == 1, 1)) && ~isempty(find(rlung_slice == 1, 1))
    			first_detected = 'both';
    			llung_start_idx = z;
    			rlung_start_idx = z;
    			break;
			elseif ~isempty(find(llung_slice == 1, 1))
				first_detected = 'left';
				llung_start_idx = z;
			elseif ~isempty(find(rlung_slice == 1, 1))
				first_detected = 'right';
				rlung_start_idx = z;
			end
		end
		if strcmp(first_detected, 'left') && isnan(rlung_start_idx) && ~isempty(find(rlung_slice == 1, 1))
			rlung_start_idx = z;
			break;
		end
		if strcmp(first_detected, 'right') && isnan(llung_start_idx) && ~isempty(find(llung_slice == 1, 1))
			llung_start_idx = z;
			break;
		end
    end

    % Get ending indices for the left and right lung (diaphragm)
    % Get XZ slice with most lung volume (appx center of XY)
    lungs2 = lungs;
    for z = 1:size(lungs2,3)
    	lungs2(:,:,z) = imfill(lungs2(:,:,z));
    	llung_slice = lungs2(:,:,z);
    	llung_slice(llung_slice == 2) = 0;
    	llung_slice = imbinarize(llung_slice);
    	llung_slice = bwareafilt(llung_slice, 1);
    	rlung_slice = lungs2(:,:,z);
    	rlung_slice(rlung_slice == 1) = 0;
    	rlung_slice(rlung_slice == 2) = 1;
    	rlung_slice = imbinarize(rlung_slice);
    	rlung_slice = bwareafilt(rlung_slice, 1);
        rlung_slice = double(rlung_slice);
    	rlung_slice(rlung_slice == 1) = 2;
    	lungs2(:,:,z) = double(llung_slice) + rlung_slice;
    end
    % Get XZ slice with greatest volume, assuming fullest view of diaphragm
	XZ_max = -Inf;
	XZ_max_val = -Inf;
	for x = 1:size(lungs2,1)
		if isequal(unique(lungs2(x,:,:)), [0; 1; 2])
			lungs2_slice = lungs2(x,:,:);
			lungs2_slice(lungs2_slice == 2) = 1;
			if sum(lungs2_slice, 'all') > XZ_max_val
				XZ_max = x;
                XZ_max_val = sum(lungs2_slice, 'all');
			end
		end
    end
    llung = lungs2;
    rlung = lungs2;
    llung(llung == 2) = 0;
    rlung(rlung == 1) = 0;
    rlung(rlung == 2) = 1;
    %rlung = lungs2(lungs2 == 2);
    %rlung(rlung == 2) = 1;
    %lungs2(lungs2 == 2) = 1;
    XZ_left = squeeze(llung(XZ_max,:,:));
    XZ_right = squeeze(rlung(XZ_max,:,:));
    XZ_left = imfill(XZ_left);
    XZ_right = imfill(XZ_right);
    llung_end_idx = Inf;
    rlung_end_idx = Inf;
    
    left_start = false;
    right_start = false;
    for z = size(lungs, 3):-1:1
        left_conn = bwconncomp(XZ_left(:,z));
        right_conn = bwconncomp(XZ_right(:,z));
        if left_conn.NumObjects == 2 && left_start == false
            left_start = true;
        end
        if right_conn.NumObjects == 2 && right_start == false
            right_start = true;
        end
        if left_conn.NumObjects == 1 && left_start == true && isinf(llung_end_idx)
            llung_end_idx = z;
        end
        if right_conn.NumObjects == 1 && right_start == true && isinf(rlung_end_idx)
            rlung_end_idx = z;
        end
    end

	%llung_end_idx = min(left_bound(:,2));
	%rlung_end_idx = min(right_bound(:,2));

    % Get central portion of mediastinum
    ms = central_body(llung, rlung, [llung_start_idx, llung_end_idx], [rlung_start_idx, rlung_end_idx]);

    % Get top portion of mediastinum
    if strcmp(first_detected, 'left')
    	ms = ms + top_body_left(llung, rlung, [llung_start_idx, llung_end_idx], [rlung_start_idx, rlung_end_idx]);
	elseif strcmp(first_detected, 'right')
		ms = ms + top_body_right(llung, rlung, [llung_start_idx, llung_end_idx], [rlung_start_idx, rlung_end_idx]);
	end
	
	% For viewing, flip upside down (original has top of lungs as 1 index)
	ms = flip(ms,3);
end

% Get mediastinum where left and right lungs are both present in slices
function [ms_w_middle] = central_body (llung, rlung, llung_idx, rlung_idx)
	ms_w_middle = zeros(size(llung));

	if llung_idx(1) < rlung_idx(1)
		start_idx = rlung_idx(1);
	else
		start_idx = llung_idx(1);
	end
	if llung_idx(2) < rlung_idx(2)
		end_idx = llung_idx(2) - 1;
	else
		end_idx = rlung_idx(2) - 1;
	end

	% Get areas for slices containing both left and right lungs
	for z = start_idx:end_idx
    	% Final slice and template slice
    	ms_slice = zeros(size(llung(:,:,z)));

    	% Slice left lung
    	llung_slice = llung(:,:,z);

    	% Slice right lung
    	rlung_slice = rlung(:,:,z);

		% Fill in holes
		llung_slice = imfill(llung_slice);
		llung_slice = imbinarize(llung_slice);
		llung_slice = bwareafilt(llung_slice, 1);

    	% Get top-right most left lung
    	[llungr, ~] = find(llung_slice == 1);
    	[ltop_rightr, ~] = min(llungr);
        [~, ltop_rightc] = find(llung_slice(ltop_rightr,:) == 1, 1, 'last');

    	% Get bottom-right most left lung
    	[lbottom_rightr, ~] = max(llungr);
    	[~, lbottom_rightc] = find(llung_slice(lbottom_rightr,:) == 1, 1, 'last');

    	% Trace boundary
		lx = ltop_rightr; ly = ltop_rightc;
		Bl = bwtraceboundary(llung_slice, [lx ly], 'SE', 8, Inf, 'clockwise');	   
		interior_leftx = find(Bl(:,1) == lbottom_rightr);

		if isempty(interior_leftx)
			break;
		end

		interior_lefty = find(Bl(interior_leftx,2) == lbottom_rightc);
		interior_lefty = interior_lefty + interior_leftx(interior_lefty);
		interior_boundaryl = Bl(1:interior_lefty,:);

		% Fill in holes
		rlung_slice = imfill(rlung_slice);
		rlung_slice = imbinarize(rlung_slice);
		rlung_slice = bwareafilt(rlung_slice, 1);

    	% Get top-left most right lung
    	[rlungr, ~] = find(rlung_slice == 1);
    	[rtop_leftr, ~] = min(rlungr);
        [~, rtop_leftc] = find(rlung_slice(rtop_leftr,:) == 1, 1);

    	% Get bottom-left most right lung
    	[rbottom_leftr, ~] = max(rlungr);
    	[~, rbottom_leftc] = find(rlung_slice(rbottom_leftr,:) == 1, 1);

		% Traverse interior of right lung wall
		rx = rbottom_leftr; ry = rbottom_leftc;
		Br = bwtraceboundary(rlung_slice, [rx ry], 'NW', 8, Inf, 'clockwise');	   
		interior_rightx = find(Br(:,1) == rtop_leftr);

		if isempty(interior_rightx)
			break;
		end

		interior_righty = find(Br(interior_rightx,2) == rtop_leftc);
		interior_righty = interior_righty + interior_rightx(interior_righty);
		interior_boundaryr = Br(1:interior_righty,:);

		% Get four corner points of mediastinum
		Px = [interior_boundaryl(:,1); interior_boundaryr(:,1); interior_boundaryl(1,1)];
		Py = [interior_boundaryl(:,2); interior_boundaryr(:,2); interior_boundaryl(1,2)];

		% Slice template for points needed to be filled (or not)
    	unknown = zeros(size(llung,1), size(llung,2));

    	% Only need to search in bounds of borders
    	unknownxrange = min(Px):max(Px);
    	unknownyrange = min(Py):max(Py);

    	% Mark rectangular range in matrix to search
    	unknown(unknownxrange, unknownyrange) = 1;

    	% Get coords of to-be-searched points
    	unknown_indices = find(unknown == 1);
    	[uix, uiy] = ind2sub(size(unknown), unknown_indices);

    	% Fill borders
    	in = insidepoly(uix, uiy, Px, Py);

    	% Store filled MS
    	known = zeros(size(llung,1), size(llung,2));
    	known(unknown_indices) = in;
    	ms_w_middle(:,:,z) = known;
        % figure(1); imagesc(known);
	end
end

% Left lung tilted above right lung
function [ms_w_top] = top_body_left (llung, rlung, llung_idx, rlung_idx)
	ms_w_top = zeros(size(llung));

	start_idx = llung_idx(1);
	end_idx = rlung_idx(1);

	rlung_slice = rlung(:,:,end_idx);
	rlung_slice = imfill(rlung_slice);
	rlung_slice = imbinarize(rlung_slice);
	rlung_slice = bwareafilt(rlung_slice, 1);

	% Get top-left most right lung
	[rlungr, ~] = find(rlung_slice == 1);
	[rtop_leftr, ~] = min(rlungr);
    [~, rtop_leftc] = find(rlung_slice(rtop_leftr,:) == 1, 1);

	% Get bottom-left most right lung
	[rbottom_leftr, ~] = max(rlungr);
	[~, rbottom_leftc] = find(rlung_slice(rbottom_leftr,:) == 1, 1);

	% Traverse interior of right lung wall
	rx = rbottom_leftr; ry = rbottom_leftc;
	Br = bwtraceboundary(rlung_slice, [rx ry], 'NW', 8, Inf, 'clockwise');	   
	interior_rightx = find(Br(:,1) == rtop_leftr);
	interior_righty = find(Br(interior_rightx,2) == rtop_leftc);
	interior_righty = interior_righty + interior_rightx(interior_righty);
	interior_boundaryr = Br(1:interior_righty,:);

	for z = start_idx:end_idx
		llung_slice = llung(:,:,z);

		% Fill in holes
		llung_slice = imfill(llung_slice);
		llung_slice = imbinarize(llung_slice);
		llung_slice = bwareafilt(llung_slice, 1);

    	% Get top-right most left lung
    	[llungr, ~] = find(llung_slice == 1);
    	[ltop_rightr, ~] = min(llungr);
        [~, ltop_rightc] = find(llung_slice(ltop_rightr,:) == 1, 1, 'last');

    	% Get bottom-right most left lung
    	[lbottom_rightr, ~] = max(llungr);
    	[~, lbottom_rightc] = find(llung_slice(lbottom_rightr,:) == 1, 1, 'last');

    	% Trace boundary
		lx = ltop_rightr; ly = ltop_rightc;
		Bl = bwtraceboundary(llung_slice, [lx ly], 'SE', 8, Inf, 'clockwise');	   
		interior_leftx = find(Bl(:,1) == lbottom_rightr);
		interior_lefty = find(Bl(interior_leftx,2) == lbottom_rightc);
		interior_lefty = interior_lefty + interior_leftx(interior_lefty);
		interior_boundaryl = Bl(1:interior_lefty,:);

		for idx1 = 1:size(interior_boundaryl,1)
			for idx2 = 1:size(interior_boundaryr,1)
				[X, Y, Z] = bresenham_line3d([interior_boundaryl(idx1, 1) interior_boundaryl(idx1, 2) z], [interior_boundaryr(idx2, 1) interior_boundaryr(idx2, 2) end_idx]);
				linearIdxs = sub2ind(size(ms_w_top), X, Y, Z);
				ms_w_top(linearIdxs) = 1;
			end
		end
	end
	for z = start_idx:end_idx
		ms_top_slice = ms_w_top(:,:,z);
        [x, y] = find(ms_top_slice == 1);
        k = boundary(x, y, 1);
        unknown = zeros(size(ms_top_slice));
        unknown(min(x(k)):max(x(k)), min(y(k)):max(y(k))) = 1;
        unknown_indices = find(unknown == 1);
    	[uix, uiy] = ind2sub(size(unknown), unknown_indices);
        in = insidepoly(uix, uiy, x(k), y(k));
        known = zeros(size(ms_top_slice));
        known(unknown_indices) = in;
        known = imfill(known);
        ms_w_top(:,:,z) = known;
	end
end

% Right lung tilted above left lung
function [ms_w_top] = top_body_right (llung, rlung, llung_idx, rlung_idx)
	ms_w_top = zeros(size(llung));

	start_idx = rlung_idx(1);
	end_idx = llung_idx(1);

	llung_slice = llung(:,:,end_idx);
	llung_slice = imfill(llung_slice);
	llung_slice = imbinarize(llung_slice);
	llung_slice = bwareafilt(llung_slice, 1);

	% Get top-right most left lung
	[llungr, ~] = find(llung_slice == 1);
	[ltop_rightr, ~] = min(llungr);
    [~, ltop_rightc] = find(llung_slice(ltop_rightr,:) == 1, 1, 'last');

	% Get bottom-right most left lung
	[lbottom_rightr, ~] = max(llungr);
	[~, lbottom_rightc] = find(llung_slice(lbottom_rightr,:) == 1, 1, 'last');

	% Trace boundary
	lx = ltop_rightr; ly = ltop_rightc;
	Bl = bwtraceboundary(llung_slice, [lx ly], 'SE', 8, Inf, 'clockwise');	   
	interior_leftx = find(Bl(:,1) == lbottom_rightr);
	interior_lefty = find(Bl(interior_leftx,2) == lbottom_rightc);
	interior_lefty = interior_lefty + interior_leftx(interior_lefty);
	interior_boundaryl = Bl(1:interior_lefty,:);

	for z = start_idx:end_idx
		rlung_slice = rlung(:,:,end_idx);

		% Fill in holes
		rlung_slice = imfill(rlung_slice);
		rlung_slice = imbinarize(rlung_slice);
		rlung_slice = bwareafilt(rlung_slice, 1);

		% Get top-left most right lung
		[rlungr, ~] = find(rlung_slice == 1);
		[rtop_leftr, ~] = min(rlungr);
	    [~, rtop_leftc] = find(rlung_slice(rtop_leftr,:) == 1, 1);

		% Get bottom-left most right lung
		[rbottom_leftr, ~] = max(rlungr);
		[~, rbottom_leftc] = find(rlung_slice(rbottom_leftr,:) == 1, 1);

		% Traverse interior of right lung wall
		rx = rbottom_leftr; ry = rbottom_leftc;
		Br = bwtraceboundary(rlung_slice, [rx ry], 'NW', 8, Inf, 'clockwise');	   
		interior_rightx = find(Br(:,1) == rtop_leftr);
		interior_righty = find(Br(interior_rightx,2) == rtop_leftc);
		interior_righty = interior_righty + interior_rightx(interior_righty);
		interior_boundaryr = Br(1:interior_righty,:);

		for idx1 = 1:size(interior_boundaryr,1)
			for idx2 = 1:size(interior_boundaryl,1)
				[X, Y, Z] = bresenham_line3d([interior_boundaryr(idx1, 1) interior_boundaryr(idx1, 2) z], [interior_boundaryl(idx2, 1) interior_boundaryl(idx2, 2) end_idx]);
				linearIdxs = sub2ind(size(ms_w_top), X, Y, Z);
				ms_w_top(linearIdxs) = 1;
			end
		end
	end 
	for z = start_idx:end_idx
		ms_top_slice = ms_w_top(:,:,z);
	    [x, y] = find(ms_top_slice == 1);
	    k = boundary(x, y, 1);
	    unknown = zeros(size(ms_top_slice));
	    unknown(min(x(k)):max(x(k)), min(y(k)):max(y(k))) = 1;
	    unknown_indices = find(unknown == 1);
		[uix, uiy] = ind2sub(size(unknown), unknown_indices);
	    in = insidepoly(uix, uiy, x(k), y(k));
	    known = zeros(size(ms_top_slice));
	    known(unknown_indices) = in;
	    known = imfill(known);
	    ms_w_top(:,:,z) = known;
	end
end

%{ 
OLD CODE:
function [validRoute] = is_adjacent(x, y, border)
	if border(x + 1, y) == 0 || border(x, y + 1) == 0 || border(x - 1, y) == 0 || border(x, y - 1) == 0
		validRoute = true;
	else
		validRoute = false;
	end
end

function [validRoute] = is_diagonal(x, y, border)
	if border(x - 1, y - 1) == 0 || border(x - 1, y + 1) == 0 || border(x + 1, y - 1) == 0 || border(x + 1, y + 1) == 0
		validRoute = true;
	else
		validRoute = false;
	end
end
   %             if facing == 1
    %             	if llung_slice(lx, ly + 1) == 100 && is_adjacent(lx, ly + 1, llung_slice) % Right
    %             		ly = ly + 1;
    %         		elseif llung_slice(lx, ly + 1) == 100 && is_diagonal(lx, ly + 1, llung_slice)
    %         			ly = ly + 1;
    %         			temp_canvas(lx, ly) = 100;
    %         			lx = lx - 1;
    %         			facing = 4;
    %     			else
    %     				facing = 2;
    % 				end
				% elseif facing == 2
				% 	if llung_slice(lx + 1, ly) == 100 && is_adjacent(lx + 1, ly, llung_slice) % Down
				% 		lx = lx + 1;
				% 	elseif llung_slice(lx + 1, ly) == 100 && is_diagonal(lx + 1, ly, llung_slice)
    %         			lx = lx + 1;
    %         			temp_canvas(lx, ly) = 100;
    %         			ly = ly + 1;
    %         			facing = 1;
				% 	else
				% 		facing = 3;
				% 	end
				% elseif facing == 3
				% 	if llung_slice(lx, ly - 1) == 100 && is_adjacent(lx, ly - 1, llung_slice) % Left
				% 		ly = ly - 1;
				% 	elseif llung_slice(lx, ly - 1) == 100 && is_diagonal(lx, ly - 1, llung_slice)
    %         			ly = ly - 1;
    %         			temp_canvas(lx, ly) = 100;
    %         			lx = lx + 1;
    %         			facing = 2;
				% 	else
				% 		facing = 4;
				% 	end
				% elseif facing == 4
				% 	if llung_slice(lx - 1, ly) == 100 && is_adjacent(lx - 1, ly, llung_slice) % Up
				% 		lx = lx - 1;
				% 	elseif llung_slice(lx - 1, ly) == 100 && is_diagonal(lx - 1, ly, llung_slice)
    %         			lx = lx - 1;
    %         			temp_canvas(lx, ly) = 100;
    %         			ly = ly - 1;
    %         			facing = 3;
				% 	else
				% 		facing = 1;
				% 	end
				% end 



				% if llung_slice(lx, ly + 1) == 100 && temp_canvas(lx, ly + 1) ~= 100 && is_adjacent(lx, ly + 1, llung_slice) % Right
				% 	ly = ly + 1;
				% elseif llung_slice(lx + 1, ly + 1) == 100 && temp_canvas(lx + 1, ly + 1) ~= 100 && is_adjacent(lx + 1, ly + 1, llung_slice) % Down right
				% 	lx = lx + 1; ly = ly + 1;
				% elseif llung_slice(lx + 1, ly) == 100 && temp_canvas(lx + 1, ly) ~= 100 && is_adjacent(lx + 1, ly, llung_slice) % Down
				% 	lx = lx + 1;
				% elseif llung_slice(lx + 1, ly - 1) == 100 && temp_canvas(lx + 1, ly - 1) ~= 100 && is_adjacent(lx + 1, ly - 1, llung_slice) % Down left
				% 	lx = lx + 1; ly = ly - 1;
				% elseif llung_slice(lx, ly - 1) == 100 && temp_canvas(lx, ly - 1) ~= 100 && is_adjacent(lx, ly - 1, llung_slice) % Left
				% 	ly = ly - 1;
				% elseif llung_slice(lx - 1, ly - 1) == 100 && temp_canvas(lx - 1, ly - 1) ~= 100 && is_adjacent(lx - 1, ly - 1, llung_slice) % Up left
				% 	lx = lx - 1; ly = ly - 1;
				% elseif llung_slice(lx - 1, ly) == 100 && temp_canvas(lx - 1, ly) ~= 100 && is_adjacent(lx - 1, ly, llung_slice) % Up
				% 	lx = lx - 1;
				% elseif llung_slice(lx - 1, ly + 1) == 100 && temp_canvas(lx - 1, ly + 1) ~= 100 && is_adjacent(lx - 1, ly + 1, llung_slice) % Up right
				% 	lx = lx - 1; ly = ly + 1;
				% else
				% 	disp('Something is fucked!');
				% 	break;		
    %             end
			while lx ~= lbottom_rightr && ly ~= lbottom_rightc
				temp_canvas(lx, ly) = 100;
				if llung_slice(lx + 1, ly) == 1 && temp_canvas(lx + 1, ly) ~= 100 && adjacent_to_white(lx + 1, y, llung_slice) % Right
					lx = lx + 1;
				elseif llung_slice(lx + 1, ly + 1) == 100 && temp_canvas(lx + 1, ly + 1) ~= 100 && adjacent_to_white(lx + 1, ly + 1, llung_slice) % Down right
					lx = lx + 1; ly = ly + 1;
				elseif llung_slice(lx, ly + 1) == 100 && temp_canvas(lx, ly + 1) ~= 100 && adjacent_to_white(lx, ly + 1, llung_slice) % Down
					ly = ly + 1;
				elseif llung_slice(lx - 1, ly + 1) == 100 && temp_canvas(lx - 1, ly + 1) ~= 100 && adjacent_to_white(lx - 1, ly + 1, llung_slice) % Down left
					lx = lx - 1; ly = ly + 1;
				elseif llung_slice(lx - 1, ly) == 100 && temp_canvas(lx - 1, ly) ~= 100 && adjacent_to_white(lx - 1, ly, llung_slice) % Left
					lx = lx - 1;
				elseif llung_slice(lx - 1, ly - 1) == 100 && temp_canvas(lx - 1, ly - 1) ~= 100 && adjacent_to_white(lx - 1, ly - 1, llung_slice) % Up left
					lx = lx - 1; ly = ly - 1;
				elseif llung_slice(lx, ly - 1) == 100 && temp_canvas(lx, ly - 1) ~= 100 && adjacent_to_white(lx, ly - 1, llung_slice) % Up
					ly = ly - 1;
				elseif llung_slice(lx + 1, ly - 1) == 100 && temp_canvas(lx + 1, ly - 1) ~= 100 && adjacent_to_white(lx + 1, ly - 1, llung_slice) % Up right
					lx = lx + 1; ly = ly - 1;
				elseif lx == ltop_rightr && ly == ltop_rightc
					is_full_lung = false;
					temp_canvas(temp_canvas == 100) = 0;
					break;
				end
			end
% Pseudocode for calculating mediastinal volume
% Assume 3D binary matrix with left and right lung segmented
% lungs = matrix
for each slice
	get left lung region
	get right lung region

	get top right most left lung position
	get top left most right lung position

	save <- calculate top line boundary

	get bottom right most left lung position
	get bottom left most right lung position

	save <- calculate bottom line boundary

	save <- subset internal boundary of left lung bounded by top and bottom most left lung
	save <- subset internal boundary of right lung bounded by top and bottom most right lung
end
fill in boundaries set by interior left lung wall, interior right lung wall, calculated top, and calculated bottom

convert from binary to volume using metadata from DICOM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	k = find(ms == 1);
	[X1,Y1,Z1] = ind2sub(size(lungs), k);
    scatter3(X1,Y1,Z1,'MarkerEdgeColor','r');
    hold on;
    k = find(ms == 2);
	[X2,Y2,Z2] = ind2sub(size(lungs), k);
    scatter3(X2,Y2,Z2,'MarkerEdgeColor','b');

    % plot_all = find(ms == 1);
    % [plotx, ploty, plotz] = ind2sub(size(ms), plot_all);
    % P = [plotx ploty plotz];
    % volumeViewer(ms);
    % k = boundary(P,1);
    % trisurf(k,P(:,1),P(:,2),P(:,3),'Facecolor','red','FaceAlpha',0.1);

    % [X,Y] = meshgrid(1:512,1:512);
    % mesh(X,Y,ms_sum);
			% if ~isempty(potential_l) && ~isempty(potential_r)
			% 	ms(x,potential_l:potential_r,z) = 3;
			% end

    % [llt, llb, rlt, rlb] = get_all_lung_borders(lungs);
    % llt = rm_interp_outliers(llt);
    % llb = rm_interp_outliers(llb);
    % rlt = rm_interp_outliers(rlt);
    % rlb = rm_interp_outliers(rlb);    
    % % ms = fill_cloud(llt, llb, rlt, rlb, lungs);

    % scatter3(llt(:,2), llt(:,1), llt(:,3), 'MarkerEdgeColor', 'r');
    % % hold on;
    % scatter3(llb(:,2), llb(:,1), llb(:,3), 'MarkerEdgeColor', 'm');
    % scatter3(rlt(:,2), rlt(:,1), rlt(:,3), 'MarkerEdgeColor', 'b');
    % scatter3(rlb(:,2), rlb(:,1), rlb(:,3), 'MarkerEdgeColor', 'c');

function [ms] = fill_cloud(llt, llb, rlt, rlb, lungs)
	ms = zeros(size(lungs));
	nslices = size(llt,1);
	[xq, yq] = find(ms(:,:,1) == 0);

	% Only look at left lung
	llung = lungs;
	llung(llung ~= 1) = 0;
	rlung = lungs;
	rlung(rlung ~= 2) = 0;

	for nslice = 1:nslices
		x = [];
		y = [];
		llung_slice = llung(:,:,llt(nslice,3));
		for ind = llt(nslice,1):llb(nslice,1)
			c = find(llung_slice(:,ind) == 1, 1, 'last');
			if ~isempty(c)
				x = [x; c];
				y = [y; ind];
			end
		end
		rlung_slice = rlung(:,:,rlt(nslice,3));
		for ind = rlt(nslice,1):rlb(nslice,1)
			c = find(rlung_slice(:,ind) == 1, 1);
			if ~isempty(c)
				x = [x; c];
				y = [y; ind];
			end
		end
		x = [x; llt(nslice,2); llb(nslice,2); rlt(nslice,2); rlb(nslice,2)];
		y = [y; llt(nslice,1); llb(nslice,1); rlt(nslice,1); rlb(nslice,1)];
		k = boundary(x, y, 0);
		xv = x(k);
		yv = y(k);
		ms(:,:,llt(nslice,3)) = reshape(inpolygon(xq, yq, xv, yv), size(lungs,1), size(lungs,2));
	end
end

function [lung_border] = rm_interp_outliers(lung_border)
	% Distance from each point to every other point
	lung_dist = squareform(pdist(lung_border));
	mark_for_deletion = [];
    parfor ind = 1:length(lung_dist)
    	% Remove points with fewer than 3 points a distance of 40 away
    	if length(find(lung_dist(ind,:) < 40)) < 3
    		mark_for_deletion = [mark_for_deletion ind];
		end
	end

	% Replace points with inter/extra-polation
	lung_border_temp = lung_border;
	lung_border_temp(mark_for_deletion,:) = [];
	interpolated_x = round(interp1(lung_border_temp(:,3), lung_border_temp(:,2), lung_border(:,3) ,'pchip'));
	interpolated_y = round(interp1(lung_border_temp(:,3), lung_border_temp(:,1), lung_border(:,3), 'pchip'));
	lung_border(mark_for_deletion, :) = [interpolated_y(mark_for_deletion), interpolated_x(mark_for_deletion), lung_border(mark_for_deletion,3)];
end

function [llt, llb, rlt, rlb] = get_all_lung_borders(lungs)
	% Get image size for corrections
	lungs_size = size(lungs);

	% Get top left lung
	llt = get_llung_border(lungs, 'l');
	% Correct height
	llt(:,3) = lungs_size(3) - llt(:,3);

	% Get bottom left lung
	llb = get_llung_border(flipud(lungs), 'l');
	% Correct rows and height
	llb(:,1) = lungs_size(1) - llb(:,1);
	llb(:,3) = lungs_size(3) - llb(:,3);

	% Get top right lung
	rlt = get_llung_border(fliplr(lungs), 'r');
	% Correct columns and height
	rlt(:,2) = lungs_size(2) - rlt(:,2);
	rlt(:,3) = lungs_size(3) - rlt(:,3);

	% Get bottom right lung
	rlb = get_llung_border(rot90(lungs,2), 'r');
	% Correct columns and rows and height
	rlb(:,2) = lungs_size(2) - rlb(:,2);
	rlb(:,1) = lungs_size(1) - rlb(:,1);
	rlb(:,3) = lungs_size(3) - rlb(:,3);
end

function [tmll] = get_llung_border(lungs, l_or_r)
	% Which lung?
	if l_or_r == 'l'
		look_for = 1;
	elseif l_or_r == 'r'
		look_for = 2;
	else
		return;
	end

	% Only look at left lung
	llung = lungs;
	llung(llung ~= look_for) = 0;

	% Get coordinates of all left lung
	tmll_find = find(llung == look_for);
	[tmll_y, tmll_x, tmll_z] = ind2sub(size(llung), tmll_find);

	% Array of all top middle most left lung coordinates
	tmll = zeros(range(tmll_z), 3);
	min_slice = min(tmll_z);
	max_slice = max(tmll_z);

	% For all slices with left lung present
	for z_slice = min_slice:max_slice
		% Get all possible coordinates for current slice
		tmll_ind = find(tmll_z == z_slice);

		% First row containing 1
		y = min(tmll_y(tmll_ind));

		% Last column with 1 in first row containing 1
		x = find(llung(y,:,z_slice) == look_for, 1, 'last');

		% Store in top middle most left lung
		tmll(z_slice - min_slice + 1, :) = [y, x, z_slice];
	end
end
%     left_straight = zeros(size(lungs,2),2);
%     right_straight = zeros(size(lungs,2),2);
%     left_side = zeros(size(lungs,3),2);
%     right_side = zeros(size(lungs,3),2);
%     for y = 1:size(lungs,2)
%     	[ls_idxr, ls_idxc] = find(XZ_left(y,:) == 1, 1, 'last');
%     	if ~isempty(ls_idxr)
%     		left_straight(y,:) = [y, ls_idxc];
% 		end
% 		[rs_idxr, rs_idxc] = find(XZ_right(y,:) == 1, 1, 'last');
% 		if ~isempty(rs_idxr)
% 			right_straight(y,:) = [y, rs_idxc];
% 		end
% 	end
% 	for z = 1:size(lungs,3)
% 		[ls_idxr, ls_idxc] = find(XZ_left(:,z) == 1, 1, 'last');
%     	if ~isempty(ls_idxr)
%     		left_side(y,:) = [ls_idxr, z];
% 		end
% 		[rs_idxr, rs_idxc] = find(XZ_right(:,z) == 1, 1);
%     	if ~isempty(rs_idxr)
%     		right_side(z,:) = [rs_idxr, z];
% 		end
% 	end
% 	left_bound = setdiff(left_straight, left_side, 'rows');
% 	right_bound = setdiff(right_straight, right_side, 'rows');
%}