function mask_fill          = smallholefill(mask, num_px, threshold, eq)
% Fill small holes in a mask.
% 
% Joe MacGregor (NASA)
% Last updated: 11/03/20

mask_fill                   = mask; % copy initial mask
mask_tmp1                   = false(size(mask)); % logical nowhere start
eval(['mask_tmp1(mask_fill ' eq ' threshold) = true;']) % start at threshold
mask_tmp2                   = imfill(mask_tmp1, 'holes'); % fill holes
holes_all                   = ~mask_tmp1 & mask_tmp2; % where holes were filled that weren't there before
holes_small                 = holes_all & ~bwareaopen(holes_all, num_px); % location of small holes only
mask_fill(holes_small)      = threshold;