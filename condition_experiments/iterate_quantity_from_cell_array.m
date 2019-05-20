function [iter_nums, dists] =  iterate_quantity_from_cell_array(ret_cell, iterate_fn)
iter_nums = cell(size(ret_cell));
dists = cell(size(ret_cell));
for i=1:size(ret_cell, 1)
    for j=1:size(ret_cell, 2)
        ret = ret_cell{i, j};
        [iter_nums{i, j}, dists{i, j}] = iterate_fn(ret);
    end
end