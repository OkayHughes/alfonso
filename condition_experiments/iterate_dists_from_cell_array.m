function [iter_nums, dists] =  iterate_dists_from_cell_array(ret_cell)
iter_nums = cell(size(ret_cell));
dists = cell(size(ret_cell));
for i=1:size(ret_cell, 1)
    for j=1:size(ret_cell, 2)
        ret = ret_cell{i, j};
        iter_nums{i, j} = cell2mat(ret.solnIterNum);
        dists{i, j} = iterate_dists(ret);
    end
end