function [iter_nums, dists] = iterate_dists(ret)
final_soln = ret.solns{end};
final_iterate = [final_soln.x; final_soln.s];
dists = zeros(size(ret.solnIterNum));
for i=1:size(ret.solnIterNum, 1)
    soln_i = ret.solns{i};
    iterate_i =[soln_i.x; soln_i.s];
    dists(i) = norm(iterate_i - final_iterate);
end
iter_nums = cell2mat(ret.solnIterNum);
end