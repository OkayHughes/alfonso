function [iter_nums, iprods] = iterate_iprod(ret)
iprods = zeros(size(ret.solnIterNum, 1)-2, 1);
iter_nums = cell2mat(ret.solnIterNum);
iter_nums = iter_nums(2:end-1);

for i=2:size(ret.solnIterNum, 1)-1
    soln_prev = ret.solns{i-1};
    soln_i = ret.solns{i};
    soln_next = ret.solns{i+1};
    iterate_prev =[soln_prev.x; soln_prev.s];
    iterate_i =[soln_i.x; soln_i.s];
    iterate_next =[soln_next.x; soln_next.s];
    to_prev = iterate_prev - iterate_i;
    to_next = iterate_next - iterate_i;
    iprods(i-1) = norm(to_prev'*to_next);
end
end