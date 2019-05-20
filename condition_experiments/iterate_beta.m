function [iter_nums, beta] = iterate_beta(ret)
iter_nums = 1:ret.nIterations;
beta = ret.betaPred;
end