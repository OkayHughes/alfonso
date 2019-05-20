function [iter_nums, eta] = iterate_eta(ret)
iter_nums = 1:ret.nIterations;
eta = ret.etaCorr;
end