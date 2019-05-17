function ret = testPolyOpt(frac)
degree = 6;

vars = msspoly('x', 3);
costPoly = vars(1)^5 + vars(2)^2 - vars(3)^1 + 1;
ub = [1; 1; 1];
lb = [-1; -1; -1];
polyWeights = (vars - ub) .* (lb - vars);

intParams = FeketeBasis(size(vars, 1), degree, 1);
intParams = cubeFromBasis(intParams);
ret = polyOpt(vars, costPoly, polyWeights, ub, lb, intParams)

scatter(cell2mat(ret.hessIterNum), cell2mat(ret.hessEigArea))

end