function ans = testPolyOpt(frac)
degree = 20;

vars = msspoly('x', 2);
costPoly = vars(1)^2 + vars(2)^2 - 1;
ub = [1; 1];
lb = [-1; -1];
polyWeights = (vars - ub) .* (lb - vars);

intParams = FeketeCubeRand(size(vars, 1), degree, frac)
ans = polyOpt(vars, costPoly, polyWeights, ub, lb, intParams)

scatter(cell2mat(ans.hessIterNum), cell2mat(ans.hessEigArea))

end