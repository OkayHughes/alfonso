function ret = testPolyOpt(nvars, deg, scalefac)
degree = deg;

vars = msspoly('x', nvars);
costPoly = ones(1, nvars) * (vars - 1).^2 - 1;
ub = ones(nvars, 1);
lb = -1 * ones(nvars, 1);
polyWeights = (vars - ub) .* (lb - vars);
ub = ub * scalefac;

intParams = FeketeBasis(size(vars, 1), degree, 1);
intParams = cubeFromBasis(intParams);
ret = polyOpt(vars, costPoly, polyWeights, ub, lb, intParams)

scatter(cell2mat(ret.hessIterNum), cell2mat(ret.hessEigArea))

end
