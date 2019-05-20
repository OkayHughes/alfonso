function ret = testPolyOpt(nvars, deg, scalefac)
degree = deg;

vars = msspoly('x', nvars);
costPoly = ones(1, nvars) * (vars + 1).^2 - 1;
ub = ones(nvars, 1);
lb = -1 * ones(nvars, 1);
polyWeights = (vars - ub) .* (lb - vars);
ub = ub * scalefac;
lb = lb * scalefac;

intParams = FeketeBasis(size(vars, 1), degree, 1);
intParams = cubeFromBasis(intParams);
ret = polyOpt(vars, costPoly, polyWeights, ub, lb, intParams);

hessIterNum = ret.solnIterNum;
hessCond = zeros(size(ret.solnIterNum));
hessEigArea = zeros(size(ret.solnIterNum));

for i=1:size(hessIterNum, 1)
    soln = ret.solns{i};
    hessCond(i) = cond(soln.H);
    hessEigArea = areaHessEigs(soln.H);
end

ret.hessEigArea = hessEigArea;
ret.hessCond = hessCond;

%scatter(cell2mat(ret.solnIterNum), ret.hessEigArea)

end
