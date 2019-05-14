function area = areaHessEigs(hessian)
    eigenvs = eig(hessian);
    eigenvs = eigenvs/max(eigenvs) / size(hessian, 1);
    area = trapz(sort(eigenvs));
end