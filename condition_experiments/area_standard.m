



cond_res = load("area_standard.mat");
poly_opt_res = cond_res.polyOptResults;

for nvar=1:5
    figure(); hold on
    title(sprintf("n = %d", nvar))
    xlabel("Iteration Number")
    axis([0, 40, 0, 1])
    ylabel("Area Under Scree Plot of Hessian")
    for hdeg=1:5
        res = poly_opt_res{nvar, hdeg};
        hess_iter_num = cell2mat(res.hessIterNum);
        hess_eig_area = cell2mat(res.hessEigArea);

        plot(hess_iter_num, hess_eig_area, "-o", 'LineStyle', 'none', 'DisplayName',sprintf('d = %d', hdeg*2))
    end
    legend()
    saveas(gcf, sprintf("plots/area_standard_nvar=%d.png", nvar))
end
