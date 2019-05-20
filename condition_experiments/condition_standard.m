



cond_res = load("condition_standard.mat");
poly_opt_res = cond_res.polyOptResults;

for nvar=1:5
    figure(); hold on
    title(sprintf("n = %d", nvar))
    xlabel("Iteration Number")
    axis([0, 60, 0, 50])
    ylabel("log(cond(Hessian))")
    for hdeg=1:5
        res = poly_opt_res{nvar, hdeg};
        hess_iter_num = cell2mat(res.hessIterNum);
        hess_eig_area = log(cell2mat(res.hessCond));

        plot(hess_iter_num, hess_eig_area, "-o", 'LineStyle', 'none', 'DisplayName',sprintf('d = %d', hdeg*2))
    end
    legend()
    saveas(gcf, sprintf("plots/cond_standard_nvar=%d.png", nvar))
end
