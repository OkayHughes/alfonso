



cond_res = load("condition_scale_pts.mat");
poly_opt_res = cond_res.scalingPolyResults;

for nvar=1:5
    figure(); hold on
    title(sprintf("n = %d", nvar))
    xlabel("Iteration Number")
    axis([0, 200, 0, 50])
    ylabel("log(cond(Hessian))")
    for distort_fac=1:5
        res = poly_opt_res{nvar, distort_fac};
        hess_iter_num = cell2mat(res.hessIterNum);
        hess_eig_area = log(cell2mat(res.hessCond));

        plot(hess_iter_num, hess_eig_area, "-o", 'LineStyle', 'none', 'DisplayName',sprintf('distortion = %0.2f', 0.2*distort_fac))
    end
    legend()
    saveas(gcf, sprintf("plots/cond_scale_nvar=%d.png", nvar))
end
