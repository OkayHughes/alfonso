
cond_res = load("area_scale_pts.mat");
poly_opt_res = cond_res.scalingPolyResults;

for nvar=1:5
    figure(); hold on
    title(sprintf("n = %d, d = 8", nvar))
    xlabel("Iteration Number")
    axis([0, 120, 0, 1])
    ylabel("Area Under Scree Plot of Hessian")
    for distort_fac=1:5
        res = poly_opt_res{nvar, distort_fac};
        hess_iter_num = cell2mat(res.hessIterNum);
        hess_eig_area = cell2mat(res.hessEigArea);

        plot(hess_iter_num, hess_eig_area, "-o", 'LineStyle', 'none', 'DisplayName',sprintf('distortion factor = %0.2f', 0.2 * distort_fac))
    end
    legend()
    saveas(gcf, sprintf("plots/scale_nvar=%d.png", nvar))
end