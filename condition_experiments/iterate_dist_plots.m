function iterate_dist_plots(rets, labels, figname)
[iter_nums, dists] = iterate_dists_from_cell_array(rets);
figure(); hold on;
title("Distance from iterate $i$ to final iterate");
xlabel("iterate");
ylabel("distance in L^2 norm")
for i=1:size(labels, 1)
    
    plot(iter_nums{i}, dists{i}, "-o", 'LineStyle', 'none', 'DisplayName',labels{i})

end
legend()

if exist("figname", 'var')
saveas(gcf, sprintf("plots/%s.png", figname));
end

end