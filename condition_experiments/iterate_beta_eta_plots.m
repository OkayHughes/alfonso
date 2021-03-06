function iterate_beta_eta_plots(rets, labels, figname)

[iter_nums, iprods] = iterate_quantity_from_cell_array(rets, @iterate_eta);
figure(); hold on;
title("Distance from iterate $i$ to final iterate");
xlabel("iterate");
ylabel("\eta")
for i=1:size(labels, 1)
    plot(iter_nums{i}, iprods{i}, "-o", 'LineStyle', 'none', 'DisplayName',labels{i})

end
legend()

if exist("figname", 'var')
saveas(gcf, sprintf("plots/%s.png", figname));
end

end