function p = plot_method(inv, tit, xlb, ylb, leg)
    plot(inv(1,:), inv(2,:));
    title(tit);
    xlabel(xlb);
    ylabel(ylb);
    smallerror = inv(2,end)
    hline = refline([0 smallerror]);
    hline.Color = 'r';
    legsmall = sprintf('error = %d',smallerror)
    legend(leg,legsmall);
end