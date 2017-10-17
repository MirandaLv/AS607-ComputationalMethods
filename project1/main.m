% Author: Miranda Lv
% Date: 09/22/2017
% Course: AS607 
% Project 1

% In order to run this script, change to local folder directory
foldername = '/Users/miranda/Documents/Learning/computational_method/AS607/Matlab/project1_MirandaLv/outputs';


% Input functions
% tol = 0.001, 0.00001, 0.00000001
% fcta = @(x) exp(x) + 2^(-x) + 2*cos(x) - 6 1<=x<=2
% fctb = @(x) log(x-1) + cos(x-1) 1.3<=x<=2
% fctc = @(x) 2*x*cos(2*x) - (x-2)^2 2<=x<=3 and 3<=x<=4
% fctd = @(x) (x-2)^2 - log(x) 1<=x<=2, e<=x<=4
% fcte = @(x) exp(x) - 3*x^2  0<=x<=1 and 3<=x<=5
% fctf = @(x) sin(x) - exp(-x) 0<=x<=1, 3<=x<=4, 6<=x<=7


% call functions

% initial setting

legs = [0.001, 0.00001, 0.00000001];
legs_legs = ["TOL==0.001", "TOL==0.00001", "TOL==0.00000001"]
k = [3 5 8]

maxiter = 100

xlb = "iteration";
ylb = "error";

% function a: fcta = @(x) exp(x) + 2^(-x) + 2*cos(x) - 6 
% interval: 1<=x<=2

fcta = @(x) exp(x) + 2^(-x) + 2*cos(x) - 6;
figure;
fplot(fcta,[1 2]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function A')
xlabel('1<=x<=2');
ylabel('y: exp(x) + 2^(-x) + 2*cos(x) - 6;');
legend('Function A', 'Red line: y = 0');
fullFileName = fullfile(foldername,'functionA_chart');
saveas(gcf,fullFileName);

for i = 1:3
    leg = legs(i)
    figure;
    
    % bisection
    fcta_bisection = bisection_method(fcta,1,2,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fcta_bisection(:,2),'-o');
    hold on
    
    % newton
    fcta_newton = newton_method(fcta,1,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fcta_newton(:,2),'-o');
    hold on
    
    % secant
    fcta_secant = secant_method(fcta, 1, 2, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fcta_secant(:,2),'-o');
    
    title(['Graph of log-error of function a per iteration' legs_legs(i) '1<=x<=2']);
    xlabel('Iteration');
    ylabel('log-error');
    
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('fa_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3
        iterat = [fcta_bisection(end,1),fcta_newton(end,1), fcta_secant(end,1)];
        errors = [fcta_bisection(end,2),fcta_newton(end,2), fcta_secant(end,2)];
        roots = [fcta_bisection(end,3),fcta_newton(end,3), fcta_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tablea = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tablea);
    end
end



%function b: fctb = @(x) log(x-1) + cos(x-1)
%interval: 1.3<=x<=2

fctb = @(x) log(x-1) + cos(x-1);
figure;
fplot(fctb,[1 2]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function B')
xlabel('1.3<=x<=2');
ylabel('y: exp(x) + 2^(-x) + 2*cos(x) - 6;');
legend('Function B','Red line: y = 0');
fullFileName = fullfile(foldername,'functionB_chart');
saveas(gcf,fullFileName)

for i = 1:3
    leg = legs(i)

    figure;
    % bisection
    fctb_bisection = bisection_method(fctb,1.3,2,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctb_bisection(:,2),'-o');
    hold on
    

    % newton
    fctb_newton = newton_method(fctb,1.3,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctb_newton(:,2),'-o');
    hold on
    
    % secant
    fctb_secant = secant_method(fctb, 1.3, 2, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fctb_secant(:,2),'-o');
    
    title(['Graph of log-error of function b per iteration' legs_legs(i) '1.3<=x<=2']);
    xlabel('Iteration');
    ylabel('log-error');
    
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('fb_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3       
        iterat = [fctb_bisection(end,1),fctb_newton(end,1), fctb_secant(end,1)];
        errors = [fctb_bisection(end,2),fctb_newton(end,2), fctb_secant(end,2)];
        roots = [fctb_bisection(end,3),fctb_newton(end,3), fctb_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tableb = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tableb);
    end
    
    
    
end


% function c: fctc = @(x) 2*x*cos(2*x) - (x-2)^2
% interval: 2<=x<=3 and 3<=x<=4
fctc = @(x) 2*x*cos(2*x) - (x-2)^2;

figure;
fplot(fctc,[2 3]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function C: 2<=x<=3')
xlabel('2<=x<=3');
ylabel('y: 2*x*cos(2*x) - (x-2)^2;');
legend('Function C','Red line: y = 0');
fullFileName = fullfile(foldername,'functionC_chart_iter1');
saveas(gcf,fullFileName)

for i = 1:3
    leg = legs(i)
    
    figure;
    % bisection
    fctc_bisection = bisection_method(fctc,2,3,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctc_bisection(:,2),'-o');
    hold on

    % newton
    fctc_newton = newton_method(fctc,2,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctc_newton(:,2),'-o');
    hold on
    
    % secant
    fctc_secant = secant_method(fctc, 2, 3, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fctc_secant(:,2),'-o');
    
    title(['Graph of log-error of function c per iteration' legs_legs(i) '2<=x<=3']);
    xlabel('Iteration');
    ylabel('log-error');
    
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('fc_inter1_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3
        iterat = [fctc_bisection(end,1),fctc_newton(end,1), fctc_secant(end,1)];
        errors = [fctc_bisection(end,2),fctc_newton(end,2), fctc_secant(end,2)];
        roots = [fctc_bisection(end,3),fctc_newton(end,3), fctc_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tablec1 = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tablec1);
    end
    
end

%main function title, xlabel, ylabel,legend, export name
%bisection: f in function, range, title, savefilename

%function c: fctc = @(x) 2*x*cos(2*x) - (x-2)^2
%interval: 3<=x<=4

fctc = @(x) 2*x*cos(2*x) - (x-2)^2;

figure;
fplot(fctc,[3 4]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function C: 3<=x<=4')
xlabel('3<=x<=4');
ylabel('y: 2*x*cos(2*x) - (x-2)^2;');
legend('Function C','Red line: y = 0');
fullFileName = fullfile(foldername,'functionC_chart_iter2');
saveas(gcf,fullFileName)

for i = 1:3
    leg = legs(i)
    
    figure;
    % bisection
    fctc_bisection = bisection_method(fctc,3.5,4,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctc_bisection(:,2),'-o');
    hold on

    % newton
    fctc_newton = newton_method(fctc,3.5,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctc_newton(:,2),'-o');
    hold on
    
    % secant
    fctc_secant = secant_method(fctc, 3.5, 4, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fctc_secant(:,2),'-o');
    
    title(['Graph of log-error of function c per iteration' legs_legs(i) '3<=x<=4']);
    xlabel('Iteration');
    ylabel('log-error');
    
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('fc_inter2_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3
        iterat = [fctc_bisection(end,1),fctc_newton(end,1), fctc_secant(end,1)];
        errors = [fctc_bisection(end,2),fctc_newton(end,2), fctc_secant(end,2)];
        roots = [fctc_bisection(end,3),fctc_newton(end,3), fctc_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tablec2 = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tablec2);
    end
    
end




% function d: fctd = @(x) (x-2)^2 - log(x)
% interval:  1<=x<=2, e<=x<=4
fctd = @(x) (x-2)^2 - log(x);

figure;
fplot(fctd,[1 2]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function D: 1<=x<=2')
xlabel('1<=x<=2');
ylabel('y: @(x) (x-2)^2 - log(x);');
legend('Function D','Red line: y = 0');
fullFileName = fullfile(foldername,'functionD_chart_iter1');
saveas(gcf,fullFileName)

for i = 1:3
    leg = legs(i)
    
    figure;
    % bisection
    fctd_bisection = bisection_method(fctd,1,2,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctd_bisection(:,2),'-o');
    hold on

    % newton
    fctd_newton = newton_method(fctd,1,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctd_newton(:,2),'-o');
    hold on
    
    
    % secant
    fctd_secant = secant_method(fctd, 1, 2, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fctd_secant(:,2),'-o');
    
    title(['Graph of log-error of function d per iteration' legs_legs(i) '1<=x<=2']);
    xlabel('Iteration');
    ylabel('log-error');
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('fd_inter1_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3
        iterat = [fctd_bisection(end,1),fctd_newton(end,1), fctd_secant(end,1)];
        errors = [fctd_bisection(end,2),fctd_newton(end,2), fctd_secant(end,2)];
        roots = [fctd_bisection(end,3),fctd_newton(end,3), fctd_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tabled1 = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tabled1);
    end
    
end


% function d: fctd = @(x) (x-2)^2 - log(x)
% interval:  e<=x<=4
% fctd = @(x) (x-2)^2 - log(x);

figure;
fplot(fctd,[exp(1) 4]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function D: e<=x<=4')
xlabel('e<=x<=4');
ylabel('y: @(x) (x-2)^2 - log(x);');
legend('Function D','Red line: y = 0');
fullFileName = fullfile(foldername,'functionD_chart_iter2');
saveas(gcf,fullFileName)

for i = 1:3
    leg = legs(i)

    % bisection
    fctd_bisection = bisection_method(fctd,exp(1),4,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctd_bisection(:,2),'-o');
    hold on
    
    % newton
    fctd_newton = newton_method(fctd,exp(1),leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctd_newton(:,2),'-o');
    hold on
    
    
    % secant
    fctd_secant = secant_method(fctd, exp(1),4, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fctd_secant(:,2),'-o');
    
    title(['Graph of log-error of function d per iteration' legs_legs(i) 'e<=x<=4']);
    xlabel('Iteration');
    ylabel('log-error');
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('fd_inter2_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3
        iterat = [fctd_bisection(end,1),fctd_newton(end,1), fctd_secant(end,1)];
        errors = [fctd_bisection(end,2),fctd_newton(end,2), fctd_secant(end,2)];
        roots = [fctd_bisection(end,3),fctd_newton(end,3), fctd_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tabled2 = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tabled2);
    end
    
    
end


% function e: fcte = @(x) e^x - 3*x^2  
% interval: 0<=x<=1      % and 3<=x<=5
fcte = @(x) exp(x) - 3*(x^2);

figure;
fplot(fcte,[0 1]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function E: 0<=x<=1')
xlabel('0<=x<=1');
ylabel('y: (x) exp(x) - 3*x^2;');
legend('Function E','Red line: y = 0');
fullFileName = fullfile(foldername,'functionE_chart_iter1');
saveas(gcf,fullFileName)

for i = 1:3
    leg = legs(i)
    
    figure;
    % bisection
    fcte_bisection = bisection_method(fcte,0.5,1,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fcte_bisection(:,2),'-o');
    hold on

    % newton
    fcte_newton = newton_method(fcte,0.5,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fcte_newton(:,2),'-o');
    hold on
    
    
    % secant
    fcte_secant = secant_method(fcte, 0.5,1, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fcte_secant(:,2),'-o');
    
    title(['Graph of log-error of function e per iteration' legs_legs(i) '0<=x<=1']);
    xlabel('Iteration');
    ylabel('log-error');
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('fe_inter1_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3
        iterat = [fcte_bisection(end,1),fcte_newton(end,1), fcte_secant(end,1)];
        errors = [fcte_bisection(end,2),fcte_newton(end,2), fcte_secant(end,2)];
        roots = [fcte_bisection(end,3),fcte_newton(end,3), fcte_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tablee1 = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tablee1);
    end
    
    
end

%function e: fcte = @(x) e^x - 3*x^2  
%interval: 3<=x<=5

figure;
fplot(fcte,[3 5]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function E: 3<=x<=5')
xlabel('3<=x<=5');
ylabel('y: (x) exp(x) - 3*x^2;');
legend('Function E','Red line: y = 0');
fullFileName = fullfile(foldername,'functionE_chart_iter2');
saveas(gcf,fullFileName)

for i = 1:3
    leg = legs(i)

    % bisection
    fcte_bisection = bisection_method(fcte,3.5,5,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fcte_bisection(:,2),'-o');
    hold on
    

    % newton
    fcte_newton = newton_method(fcte,3.5,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fcte_newton(:,2),'-o');
    hold on
    
    
    % secant
    fcte_secant = secant_method(fcte, 3.5,5, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fcte_secant(:,2),'-o');
    
    title(['Graph of log-error of function e per iteration' legs_legs(i) '3<=x<=5']);
    xlabel('Iteration');
    ylabel('log-error');
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('fe_inter2_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3    
        iterat = [fcte_bisection(end,1),fcte_newton(end,1), fcte_secant(end,1)];
        errors = [fcte_bisection(end,2),fcte_newton(end,2), fcte_secant(end,2)];
        roots = [fcte_bisection(end,3),fcte_newton(end,3), fcte_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tablee2 = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tablee2);
    end

    
    
end


% fuction f: fctf = @(x) sin(x) - e^(-x) 
% intervals: 0<=x<=1   %3<=x<=4, 6<=x<=7
fctf = @(x) sin(x) - exp(-x);

figure;
fplot(fctf,[0 1]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function F: 0<=x<=1')
xlabel('0<=x<=1');
ylabel('y: @(x) sin(x) - e^(-x);');
legend('Function F','Red line: y = 0');
fullFileName = fullfile(foldername,'functionF_chart_iter1');
saveas(gcf,fullFileName)

for i = 1:3
    leg = legs(i)
    
    figure;
    % bisection
    fctf_bisection = bisection_method(fctf,0,1,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctf_bisection(:,2),'-o');
    hold on
    

    % newton
    fctf_newton = newton_method(fctf,0,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctf_newton(:,2),'-o');
    hold on
    
    % secant
    fctf_secant = secant_method(fctf, 0,1, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fctf_secant(:,2),'-o');
    
    title(['Graph of log-error of function f per iteration' legs_legs(i) '0<=x<=1']);
    xlabel('Iteration');
    ylabel('log-error');
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('ff_inter1_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3      
        iterat = [fctf_bisection(end,1),fctf_newton(end,1), fctf_secant(end,1)];
        errors = [fctf_bisection(end,2),fctf_newton(end,2), fctf_secant(end,2)];
        roots = [fctf_bisection(end,3),fctf_newton(end,3), fctf_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tablef1 = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tablef1);
    end

    
    
end

% fuction f: fctf = @(x) sin(x) - e^(-x) 
% intervals: 3<=x<=4

figure;
fplot(fctf,[3 4]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function F: 3<=x<=4')
xlabel('3<=x<=4');
ylabel('y: @(x) sin(x) - e^(-x);');
legend('Function F','Red line: y = 0');
fullFileName = fullfile(foldername,'functionF_chart_iter2');
saveas(gcf,fullFileName)

for i = 1:3
    leg = legs(i)
    
    figure;
    % bisection
    fctf_bisection = bisection_method(fctf,3,4,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctf_bisection(:,2),'-o');
    hold on

    % newton
    fctf_newton = newton_method(fctf,3,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctf_newton(:,2),'-o');
    hold on
    
    % secant
    fctf_secant = secant_method(fctf, 3,4, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fctf_secant(:,2),'-o');
    
    title(['Graph of log-error of function f per iteration' legs_legs(i) '3<=x<=4']);
    xlabel('Iteration');
    ylabel('log-error');
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('ff_inter2_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3        
        iterat = [fctf_bisection(end,1),fctf_newton(end,1), fctf_secant(end,1)];
        errors = [fctf_bisection(end,2),fctf_newton(end,2), fctf_secant(end,2)];
        roots = [fctf_bisection(end,3),fctf_newton(end,3), fctf_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tablef2 = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tablef2);
    end
    
    
end


% fuction f: fctf = @(x) sin(x) - e^(-x) 
% intervals: 6<=x<=7

figure;
fplot(fctf,[6 7]);
hline = refline([0 0]);
hline.Color = 'r';
title('Function F: 6<=x<=7')
xlabel('6<=x<=7');
ylabel('y: @(x) sin(x) - e^(-x);');
legend('Function F','Red line: y = 0');
fullFileName = fullfile(foldername,'functionF_chart_iter3');
saveas(gcf,fullFileName)

for i = 1:3
    leg = legs(i)

    % bisection
    fctf_bisection = bisection_method(fctf,6,7,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctf_bisection(:,2),'-o');
    hold on
    % newton
    fctf_newton = newton_method(fctf,6,leg,maxiter);
    set(gca, 'YScale', 'log');
    plot(fctf_newton(:,2),'-o');
    hold on
    
    % secant
    fctf_secant = secant_method(fctf, 6,7, leg, maxiter);
    set(gca, 'YScale', 'log');
    plot(fctf_secant(:,2),'-o');
    
    title(['Graph of log-error of function f per iteration' legs_legs(i) '6<=x<=7']);
    xlabel('Iteration');
    ylabel('log-error');
    
    legend('Bisection', 'Newton', 'Secant')
    hold off
    
    baseFilename = sprintf('ff_inter3_tol-%d.png',k(i));
    fullFileName = fullfile(foldername,baseFilename);
    saveas(gcf,fullFileName);
    
    if i == 3      
        iterat = [fctf_bisection(end,1),fctf_newton(end,1), fctf_secant(end,1)];
        errors = [fctf_bisection(end,2),fctf_newton(end,2), fctf_secant(end,2)];
        roots = [fctf_bisection(end,3),fctf_newton(end,3), fctf_secant(end,3)];
        %rownames = ["iteration","error", "root"];
        combine = [iterat; errors; roots];
        tablef3 = array2table(combine, 'VariableNames', {'Bisection', 'Newton', 'Secant'})
        writetable(tablef3);
    end
    
    
end









