clear
close all

addpath(genpath('../utilities'))

% make path to write figures
figPath = '../../fig/memory_inference/';
mkdir(figPath)

% load 
load('mem_data_0.mat')
mem_data_0 = mem_data;
clear mem_data
load('mem_data_1.mat')
mem_data_1 = mem_data;
clear mem_data
%% Make figures 
n_plot = 1; % plot best N for each mem value

close all
bkg_color = [228,221,209]/255;
err_sig = 1;
lb_prct = 100*normcdf(-err_sig,0,1);
ub_prct = 100*normcdf(err_sig,0,1);

cmap_bu = brewermap(n_plot+2,'Blues');
cmap_gra = brewermap(n_plot+2,'Greys');
cmap_pu = brewermap(n_plot+2,'Purples');
cmap_gr = brewermap(n_plot+2,'Greens');
cmap_rd = brewermap(n_plot+2,'Reds');
color_cell = {cmap_gra cmap_pu cmap_gr cmap_bu cmap_rd};

% v and A inference over time
mem_vec = mem_data_1.mem_vec;


for temp = 0:1
    if ~temp
        mem_data = mem_data_0;
    else
        mem_data = mem_data_1;
    end
    step_axis = 1:size(mem_data.mem_array,1);
    mem_fig = figure;
    hold on
    for m = 1:length(mem_vec)
        [~, si] = sort(mem_data.logL_array(m,:), 'descend');        
        p1 = plot(step_axis, repelem(mem_vec(m),length(step_axis)),'-.','Color',color_cell{m}(end,:),'LineWidth',2);
        for n = 1:n_plot
            p2 = plot(step_axis, mem_data.mem_array(:,si(n),m),'-','Color',color_cell{m}(n+1,:),'LineWidth',2);
        end
    end

    xlabel('inference step')
    ylabel('memory (w)')

    legend([p1 p2], 'ground truth', 'model inference','Location','southeast')
    set(gca,'Fontsize',14)
    ylim([2 12])
    xlim([0 max(step_axis)])
    set(gca,'Color',bkg_color) 
    box on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    mem_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');

    saveas(mem_fig,[figPath 'mem_fig_' num2str(temp) '.png'])
    saveas(mem_fig,[figPath 'mem_fig_' num2str(temp) '.pdf'])
end

%% Make temperature grad plot
close all

stepVec = 1:750;
sg = 250;
gradientVec = 2*exp(-stepVec/sg*5)+1;

        
grad_fig = figure;
cmap = brewermap([],'Set2');
hold on
plot(stepVec, gradientVec, 'Color', cmap(2,:), 'LineWidth', 3)
plot(stepVec, repelem(1,length(gradientVec)), '-.k', 'LineWidth', 3)

xlabel('inference step')
ylabel('tempurature (\sigma^*/\sigma)')

legend('noise-based annealing','standard inference')
set(gca,'Fontsize',14)
ylim([0 3])
xlim([0 750])
% xlim([0 max(step_axis)])
set(gca,'Color',bkg_color) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

grad_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(grad_fig,[figPath 'temp_fig.png'])
saveas(grad_fig,[figPath 'temp_fig.pdf'])