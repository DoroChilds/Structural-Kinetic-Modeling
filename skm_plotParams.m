function figure_handles = skm_plotParams(modelParam_values, eigenvalues, stability, modelParam_names)
% Plotting the distributions of model parameters and resulting eigenvalues for stable 
% and unstable models.

% Extract model parameters:
params_s = modelParam_values.Enzyme_Substrates;
params_p = modelParam_values.Enzyme_Products;
params_r = modelParam_values.Regulators;
params_f = modelParam_values.FurtherParams;

% Check which types of parameters have been used and which parameter lists are empty:
plot_s = ~isempty(params_s);
plot_p = ~isempty(params_p);
plot_r = ~isempty(params_r);
plot_f = ~isempty(params_f);

% Create cells with parameter names:
theta_names = struct;
if isfield(modelParam_names, 'Enzyme_Substrates')
    theta_names.Enzyme_Substrates = strcat(modelParam_names.Enzyme_Substrates(:,1), '-',modelParam_names.Enzyme_Substrates(:,2));
else
    modelParam_names.Enzyme_Substrates = 1:size(params_s, 1);
end
if isfield(modelParam_names, 'Enzyme_Products')
    theta_names.Enzyme_Products = strcat(modelParam_names.Enzyme_Products(:,1), '-',modelParam_names.Enzyme_Products(:,2));
else
    modelParam_names.Enzyme_Products = 1:size(params_p, 1);
end
if isfield(modelParam_names, 'Regulators')
    theta_names.Regulators = strcat(modelParam_names.Regulators(:,1), '-',modelParam_names.Regulators(:,2));
else
    modelParam_names.Regulators = 1:size(params_r, 1);
end
if isfield(modelParam_names, 'FurtherParams')
    theta_names.FurtherParams = strcat(modelParam_names.FurtherParams(:,1), '-',modelParam_names.FurtherParams(:,2));
else
    modelParam_names.FurtherParams = 1:size(params_f, 1);
end


% Check if global parameters for substrates and products were applied (as in Steuer et al. (2006):
all_equal_random = 0;
if plot_s && plot_p
    if max(abs(std([params_s, params_p]')))<2*eps
        all_equal_random = 1;
    end
end

% Determine row indices of the stable, unstable, and undefined models in the parameter and eigenvalue matrices:
S_index = stability==1;
I_index = stability==0;
U_index = isnan(stability); % Models for which the stability could not be determined because they are too close to zero (marked as 'nan' in the labels vector):

% Check if there any models exhibiting each of the different types of dynamical properties:
plot_STABLE     = sum(S_index)>0;
plot_UNSTABLE   = sum(I_index)>0;
plot_UNCERTAIN  = sum(U_index)>0;

figure_handles = [];

%% Plot histogram of largest real parts of eigenvalues
h_temp = figure;
figure_handles = [figure_handles, h_temp];

plot_rows = 1 + plot_STABLE + plot_UNSTABLE + plot_UNCERTAIN;

% ALL models (stable and unstable):
subplot(plot_rows, 1, 1)
hist(real(eigenvalues(:,1)), 100);
title('Distribution of largest real parts of the eigenvalues', 'FontSize', 14)
xlabel('\Lambda^{max}_{R}', 'FontSize', 12)
set(gca, 'FontSize',12)

plot_iter = 2;
% STABLE models only:
if plot_STABLE
    subplot(plot_rows, 1, plot_iter)
    hist(real(eigenvalues(S_index,1)), 100)
    title('Distribution of negative largest real parts only', 'FontSize', 14)
    xlabel('\Lambda^{max, -}_{R}', 'FontSize', 12)
    set(gca, 'FontSize',12)
    plot_iter = plot_iter + 1;
end

% UNSTABLE models only:
if plot_UNSTABLE
    subplot(plot_rows, 1, plot_iter)
    hist(real(eigenvalues(I_index,1)), 100)
    title('Distribution of positive largest real parts only', 'FontSize', 14)
    xlabel('\Lambda^{max, +}_{R}', 'FontSize', 12)
    set(gca, 'FontSize',12)
    plot_iter = plot_iter + 1;
end

% UNCERTAIN models only:
if plot_UNCERTAIN
    subplot(plot_rows, 1, plot_iter)
    hist(real(eigenvalues(U_index,1)), 100)
    title('Distribution of largest real parts close to zero', 'FontSize', 14)
    xlabel('\Lambda^{max, 0}_{R}', 'FontSize', 12)
    set(gca, 'FontSize',12)
end


%% Plot distribution of parameters for stable, unstable, and undefined models

if ~ all_equal_random
    
    % Determine number of subfigures:
    plot_rows = plot_STABLE + plot_UNSTABLE + plot_UNCERTAIN;
    plot_cols = sum([plot_s, plot_p, plot_r, plot_f]);
    
    % Determine range of those model parameters for which the sampling intervals are not necessarily equal to [0,1]:
    if plot_r
        lim_r = [min(min(params_r))-0.1, max(max(params_r))+0.1];   % Regulatory parameters
    end
    if plot_f
        lim_f = [min(min(params_f))-0.1, max(max(params_f))+0.1];   % Further parameters of unspecified type
    end
    
    lim_pos = [-0.05,1.05];
    lim_neg = [-1.05,0.05];
    angle_labels = 45;
    
    h_temp = figure;
    figure_handles = [figure_handles, h_temp];

    fig_count = 1;
    
    %% Stable states:
    if plot_STABLE
        if plot_s
            subplot(plot_rows, plot_cols, fig_count)
            boxplot(params_s(S_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_pos)
            title(sprintf('Substrate parameters \n in %g stable models', sum(S_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.Enzyme_Substrates;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end
            hold on
            fig_count = fig_count+1;
        end
        
        if plot_p
            subplot(plot_rows,plot_cols,fig_count)
            boxplot(params_p(S_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_neg)
            title(sprintf('Product parameters \n in %g stable models', sum(S_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.Enzyme_Products;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end
            hold on
            fig_count = fig_count+1;
        end
        
        if plot_r
            subplot(plot_rows,plot_cols,fig_count)
            boxplot(params_r(S_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_r)
            title(sprintf('Regulatory parameters \n in %g stable models', sum(S_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.Regulators;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end
            hold on
            fig_count = fig_count+1;
        end
        
        if plot_f
            subplot(plot_rows,plot_cols,fig_count)
            boxplot(params_f(S_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_f)
            title(sprintf('Further parameters \n in %g stable models', sum(S_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.FurtherParams;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end
            hold on
            fig_count = fig_count+1;
        end
    end
    
    %% Unstable states:
    if plot_UNSTABLE
        if plot_s
            subplot(plot_rows,plot_cols, fig_count)
            boxplot(params_s(I_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_pos)
            title(sprintf('Substrate parameters \n in %g unstable models', sum(I_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.Enzyme_Substrates;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                    'FontSize', 12)
            end            
            hold on
            fig_count = fig_count+1;
        end
        
        if plot_p
            subplot(plot_rows,plot_cols,fig_count)
            boxplot(params_p(I_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_neg)
            title(sprintf('Product parameters \n in %g unstable models', sum(I_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.Enzyme_Products;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end            
            hold on
            fig_count = fig_count+1;
        end
        
        if plot_r
            subplot(plot_rows,plot_cols,fig_count)
            boxplot(params_r(I_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_r)
            title(sprintf('Regulatory parameters \n in %g unstable models', sum(I_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.Regulators;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end            
            hold on
            fig_count = fig_count+1;
        end
        
        if plot_f
            subplot(plot_rows,plot_cols,fig_count)
            boxplot(params_f(I_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_f)
            title(sprintf('Further parameters \n in %g unstable models', sum(I_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.FurtherParams;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end            
            hold on
            fig_count = fig_count+1;
        end
    end
    
    %% Undetermined states:
    if plot_UNCERTAIN
        if plot_s
            subplot(plot_rows,plot_cols, fig_count)
            boxplot(params_s(U_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_pos)
            title(sprintf('Substrate parameters \n in %g undefined models', sum(U_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.Enzyme_Substrates;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end            
            fig_count = fig_count+1;
        end
        
        if plot_p
            subplot(plot_rows,plot_cols,fig_count)
            boxplot(params_p(U_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_neg)
            title(sprintf('Product parameters \n in %g undefined models', sum(U_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.Enzyme_Products;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end            
            fig_count = fig_count+1;
        end
        
        if plot_r
            subplot(plot_rows,plot_cols,fig_count)
            boxplot(params_r(U_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_r)
            title(sprintf('Regulatory parameters \n in %g undefined models', sum(U_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.Regulators;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end            
            fig_count = fig_count+1;
        end
        
        if plot_f
            subplot(plot_rows,plot_cols,fig_count)
            boxplot(params_f(U_index,:), 'labelorientation','inline')
            set(gca, 'FontSize', 12)
            ylim(lim_f)
            title(sprintf('Further parameters \n in %g undefined models', sum(U_index)), 'FontSize', 14);
            text_handle = findobj(gca, 'Type', 'text');
            label_vec = theta_names.FurtherParams;
            for cnt = 1:length(text_handle)
                set(text_handle(cnt), 'Rotation', angle_labels, ...
                    'String', label_vec{length(label_vec)-cnt+1}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                    'FontSize', 12)
            end            
            fig_count = fig_count+1;
        end
    end
end

%% If all parameters are set to equal values per SK-model, plot the stability boundaries in 2D
if all_equal_random
    t_vec = params_s(:,1);
    h_temp = figure;
    figure_handles = [figure_handles, h_temp];

    hold on
    if ~isempty(params_p)
        p_vec = params_p(:,1);
        plot(t_vec(S_index), p_vec(S_index), '.', 'color', 'red', 'LineWidth', 1)
        plot(t_vec(I_index), p_vec(I_index), '.', 'color', 'blue', 'LineWidth', 1)
        xlabel('global substrate saturation \theta^{\mu}_S', 'FontSize', 12)
        ylabel('global product saturation \theta^{\mu}_P', 'FontSize', 12)
    else
        plot(t_vec(S_index), real(max_eigen(S_index,1)), '.', 'color', 'red', 'LineWidth', 1)
        plot(t_vec(S_index), imag(max_eigen(S_index,1)), '.', 'color', 'black', 'LineWidth', 1)
        plot(t_vec(S_index), imag(max_eigen(S_index,2)), '.', 'color', 'black', 'LineWidth', 1)
        plot(t_vec(I_index), real(max_eigen(I_index,1)), '.', 'color', 'blue', 'LineWidth', 1)
        plot(t_vec(I_index), imag(max_eigen(I_index,1)), '.', 'color', 'c', 'LineWidth', 1)
        plot(t_vec(I_index), imag(max_eigen(I_index,2)), '.', 'color', 'c', 'LineWidth', 1)
        xlabel('global substrate saturation \theta^{\mu}_S', 'FontSize', 12)
        ylabel('\Lambda^{real}_{max}', 'FontSize', 12)
    end
    hold off
    legend('Stable', 'Instable')
    title('Stable and instable regions', 'FontSize', 14)
    %     saveas(gcf, fullfile(dir_results,[filename_fixed,'-globalParams.png']), 'png')
end