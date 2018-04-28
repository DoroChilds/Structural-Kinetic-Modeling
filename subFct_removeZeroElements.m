function [N, S0, v0, m_dependent, param_intervals, m_names_left, r_names_left] = subFct_removeZeroElements(N, S0, v0, m_dependent, param_intervals, m_names, r_names, rm_zero_met, rm_zero_rct, rm_zero_N, verbose, kinetics)
% Rewrite function so that it calls the subfunctions already existing for removing reactions or metabolites.

%% Determine which metabolites or reactions should be removed and store indices in a boolean array
m_exclude = false(1, length(m_names));
r_exclude = false(1, length(r_names));

% --------------------------- %
% Zero metabolites (optional) %
% --------------------------- %
if rm_zero_met;
    m_exclude(S0==0) = 1;
end

% ------------------------- %
% Zero reactions (optional) %
% ------------------------- %
if rm_zero_rct;
    if strcmp(kinetics, 'enzymatic_rev')
        zero_rct = v0==0;
        zero_forward_rct = zero_rct(1:2:end-1);
        zero_backward_rct = zero_rct(2:2:end);
        zero_both_rct = find(zero_forward_rct & zero_backward_rct);
        r_exclude = false(size(v0));
        r_exclude([2*zero_both_rct-1, 2*zero_both_rct]) = true;
    else
        r_exclude(v0==0) = 1;
    end
end

% --------------------------------- %
% Zero rows/columns in N (optional) %
% --------------------------------- %
if rm_zero_N
    remove_cols = sum(abs(N), 1) == 0;
    if strcmp(kinetics, 'enzymatic_rev')        
        zero_forward_cols = remove_cols(1:2:end-1);
        zero_backward_cols = remove_cols(2:2:end);
        zero_both_cols = find(zero_forward_cols & zero_backward_cols);
        remove_cols = false(1, size(N, 2));
        remove_cols([2*zero_both_cols-1, 2*zero_both_cols]) = true;  
    end
    remove_rows = sum(abs(N), 2) == 0;
    
    r_exclude(remove_cols) = 1;
    m_exclude(remove_rows) = 1;
end

if verbose
    if sum(m_exclude)
        fprintf('Warning: The stoichiometric matrix contains empty rows!\nRemoving the following rows:\n');
        disp(char(m_names(m_exclude)));
        fprintf('With row index:\n');
        disp(num2str(find(m_exclude)));
        fprintf('\n');
    end
    
    if sum(r_exclude)
        fprintf('Warning: The stoichiometric matrix contains empty columns!\nRemoving the following columns:\n');
        disp(char(r_names(r_exclude)));
        fprintf('With column index:\n');
        disp(num2str(find(r_exclude)));
        fprintf('\n');
    end
end

%% Remove metabolites and reactions from the model
[N, S0, m_dependent, m_names_left, param_intervals] = subFct_removeMetabolites(m_exclude, [], N, S0, m_dependent, param_intervals, m_names, verbose);
[N, v0, r_names_left, param_intervals]              = subFct_removeReactions(r_exclude,   [], N, v0, param_intervals, r_names, verbose);
