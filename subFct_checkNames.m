function [m_names, r_names] = subFct_checkNames(options, n_met, n_rct)
% If provided by the user, extract names for metabolites and reactions and check them for validity.
% Otherwise create names generically.

% Metabolite names:
generic_metabolite_names = 1;
if isfield(options, 'm_names')
    if ~isempty(options.m_names)
        if length(options.m_names) == n_met && iscell(options.m_names)
            m_names = options.m_names;
            generic_metabolite_names = 0;
        else
            if options.verbose
                fprintf('Warning: Given metabolite names cannot be used (wrong format or number). Assigning generic names instead.\n')
            end
        end
    end
end

if generic_metabolite_names
    m_names = cell(1, n_met);
    for i=1:n_met
        m_names{i} = ['S', num2str(i)];
    end
end

% Reaction names:
generic_reaction_names = 1;
if isfield(options, 'r_names')
    if ~isempty(options.r_names)
        if length(options.r_names) == n_rct && iscell(options.r_names)
            r_names = options.r_names;
            generic_reaction_names = 0;
        else
            if options.verbose
                fprintf('Warning: Given reaction names cannot be used (wrong format or number). Assigning generic names instead.\n')
            end
        end
    end
end
if generic_reaction_names
    r_names = cell(1, n_rct);
    for i=1:n_rct
        r_names{i} = ['v', num2str(i)];
    end
end


