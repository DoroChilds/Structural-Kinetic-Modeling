function options = subFct_checkOptions(options)

%% Default values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for model structure: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(options, 'm_names');    options.m_names   = cell(0);   end
if ~isfield(options, 'r_names');    options.r_names   = cell(0);   end
if ~isfield(options, 'm_dependent'); options.m_dependent   = [];   end
if ~isfield(options, 'm_exclude');  options.m_exclude = cell(0);    end     % metabolites to remove
if ~isfield(options, 'm_include');  options.m_include = cell(0);    end     % metabolites to include (the rest will be removed)
if ~isfield(options, 'r_exclude');  options.r_exclude = cell(0);    end     % reactions to remove
if ~isfield(options, 'r_include');  options.r_include = cell(0);    end     % reactions to remove (the rest will be removed)
if ~isfield(options, 'rm_zero_met');    options.rm_zero_met = false;    end     % remove metabolites with concentration = 0
if ~isfield(options, 'rm_zero_rct');    options.rm_zero_rct = false;    end     % remove reactions with velocity = 0
if ~isfield(options, 'rm_zero_N');      options.rm_zero_N   = false;    end     % remove all metabolites and reactions, which only have zero rows/columns in N


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for model parameters: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type of kinetics to be used when assigning SK-model parameters. 
% (choices: 'enzymatic_irrev', 'enzymatic_rev', 'massAction', default: 'enzymatic_irrev')
if ~isfield(options, 'kinetics');   options.kinetics  = 'enzymatic_irrev'; end
default_function = @(T, struct_indices, options) T;

% Handle to a function that performs user-defined postprocessing on the parameter matrix 
% before using it for computing the Jacobian matrix 
% (see example script ExampleScript_Girbig2012_CalvinCycle})
if ~isfield(options, 'fct_modify_params');    options.fct_modify_params = default_function;  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for program output: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print output to screen
if ~isfield(options, 'verbose');        options.verbose     = true;     end

% Plot parameter- and eigenvalue distributions
if ~isfield(options, 'plot');           options.plot        = true;     end

% Create balanced numbers of stable and unstable states. 
% Can be useful for machine-learning. (choices: TRUE/FALSE, default: FALSE)
if ~isfield(options, 'balance_stability');  options.balance_stability    = false;    end     

% Compute the full set of eigenvalues using the 'eig' function (if FALSE: compute only the 
% value with maximum real part using the faster 'eigs' function). (choices: TRUE/FALSE, default: TRUE)
if ~isfield(options, 'return_all_eigvals');  options.return_all_eigvals  = true;    end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include in future version: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sampling information
% if ~isfield(options, 'simType');    options.simType = 'MonteCarlo'; end


%% Kinetics
% Check if kinetics field contains a known value:
if ~sum(ismember({'massAction', 'enzymatic_irrev', 'enzymatic_rev'}, options.kinetics))
    %     fprintf('Warning: user-defined kinetics ("%s") not known.\nIf not already provided by the user, missing fields in theParamIntervals -\nstruct will be filled using irreversible enzyme kinetics.\n', options.kinetics)
    fprintf('Warning: user-defined kinetics (options.kinetics = "%s") not known.\n', options.kinetics)
    options.kinetics = 'enzymatic_irrev';
    fprintf('The setting options.kinetics = "%s" will be used instead.\n\n', options.kinetics)
end

%% Included or exluded components
% Allow only user-input about either EXCLUDED or INCLUDED metabolites, but not both:
if ~ ( isempty(options.m_exclude) || isempty(options.m_include) )
    fprintf('Warning: options field "exclude_met" and "include_met" both contain entries.\n\t\tPlease restrict your input to only one of these fields in order to avoid conflicts.\n\t\tContinuing with ALL metabolites for now...\n')
end

% Allow only user-input about either EXCLUDED or INCLUDED reactions, but not both:
if ~ ( isempty(options.r_exclude) || isempty(options.r_include) )
    fprintf('Warning: options field "exclude_rct" and "include_rct" both contain entries.\n\t\tPlease restrict your input to only one of these fields in order to avoid conflicts.\n\t\tContinuing with ALL reactions for now...\n')
end