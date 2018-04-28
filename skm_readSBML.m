function [m_names, r_names, N] = skm_readSBML(SBML_file)
%
% This function enables the extraction of metabolites, reactions, and stoichiometric matrix from an SBML file for use in the SKM toolbox.
% 
% Requires the LibSBML package (including MATLAB API), which can be downloaded here: http://sbml.org/Software/libSBML


%%%%%%%%%%%%%%%%%%
% Read SBML file %
%%%%%%%%%%%%%%%%%%
SBML_struct = TranslateSBML(SBML_file, 0, 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve model components %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Metabolites:
n_met = size(SBML_struct.species,2);
m_names = cell(1, n_met);
for i = 1:n_met
    m_current = SBML_struct.species(i);
    m_id = m_current.id;
    m_names{i} = m_id;
end

% Reactions:
n_rct = size(SBML_struct.reaction,2);
r_names = cell(1, n_rct);
for i = 1:n_rct
    r_current = SBML_struct.reaction(i);
    r_id = r_current.id;
    r_names{i} = r_id;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain stoichiometric matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = zeros(n_met, n_rct);
for j=1:n_rct
    r_current = SBML_struct.reaction(j);
    
    substrates  = r_current.reactant;
    products    = r_current.product;
    n_subst = size(substrates,2);
    n_prod  = size(products,2);

    for s = 1:n_subst
        substrate_id = substrates(s).species;
        substrate_n  = substrates(s).stoichiometry;
        substrate_index = strcmp(m_names, substrate_id);
        N(substrate_index, j) = -substrate_n;
    end
    for p = 1:n_prod
        product_id = products(p).species;
        product_n  = products(p).stoichiometry;
        product_index = strcmp(m_names, product_id);
        N(product_index, j) = product_n;
    end
end
