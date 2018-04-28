function struct_indices = subFct_names_to_indices(m_names, r_names)

struct_m = struct();
struct_r = struct();
struct_indices = struct();

for i=1:length(m_names)
    m_temp = m_names{i};
    struct_m.(m_temp) = i;
end

for i=1:length(r_names)
    r_temp = r_names{i};
    struct_r.(r_temp) = i;
end

struct_indices.metabolites  = struct_m;
struct_indices.reactions    = struct_r;
