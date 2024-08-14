lys_names = readlines("./gene_labels/GO_lys_strict.csv");
ubq_names = readlines("./gene_labels/GO_ubq_strict.csv");
all_names = readlines("./gene_labels/all_annotated_AHBA_GO_genes.csv");

risk_names = setxor(all_names, lys_names);
risk_names = setxor(risk_names, ubq_names);

lys_bool = ismember(lys_names, gene_expr.Properties.VariableNames);
valid_lys = lys_names(lys_bool);
lys_genes = gene_expr(:, valid_lys);

ubq_bool = ismember(ubq_names, gene_expr.Properties.VariableNames);
valid_ubq = ubq_names(ubq_bool);
ubq_genes = gene_expr(:, valid_ubq);

risk_bool = ismember(risk_names, gene_expr.Properties.VariableNames);
valid_risk = risk_names(risk_bool);
risk_genes = gene_expr(:, valid_risk);

lys_mat = table2array(lys_genes);
ubq_mat = table2array(ubq_genes);
risk_mat = table2array(risk_genes);

[lys_coeff, lys_PC] = pca(lys_mat);
[ubq_coeff, ubq_PC] = pca(ubq_mat);
[risk_coeff, risk_PC] = pca(risk_mat);

lys_expr = array2table(lys_PC);
ubq_expr = array2table(ubq_PC);
risk_expr = array2table(risk_PC);
gene_expr = [risk_expr, lys_expr, ubq_expr];

