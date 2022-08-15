# Input
Copy-number alterations data for ovarian tumors (TCGA) and ovarian cancer cell-lines (CCLE) obtained from cBioPortal. Mutation data for ovarian tumors for the most frequently mutated genes in ovarian cancer reported in the COSMIC database.

# Libraries
Uses the nnmf() function built-into MATLAB for consensus feature NNMF. Relevant study cited in the manuscript.

# Output
k: user-defined number of clusters.

W matrix: N x k matrix of (molecular) feature scores for each cluster, where N is the total number of molecular features

H matrix: k x P matrix of sample loadings, where P is the number of samples/patients/tumors.
