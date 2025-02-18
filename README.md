# Unveiling Complex Patterns: An Information-Theoretic Approach to High-Order Behaviors in Microarray Data
## Antonio Lacalamita, Alfonso Monaco*, Grazia Serino, Daniele Marinazzo, Nicola Amoroso, Loredana Bellantuono, Marianna La Rocca, Tommaso Maggipinto, Ester Pantaleo, Emanuele Piccinno, Viviana Scalavino, Sabina Tangaro, Gianluigi Giannelli, Sebastiano Stramaglia, Roberto Bellotti.

This is a codebase behind the analysis of the paper: "Unveiling Complex Patterns: An Information-Theoretic Approach to High-Order Behaviors in Microarray Data".
There are two different code examples regarding two different communities: one for HCC and one for ASD.
The Matlab script takes as input the initial community and creates the Synergy matrices (note that you have to compile this code with both the original matrix and the transposed matrix so that you can perform the analysis on both genes and subjects).
The R script analyzes the Synergy data in order to produce the differentially expressed genes or biologically enriched functions.
