# Define common ancestry groups

According to the phylogenetic tree of BA.2, we defined common ancestry groups of the BA.2 variants bearing mutations at position 452 in S as the follow procedures. First, the ancestral state of the amino acid at position 452 in S at each node was estimated using a fixed-rates continuous-time Markov model (Mk model) implemented in the R package “castor” (Louca and Doebeli, 2018). As a type of transition matrix in the Mk model, all rate different (ARD) matrix was selected. Second, we identified the branches connecting the parental-state (L) nodes and the mutated-sate (R, Q, or M) nodes. In these branches, it is expected that the mutation acquisitions in the S L452 residue occurred. Finally, we counted the descendant sequences of respective branches where the mutations in the S L452 were likely acquired. If the number of descendants is ≥10, we defined these descendant sequences as a common ancestry group of the BA.2 variants, which bears a mutation at position 452 in S.

## Contents:
*  **common_ancestry_groups.R**: scripts for tree construction
