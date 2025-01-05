# Zero-Knowledge Proofs for Isogeny Walks
This repository contains the code from the paper "Verifying Isogeny Walks with Succinct Zero-Knowledge Arguments" by Bootle, de Feo, Spjeldnæs [2025].

The code is organized into the following modules:
- [isogeny_walk.txt](isogeny_walk.txt) contains the isogeny walk we generated for the proofs
- [src/generalized_fri.rs](src/generalized_fri.rs) is a module for the generalized FRI PCS introduced in the paper. This can be used in external projects as is.
- [src/isogeny_prove.rs](src/isogeny_prove.rs) contains the prover and verifier for the zk-PolyIOP for 2-isogeny walks
- [src/main.rs](src/main.rs) contains the script which generated the benchmarks from the isogeny walk
- [src/matrix.rs](src/matrix.rs) contains the 'solve_linear_system' function which is used in the generalized FRI algorithm
- [src/merkle.rs](src/merkle.rs) sets the parameters for the Merkle trees used in the protocol
- [src/sidh_sike_p434.rs](src/sidh_sike_p434.rs) defines the field used in the protocol. 


The repository has been built using [Arkworks](https://github.com/arkworks-rs) for finite field algebra, polynomials, and hashing.
