# SliQSim Final Merge Report

## Source Project Capabilities

`SliQSim_multi` was used as the base project. It supports labeled multi-initial-state simulation with selector BDD variables, basis initial states, dense arbitrary-amplitude initial states, sparse arbitrary-amplitude initial states, repeated access to selected output labels, all-in-one versus sequential verification, observation queries per label, and BDD DOT dumps.

`SliQSim` contributes dynamic-circuit execution. The merged final project supports mid-circuit measurement by projecting and renormalizing the current BDD state, active qubit reset by measuring and conditionally applying `x`, classical-register conditions through `if`, and repeat-until-success style `while` loops.

`SliQSim-QDA` contributes additional gates and probability views. The merged final project supports `ct q[c],q[t]` as a controlled T gate, `majx q[c1],q[c2],q[c3],q[t]` as a majority-controlled X gate, and simulation types 3, 4, and 5 for basis probability listings, probability extremes, and ranged probability queries.

## Merged Method

The final implementation keeps the original SliQSim representation: each amplitude is represented by four fixed-point integer components, and each component bit is a BDD root. Multi-initial-state simulation adds selector variables after the quantum variables. Initial states are encoded by OR-combining label cubes with the BDDs for the corresponding basis amplitudes. Quantum gates still act on the quantum variables only, so the selector variables remain as symbolic labels that partition independent initial states.

Dynamic circuits mutate the current BDD state because measurement collapse and reset are state-changing operations. For a selected multi-initial-state label, the final project first cofactors the all-in-one BDD state by that label, snapshots the selected state, and restores that snapshot before each shot. This preserves correct shot semantics for mid-circuit measurement while still allowing one command to process all labels.

## Usage

Build:

```bash
make
```

Static multi-initial simulation with QDA gates:

```bash
./SliQSim --sim_qasm benchmarks/final_merge/multi_qda_static.qasm --init_states benchmarks/final_merge/init_4q_basis.txt --shots 1000 --seed 0
```

Arbitrary-amplitude multi-initial simulation with probability extremes:

```bash
./SliQSim --sim_qasm benchmarks/final_merge/multi_qda_static.qasm --init_states benchmarks/final_merge/init_4q_arbitrary.txt --type 4 --seed 0
```

Repeat-until-success style dynamic simulation with multi-initial states and QDA gates:

```bash
./SliQSim --sim_qasm benchmarks/final_merge/rus_multi_qda.qasm --init_states benchmarks/final_merge/init_4q_basis.txt --shots 100 --seed 0
```

All-in-one versus sequential verification:

```bash
./SliQSim --sim_qasm benchmarks/final_merge/multi_qda_static.qasm --init_states benchmarks/final_merge/init_4q_basis.txt --verify_init --verify_tol 0.05 --shots 1000 --seed 0
```

## Smoke-Test Results

The final project was compiled successfully after regenerating CUDD's local configuration for the installed compiler.

Executed checks:

- Static multi-initial sampling with `ct` and `majx`: passed.
- Sparse arbitrary-amplitude multi-initial probability extremes: passed.
- Dynamic repeat-until-success style simulation with mid-circuit measurement, reset, labels, `ct`, and `majx`: passed.
- All-in-one versus sequential probability verification: passed.
