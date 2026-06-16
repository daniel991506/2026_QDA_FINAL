# SliQSim - A BDD-based Quantum Circuit Simulator

## Introduction
`SliQSim` is a BDD-based quantum circuit simulator implemented in C/C++ on top of [CUDD](http://web.mit.edu/sage/export/tmp/y/usr/share/doc/polybori/cudd/cuddIntro.html) package. In `SliQSim`, a bit-slicing technique based on BDDs is used to represent quantum state vectors. For more details of the simulator, please refer to the [paper](https://arxiv.org/abs/2007.09304).

## Build
To build the simulator, one needs to first configure `CUDD`:
```commandline
cd cudd
./configure --enable-dddmp --enable-obj --enable-shared --enable-static 
cd ..
```
Next, build the binary file `SliQSim`:
```commandline
make
```

## Execution
The circuit format being simulated is `OpenQASM` used by IBM's [Qiskit](https://github.com/Qiskit/qiskit), and the gate set supported in this simulator now contains Pauli-X (x), Pauli-Y (y), Pauli-Z (z), Hadamard (h), Phase and its inverse (s and sdg), π/8 and its inverse (t and tdg), Rotation-X with phase π/2 (rx(pi/2)), Rotation-Y with phase π/2 (ry(pi/2)), Controlled-NOT (cx), Controlled-Z (cz), Toffoli (ccx and mcx), SWAP (swap), and Fredkin (cswap). One can find some example benchmarks in `examples` folder. 

For simulation types, we provide "sampling", "all_amplitude", and "query" simulation options. The help message states the details:

```commandline
$ ./SliQSim --help
Options:
--help                produce help message
--sim_qasm arg        simulate qasm file string
--seed [=arg(=1)]     seed for random number generator
--print_info          print simulation statistics such as runtime, memory, etc.
--type arg (=0)       the simulation type being executed.
                      0: sampling mode (default option), where the sampled outcomes will be provided. 
                      1: all_amplitude mode, where the final state vector will be shown. 
                      2: query mode, where only the values of properties defined in obs_file will be provided.
--shots arg (=1)      the number of outcomes being sampled in "sampling mode" .
--obs_file arg        self-defined measurement operation file string (if any).
--init_states arg     multi-initial-state file. Each non-comment line is:
                      <label> <bitstring>, <label> <amp0> <amp1> ...,
                      or sparse entries <label> <bitstring>:<amp> ...
--select_init arg (=all)
                      initial-state label(s) to access after simulation. Use
                      'all' or comma-separated labels.
--sequential_init     run selected initial states one by one instead of using
                      selector variables. Useful for verification.
--verify_init         run both all-in-one and sequential initial-state modes,
                      then compare their output.
--verify_tol arg (=1e-09)
                      probability tolerance for --verify_init counts comparison.
--dump_bdds arg       write BDD Graphviz .dot dumps using the given filename
                      prefix.
--r arg (=32)         integer bit size.
--reorder arg (=1)    allow variable reordering or not.
                      0: disable reordering.
                      1: enable reordering (default option).
--alloc arg (=1)      allocate new BDDs when overflow is detected.
                      0: do not allocate new BDDs. This may lead to numerical errors.
                      1: allocate new BDDs (default option).
./SliQSim --sim_qasm examples/bell_state.qasm --type 1 \
  --init_states init_states.txt \
  --dump_bdds /tmp/sliqsim_bdd
```
To use the sampling mode (default), it is required to have measurement operations included in the qasm file. Conversely, in all_amplitude mode, measurement operations are generally omitted, but if they are present in the qasm file, the final state vector will collapse based on the measurement result. It is important to note that all_amplitude mode is not recommended for simulations involving a large number of qubits, as it could result in a significantly long runtime.

## Multi-Initial-State Simulation

SliQSim can simulate multiple labeled initial states together through the same circuit. Each non-comment line in the initial-state file can be either a computational-basis bitstring or a full statevector written in computational-basis order:

```text
zero 00
one 01
plus 0.70710678118 0.70710678118
i_state 0 1i
sparse_plus 000:0.70710678118 111:0.70710678118
```

Then run with `--init_states`. By default, SliQSim accesses every labeled output state after the shared circuit simulation:

```commandline
./SliQSim --sim_qasm examples/bell_state.qasm --type 1 --init_states init_states.txt
```

To access only some output states, pass `--select_init` with one label or a comma-separated list:

```commandline
./SliQSim --sim_qasm examples/bell_state.qasm --type 1 --init_states init_states.txt --select_init zero
```

The same access mechanism works with sampling, all-amplitude statevector extraction, and query mode. Internally, the labeled initial states are encoded with auxiliary Boolean selector variables, evolved in one BDD simulation, and projected by label when an output state is requested.

Dense amplitude-vector entries use `2^n` amplitudes for an `n`-qubit circuit. Sparse entries use only the nonzero terms, written as `<bitstring>:<amplitude>`, and are recommended for large circuits such as `examples/bv100.qasm`. Real and complex values such as `0.5`, `-0.5`, `1i`, `-1i`, and `0.5+0.5i` are supported. Decimal amplitudes are stored in SliQSim's bit-sliced integer representation using fixed-point scaling, so values are approximate at the precision controlled by `--r`.

To verify a multi-initial-state run, add `--sequential_init`. This runs the same selected labels one by one as independent simulations instead of using selector variables:

```commandline
./SliQSim --sim_qasm examples/bv100.qasm --shots 1 --seed 0 --init_states examples/bv100_sparse_arbitrary_init.txt --select_init arb00,arb02 --sequential_init
```

The output format is the same as the combined multi-initial-state mode, so selected labels can be compared directly.

To have SliQSim perform that comparison automatically, use `--verify_init`:

```commandline
./SliQSim --sim_qasm examples/bv100.qasm --shots 1 --seed 0 --init_states examples/bv100_sparse_arbitrary_init.txt --select_init arb00,arb02 --verify_init
```

For sampling output, `--verify_init` first checks whether both modes produced the same set of selected labels and output bitstrings. It then converts counts to probabilities and compares every output bitstring for every selected label. The default tolerance is `1e-9`; override it with `--verify_tol` when comparing finite-shot samples:

```commandline
./SliQSim --sim_qasm examples/bv100.qasm --shots 1000 --seed 0 --init_states examples/bv100_sparse_arbitrary_init.txt --select_init arb00,arb02 --verify_init --verify_tol 0.05
```

This prints `{ "verification": "pass" }` when the support matches and the probability differences are within tolerance. If a label or output bitstring appears in one mode but not the other, the report includes `support_match: false`, `missing_from`, `missing_label`, and, when applicable, `missing_output`. If the outputs are not sampling counts, it falls back to exact output comparison.

To inspect the BDDs, use `--dump_bdds` with a filename prefix:

```commandline
./SliQSim --sim_qasm examples/bell_state.qasm --type 1 --init_states init_states.txt --dump_bdds /tmp/sliqsim_bdd
```

This writes Graphviz `.dot` files such as:

```text
/tmp/sliqsim_bdd_initial_multi_input_slices.dot
/tmp/sliqsim_bdd_final_multi_state_slices.dot
/tmp/sliqsim_bdd_final_merged_all_in_one_bigBDD.dot
```

For a normal single-initial-state run, the corresponding files are:

```text
/tmp/sliqsim_bdd_initial_input_slices.dot
/tmp/sliqsim_bdd_final_merged_bigBDD.dot
```

The slice files contain the `All_Bdd[component][bit]` roots (`c0_b0`, `c0_b1`, ...), and the merged file contains the `bigBDD` root used for measurement/probability traversal. For large circuits these files can be very large.

To generate report-ready comparison statistics for SliQSim's all-in-one multi-initial-state mode, SliQSim's sequential mode, and MQT DDSIM's sequential simulation, use:

```commandline
python3 tools/compare_simulators.py --qasm examples/bell_state.qasm --init-states init_state/init_states.txt --type 1 --seed 0 --ddsim-timeout 30 --json-out /tmp/sliqsim_compare.json
```

The script prints Markdown tables with runtime, observed peak memory, SliQSim's reported runtime, SliQSim's max BDD nodes, DDSIM per-label runtime/depth/operation counts, and probability-agreement statistics. MQT DDSIM is run label by label because it does not use SliQSim's selector-variable multi-initial-state encoding. Basis initial states are supported directly. Arbitrary amplitude DDSIM initial states are attempted for every selected label by expanding the state vector and prepending an initialization circuit. To keep very large cases bounded, `--ddsim-timeout` sets a hard per-label DDSIM wall-time limit; labels that hit the limit or require an infeasible dense state vector are reported as aborted.

For example, simulating `examples/bell_state_measure.qasm`, which is a 2-qubit bell state circuit with measurement gates at the end, with the sampling mode simulation option can be executed by
```commandline
./SliQSim --sim_qasm examples/bell_state_measure.qasm --type 0 --shots 1024
```

Then the sampled results will be shown:
```commandline
{ "counts": { "11": 542, "00": 482 } }
```

If option `--print_info` is used, simulation statistics such as runtime and memory usage will also be provided: 
```commandline
  Runtime: 0.014433 seconds
  Peak memory usage: 12611584 bytes
  #Applied gates: 2
  Max #nodes: 13
  Precision of integers: 32
  Accuracy loss: 2.22045e-16
```

To demonstrate the all_amplitude mode simulation, we use `examples/bell_state.qasm`, which is the same circuit as in the sampling mode simulation example except that the measurement gates are removed:
```commandline
./SliQSim --sim_qasm examples/bell_state.qasm --type 1
```

This will show the resulting state vector:
```commandline
{"statevector": ["0.707107", "0", "0", "0.707107"] }
```

`SliQSim` also supports various types of queries for user-defined properties. The users can use the expressions listed below to specify their desired properties.

* `bf {formula}`: It returns the probability that the measurement result satisfies the specified formula.
* `hweq/hwneq/hwgt/hwlt {qubits} {num}`: It returns the probability that the Hamming weight of specified qubits is equal to, non-equal to, greater than, lower than the specified number num.
* `inteq/intneq/intgt/intlt {qubits} {num}`: It returns the probability that the binary integer represented by the specified qubits is equal to, non-equal to, greater than, lower than the specified number num.
* `expt {qubits} {Pauli(-value)_string}`: It returns the expectation value of the observable defined by the specified Pauli string with optionally specified measured-values (0 or 1) for each individual Pauli operator over the specified qubits. Specifically, the returned value is calculated by the product of the probability of obtaining the specified measured values and the expectation value of the Pauli operators without measured values with respect to the post-measurement state.
* `weightedsum {weight} {expression} {...} endweighteedsum`: It returns the weighted sum of values of the specified expressions. Each weight and its corresponding expression are stated in an independent line, and the lines are clipped by keywords weightedsum and endweightedsum. The expression can be `bf`, `hweq`, `hwneq`, `hwgt`, `hwlt`, `inteq`, `intneq`, `intgt`, `intlt`, or `expt`.
* `assign {var_name} {expression}`: It does not return values, but rather stores the specified expression in the specified variable name var_name, which can be further utilized in `bf` expression. The expression can be `bf`, `hweq`, `hwneq`, `hwgt`, `hwlt`, `inteq`, `intneq`, `intgt`, or `intlt`.
* `between/outof/leq/geq {threshold}`: It can be specified before the expressions: `bf`, `hweq`, `hwneq`, `hwgt`, `hwlt`, `inteq`, `intneq`, `intgt`, `intlt`, `expt`, `weightedsum`. The predicate returns
true or false according to whether or not the probability returned by its subsequent expression function is between, out of, less than or equal to, or greater than or equal to the specified range.
* `amp {compt_basis}`: It returns the probability amplitude (as a complex number) of the specified computational basis compt_basis.
* `dist {qubits}`: It returns the exact spectrum of the probability distribution upon measuring the specified qubits.

To demonstrate the query mode simulation, we use `examples/grover_10.qasm` and `examples/demo.obs`:
```commandline
./SliQSim --sim_qasm examples/grover_10.qasm --obs_file examples/demo.obs --type 2
```

This will show the query results of each user-defined properties:
```commandline
"bf tmp_var^q[2] | tmp_var":
        0.541386
"hweq q[0] q[2] 1":
        0.5
"hwgt q[0] q[2] q[1] 2":
        0.0206929
"inteq q[0] q[1] q[2] 4":
        0.0206929
"intneq q[0] q[1] q[2] 5":
        0.979307
"expt q[0] q[1] q[2] zzz":
        0
"expt q[0] q[1] q[2] yxy":
        -5.55112e-17
"weightedsum
0.750000 expt q[0] q[1] q[2] zzz
0.150000 expt q[0] q[1] q[2] izi
-0.500000 expt q[0] q[1] q[2] zzi
endweightedsum":
        -0.125169
"0.5 <= hwgt q[0] q[2] q[1] 2 <= 0.8":
        false
"amp 1001001010":
        -0.012715
"dist q[0] q[1] q[2]":
        000: 0, 001: 0.000275976, 010: 0.437921, 011: 0.0206929, 100: 0, 101: 0.000275976, 110: 0.437921, 111: 0.0206929
```

One may also execute our simulator as a backend option of Qiskit through [SliQSim Qiskit Interface](https://github.com/NTU-ALComLab/SliQSim-Qiskit-Interface), which supports "sampling" and "all_amplitude" simulation options now.


## Citation
Please cite the following paper if you use our simulator for your research:

<summary>
  <a href="https://ieeexplore.ieee.org/document/9586191">Y.-H. Tsai, J.-H. R. Jiang, and C.-S. Jhang, “Bit-slicing the Hilbert space:  Scaling up accurate quantum circuit simulation,” in <em>Design Automation Conference (DAC)</em>, 2021, pp. 439–444.</a>
</summary>

```bibtex
@INPROCEEDINGS{9586191,
  author={Tsai, Yuan-Hung and Jiang, Jie-Hong R. and Jhang, Chiao-Shan},
  booktitle={Design Automation Conference (DAC)}, 
  title={Bit-Slicing the Hilbert Space: Scaling Up Accurate Quantum Circuit Simulation}, 
  year={2021},
  pages={439-444},
  doi={10.1109/DAC18074.2021.9586191}
}
```

<summary>
  <a href="https://link.springer.com/chapter/10.1007/978-3-031-90660-2_7">T.-F. Chen, and J.-H. R. Jiang, “SliQSim: A Quantum Circuit Simulator and Solver for Probability and Statistics Queries,” in <em>International Conference on Tools and Algorithms for the Construction and Analysis of Systems (TACAS)</em>, 2025, pp. 129–138</a>
</summary>

```bibtex
@INPROCEEDINGS{chen2025sliqsim,
  author={Chen, Tian-Fu and Jiang, Jie-Hong R.},
  booktitle={International Conference on Tools and Algorithms for the Construction and Analysis of Systems (TACAS)}, 
  title={SliQSim: A Quantum Circuit Simulator and Solver for Probability and Statistics Queries}, 
  year={2025},
  pages={129-138},
  doi={10.1007/978-3-031-90660-2_7}
}
```

## Contact

If you have any questions or suggestions, feel free to [create an issue](https://github.com/NTU-ALComLab/SliQSim/issues), or contact us through matthewyhtsai@gmail.com.
