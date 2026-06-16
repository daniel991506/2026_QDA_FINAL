# Final Merge Benchmarks

These examples exercise the merged feature set:

- `multi_qda_static.qasm`: standard circuit simulation with multi-initial states plus the QDA `ct` and `majx` gates.
- `rus_multi_qda.qasm`: repeat-until-success style dynamic simulation with mid-circuit measurement, reset, a `while` loop, and QDA gates.
- `init_4q_basis.txt`: four labeled basis initial states.
- `init_4q_arbitrary.txt`: three sparse arbitrary-amplitude initial states.

Example commands from the project root:

```bash
./SliQSim --sim_qasm benchmarks/final_merge/multi_qda_static.qasm --init_states benchmarks/final_merge/init_4q_basis.txt --shots 1000 --seed 0
./SliQSim --sim_qasm benchmarks/final_merge/multi_qda_static.qasm --init_states benchmarks/final_merge/init_4q_arbitrary.txt --type 4 --seed 0
./SliQSim --sim_qasm benchmarks/final_merge/rus_multi_qda.qasm --init_states benchmarks/final_merge/init_4q_basis.txt --shots 100 --seed 0
./SliQSim --sim_qasm benchmarks/final_merge/multi_qda_static.qasm --init_states benchmarks/final_merge/init_4q_basis.txt --verify_init --verify_tol 0.05 --shots 1000 --seed 0
```
