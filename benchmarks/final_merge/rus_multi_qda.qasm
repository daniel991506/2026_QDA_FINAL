OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[2];
measure q[0] -> c[0];
while (c==0) {
reset q[0];
x q[0];
measure q[0] -> c[0];
}
ct q[0],q[1];
majx q[3],q[0],q[1],q[2];
measure q[1] -> c[1];
