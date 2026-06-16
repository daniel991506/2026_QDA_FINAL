#!/usr/bin/env python3
"""Compare SliQSim multi-init, SliQSim sequential-init, and MQT DDSIM.

The script is intentionally report-oriented: it records runtime, peak RSS,
SliQSim's own node/gate statistics, output support, and probability error.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import queue
import re
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class InitState:
    label: str
    entries: list[str]


def parse_complex(raw: str) -> complex:
    text = raw.strip().replace("i", "j")
    if text in {"j", "+j"}:
        text = "1j"
    elif text == "-j":
        text = "-1j"
    return complex(text)


def parse_init_states(path: Path, selection: str) -> list[InitState]:
    selected = None
    if selection and selection != "all":
        selected = {item.strip() for item in selection.split(",") if item.strip()}

    states: list[InitState] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("//", 1)[0].split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        label, entries = parts[0], parts[1:]
        if selected is None or label in selected:
            states.append(InitState(label, entries))
    return states


def is_bitstring(value: str) -> bool:
    return bool(value) and all(ch in "01" for ch in value)


def state_to_probabilities(payload: dict[str, Any]) -> dict[str, float]:
    if "counts" in payload:
        counts = {str(k): int(v) for k, v in payload["counts"].items()}
        total = sum(counts.values()) or 1
        return {state: count / total for state, count in counts.items()}

    vector = payload.get("statevector")
    if vector is None:
        return {}
    probs: dict[str, float] = {}
    nstates = len(vector)
    nqubits = int(math.log2(nstates)) if nstates > 0 else 0
    for index, amp in enumerate(vector):
        if isinstance(amp, list) and len(amp) == 2:
            value = complex(amp[0], amp[1])
        else:
            value = complex(amp)
        prob = abs(value) ** 2
        if prob > 1e-15:
            probs[format(index, f"0{nqubits}b")] = prob
    return probs


def parse_sliqsim_statevector(text: str) -> list[complex]:
    values: list[complex] = []
    for raw in text.split(","):
        token = raw.strip()
        if not token:
            continue
        values.append(parse_complex(token))
    return values


def vector_to_probabilities(vector: list[complex]) -> dict[str, float]:
    probs: dict[str, float] = {}
    nstates = len(vector)
    nqubits = int(math.log2(nstates)) if nstates > 0 else 0
    for index, value in enumerate(vector):
        prob = abs(value) ** 2
        if prob > 1e-15:
            probs[format(index, f"0{nqubits}b")] = prob
    return probs


def parse_sliqsim_non_json_line(line: str) -> tuple[str, dict[str, float]] | None:
    label_match = re.search(r'"initial_state"\s*:\s*"([^"]+)"', line)
    if not label_match:
        return None
    statevector_match = re.search(r'"statevector"\s*:\s*\[([^\]]*)\]', line)
    if statevector_match:
        vector = parse_sliqsim_statevector(statevector_match.group(1))
        return label_match.group(1), vector_to_probabilities(vector)
    counts_match = re.search(r'"counts"\s*:\s*\{([^}]*)\}', line)
    if counts_match:
        counts: dict[str, int] = {}
        for state, count in re.findall(r'"([^"]+)"\s*:\s*([0-9]+)', counts_match.group(1)):
            counts[state] = int(count)
        total = sum(counts.values()) or 1
        return label_match.group(1), {state: count / total for state, count in counts.items()}
    return None


def parse_sliqsim_output(stdout: str) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    unlabeled_index = 0
    for line in stdout.splitlines():
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError:
            parsed = parse_sliqsim_non_json_line(line)
            if parsed:
                label, probabilities = parsed
                result[label] = probabilities
            continue
        if "initial_state" in payload:
            label = str(payload["initial_state"])
            body = payload.get("result", {})
        else:
            label = f"run{unlabeled_index}"
            unlabeled_index += 1
            body = payload
        result[label] = state_to_probabilities(body)
    return result


def parse_sliqsim_stats(stdout: str) -> dict[str, float]:
    stats: dict[str, float] = {}
    patterns = {
        "reported_runtime_sec": r"Runtime:\s+([0-9.eE+-]+)\s+seconds",
        "reported_peak_memory_bytes": r"Peak memory usage:\s+([0-9.eE+-]+)\s+bytes",
        "applied_gates": r"#Applied gates:\s+([0-9.eE+-]+)",
        "max_nodes": r"Max #nodes:\s+([0-9.eE+-]+)",
        "initial_states": r"#Initial states:\s+([0-9.eE+-]+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, stdout)
        if match:
            stats[key] = float(match.group(1))
    return stats


def monitor_process_peak(proc: subprocess.Popen[str], stop: threading.Event, peak: dict[str, int]) -> None:
    try:
        import psutil
    except ImportError:
        return
    try:
        ps_proc = psutil.Process(proc.pid)
    except psutil.Error:
        return
    while not stop.is_set():
        try:
            rss = ps_proc.memory_info().rss
            for child in ps_proc.children(recursive=True):
                try:
                    rss += child.memory_info().rss
                except psutil.Error:
                    pass
            peak["bytes"] = max(peak.get("bytes", 0), rss)
        except psutil.Error:
            pass
        time.sleep(0.01)


def run_command(cmd: list[str], cwd: Path) -> dict[str, Any]:
    start = time.perf_counter()
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    stop = threading.Event()
    peak = {"bytes": 0}
    thread = threading.Thread(target=monitor_process_peak, args=(proc, stop, peak), daemon=True)
    thread.start()
    stdout, stderr = proc.communicate()
    stop.set()
    thread.join(timeout=0.2)
    wall = time.perf_counter() - start
    return {
        "cmd": cmd,
        "returncode": proc.returncode,
        "stdout": stdout,
        "stderr": stderr,
        "wall_time_sec": wall,
        "observed_peak_rss_bytes": peak["bytes"] or None,
    }


def run_sliqsim(args: argparse.Namespace, sequential: bool) -> dict[str, Any]:
    cmd = [
        str(args.sliqsim),
        "--sim_qasm",
        str(args.qasm),
        "--type",
        str(args.type),
        "--shots",
        str(args.shots),
        "--seed",
        str(args.seed),
        "--r",
        str(args.r),
        "--init_states",
        str(args.init_states),
        "--select_init",
        args.select_init,
        "--print_info",
    ]
    if sequential:
        cmd.append("--sequential_init")
    run = run_command(cmd, args.repo)
    run["probabilities"] = parse_sliqsim_output(run["stdout"])
    run["reported_stats"] = parse_sliqsim_stats(run["stdout"])
    return run


def compare_probabilities(
    left: dict[str, dict[str, float]],
    right: dict[str, dict[str, float]],
) -> dict[str, Any]:
    labels = sorted(set(left) | set(right))
    missing_labels = {
        "left": sorted(set(right) - set(left)),
        "right": sorted(set(left) - set(right)),
    }
    max_abs_diff = 0.0
    max_label = ""
    max_state = ""
    total_variation_by_label: dict[str, float] = {}
    support_mismatches: list[dict[str, str]] = []

    for label in labels:
        ldist = left.get(label, {})
        rdist = right.get(label, {})
        states = sorted(set(ldist) | set(rdist))
        tvd = 0.0
        for state in states:
            in_left = state in ldist
            in_right = state in rdist
            if in_left != in_right:
                support_mismatches.append({"label": label, "state": state})
            diff = abs(ldist.get(state, 0.0) - rdist.get(state, 0.0))
            tvd += diff
            if diff > max_abs_diff:
                max_abs_diff = diff
                max_label = label
                max_state = state
        total_variation_by_label[label] = 0.5 * tvd

    return {
        "same_label_set": not missing_labels["left"] and not missing_labels["right"],
        "support_match": not support_mismatches,
        "missing_labels": missing_labels,
        "support_mismatches": support_mismatches[:20],
        "max_abs_probability_diff": max_abs_diff,
        "max_diff_label": max_label,
        "max_diff_state": max_state,
        "total_variation_by_label": total_variation_by_label,
    }


def dense_vector_from_state(state: InitState, nqubits: int) -> tuple[list[complex] | None, str | None]:
    entries = state.entries
    if len(entries) == 1 and is_bitstring(entries[0]):
        return None, None

    size = 1 << nqubits
    vector = [0j] * size
    if any(":" in entry for entry in entries):
        for entry in entries:
            if ":" not in entry:
                return None, f"mixed sparse/dense entry '{entry}'"
            bitstring, raw_amp = entry.split(":", 1)
            if not is_bitstring(bitstring) or len(bitstring) != nqubits:
                return None, f"invalid sparse basis '{bitstring}'"
            vector[int(bitstring, 2)] = parse_complex(raw_amp)
    else:
        if len(entries) != size:
            return None, f"dense state has {len(entries)} amplitudes, expected {size}"
        vector = [parse_complex(entry) for entry in entries]
    norm = math.sqrt(sum(abs(value) ** 2 for value in vector))
    if norm == 0:
        return None, "zero-norm state"
    vector = [value / norm for value in vector]
    return vector, None


def build_ddsim_circuit(qasm_path: Path, state: InitState):
    from qiskit import QuantumCircuit, transpile

    original = QuantumCircuit.from_qasm_file(str(qasm_path))
    nqubits = original.num_qubits
    nclbits = max(original.num_clbits, 1)
    circuit = QuantumCircuit(nqubits, nclbits)

    if len(state.entries) == 1 and is_bitstring(state.entries[0]):
        bits = state.entries[0]
        if len(bits) != nqubits:
            return None, f"basis state length {len(bits)} does not match {nqubits} qubits"
        for qindex, bit in enumerate(reversed(bits)):
            if bit == "1":
                circuit.x(qindex)
    else:
        vector, error = dense_vector_from_state(state, nqubits)
        if error:
            return None, error
        circuit.initialize(vector, list(range(nqubits)))

    circuit.compose(original, qubits=list(range(nqubits)), clbits=list(range(original.num_clbits)), inplace=True)
    basis = ["id", "x", "y", "z", "h", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "u1", "u2", "u3", "cx", "cz", "ccx", "swap"]
    circuit = transpile(circuit, basis_gates=basis, optimization_level=0)
    return circuit, None


def run_ddsim_label_worker(
    qasm_path: str,
    state: InitState,
    sim_type: int,
    shots: int,
    seed: int,
    out_queue: mp.Queue,
) -> None:
    try:
        from mqt.ddsim import DDSIMProvider

        backend_name = "qasm_simulator" if sim_type == 0 else "statevector_simulator"
        backend = DDSIMProvider().get_backend(backend_name)
        circuit, error = build_ddsim_circuit(Path(qasm_path), state)
        if error:
            out_queue.put({"status": "error", "label": state.label, "reason": error})
            return

        start = time.perf_counter()
        result = backend.run(circuit, shots=shots, seed_simulator=seed).result()
        wall = time.perf_counter() - start
        if sim_type == 0:
            raw_counts = result.get_counts(circuit)
            counts = {}
            for key, value in raw_counts.items():
                state_bits = format(int(key, 16), f"0{circuit.num_clbits}b") if key.startswith("0x") else key
                counts[state_bits] = int(value)
            total = sum(counts.values()) or 1
            probabilities = {key: value / total for key, value in counts.items()}
        else:
            vector = list(result.get_statevector(circuit))
            nstates = len(vector)
            nqubits = int(math.log2(nstates)) if nstates else 0
            probabilities = {
                format(index, f"0{nqubits}b"): abs(complex(amp)) ** 2
                for index, amp in enumerate(vector)
                if abs(complex(amp)) ** 2 > 1e-15
            }

        backend_time = None
        try:
            backend_time = float(result.results[0].data.time_taken)
        except Exception:
            pass
        out_queue.put({
            "status": "ok",
            "label": state.label,
            "probabilities": probabilities,
            "stats": {
                "wall_time_sec": wall,
                "backend_time_sec": backend_time,
                "num_qubits": circuit.num_qubits,
                "num_clbits": circuit.num_clbits,
                "circuit_depth": circuit.depth(),
                "circuit_size": circuit.size(),
            },
        })
    except (MemoryError, OverflowError) as exc:
        out_queue.put({
            "status": "aborted",
            "label": state.label,
            "reason": f"DDSIM dense initial-state preparation is infeasible for this label: {exc}",
        })
    except Exception as exc:
        out_queue.put({"status": "error", "label": state.label, "reason": str(exc)})


def process_rss_bytes(pid: int) -> int:
    try:
        import psutil
    except ImportError:
        return 0
    try:
        proc = psutil.Process(pid)
        rss = proc.memory_info().rss
        for child in proc.children(recursive=True):
            try:
                rss += child.memory_info().rss
            except psutil.Error:
                pass
        return rss
    except psutil.Error:
        return 0


def run_ddsim_label_with_timeout(args: argparse.Namespace, state: InitState, backend_name: str) -> dict[str, Any]:
    ctx = mp.get_context("fork")
    out_queue = ctx.Queue()
    proc = ctx.Process(
        target=run_ddsim_label_worker,
        args=(str(args.qasm), state, args.type, args.shots, args.seed, out_queue),
    )
    start = time.perf_counter()
    proc.start()
    peak_rss = 0
    while proc.is_alive():
        if proc.pid is not None:
            peak_rss = max(peak_rss, process_rss_bytes(proc.pid))
        if time.perf_counter() - start >= args.ddsim_timeout:
            proc.terminate()
            proc.join(timeout=1.0)
            if proc.is_alive():
                proc.kill()
                proc.join(timeout=1.0)
            return {
                "status": "aborted",
                "label": state.label,
                "reason": f"hit DDSIM per-label wall-time limit of {args.ddsim_timeout} seconds",
                "wall_time_sec": time.perf_counter() - start,
                "observed_peak_rss_bytes": peak_rss,
            }
        time.sleep(0.01)

    proc.join()
    wall = time.perf_counter() - start
    if proc.pid is not None:
        peak_rss = max(peak_rss, process_rss_bytes(proc.pid))
    try:
        item = out_queue.get_nowait()
    except queue.Empty:
        item = {
            "status": "error",
            "label": state.label,
            "reason": f"DDSIM worker exited with code {proc.exitcode} without returning a result",
        }
    item["wall_time_sec"] = wall
    item["observed_peak_rss_bytes"] = peak_rss
    item["backend"] = backend_name
    return item


def monitor_self_peak(stop: threading.Event, peak: dict[str, int]) -> None:
    try:
        import psutil
    except ImportError:
        return
    proc = psutil.Process(os.getpid())
    while not stop.is_set():
        try:
            peak["bytes"] = max(peak.get("bytes", 0), proc.memory_info().rss)
        except psutil.Error:
            pass
        time.sleep(0.01)


def run_ddsim(args: argparse.Namespace, states: list[InitState]) -> dict[str, Any]:
    try:
        from mqt.ddsim import DDSIMProvider
    except ImportError as exc:
        return {"available": False, "error": str(exc), "probabilities": {}}

    backend_name = "qasm_simulator" if args.type == 0 else "statevector_simulator"
    DDSIMProvider().get_backend(backend_name)

    probabilities: dict[str, dict[str, float]] = {}
    per_label: dict[str, dict[str, Any]] = {}
    aborted: dict[str, str] = {}
    errors: dict[str, str] = {}
    start_total = time.perf_counter()
    peak = {"bytes": 0}
    stop = threading.Event()
    thread = threading.Thread(target=monitor_self_peak, args=(stop, peak), daemon=True)
    thread.start()
    try:
        for state in states:
            result = run_ddsim_label_with_timeout(args, state, backend_name)
            peak["bytes"] = max(peak["bytes"], result.get("observed_peak_rss_bytes") or 0)
            if result["status"] == "ok":
                probabilities[state.label] = result["probabilities"]
                per_label[state.label] = result["stats"]
                per_label[state.label]["total_wall_time_sec"] = result["wall_time_sec"]
                per_label[state.label]["observed_peak_rss_bytes"] = result["observed_peak_rss_bytes"]
            elif result["status"] == "aborted":
                aborted[state.label] = result["reason"]
                per_label[state.label] = {
                    "total_wall_time_sec": result["wall_time_sec"],
                    "observed_peak_rss_bytes": result["observed_peak_rss_bytes"],
                    "status": "aborted",
                }
            else:
                errors[state.label] = result["reason"]
                per_label[state.label] = {
                    "total_wall_time_sec": result["wall_time_sec"],
                    "observed_peak_rss_bytes": result["observed_peak_rss_bytes"],
                    "status": "error",
                }
    finally:
        stop.set()
        thread.join(timeout=0.2)

    return {
        "available": True,
        "backend": backend_name,
        "wall_time_sec": time.perf_counter() - start_total,
        "observed_peak_rss_bytes": peak["bytes"],
        "probabilities": probabilities,
        "per_label": per_label,
        "aborted": aborted,
        "errors": errors,
    }


def fmt(value: Any) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def print_report(report: dict[str, Any]) -> None:
    print("# Simulator Comparison Report")
    print()
    print("## Configuration")
    cfg = report["configuration"]
    for key in ["qasm", "init_states", "selected_labels", "type", "shots", "seed", "r", "ddsim_timeout_sec"]:
        print(f"- {key}: {cfg[key]}")
    print()

    print("## Runtime and Size Statistics")
    print("| mode | wall time (s) | observed peak RSS (bytes) | reported runtime (s) | reported peak memory (bytes) | applied gates | max BDD nodes |")
    print("|---|---:|---:|---:|---:|---:|---:|")
    for key, title in [
        ("sliqsim_multi", "SliQSim one-run multi-init"),
        ("sliqsim_sequential", "SliQSim sequential"),
        ("ddsim", "MQT DDSIM sequential"),
    ]:
        item = report[key]
        stats = item.get("reported_stats", {})
        print(
            f"| {title} | {fmt(item.get('wall_time_sec'))} | {fmt(item.get('observed_peak_rss_bytes'))} | "
            f"{fmt(stats.get('reported_runtime_sec'))} | {fmt(stats.get('reported_peak_memory_bytes'))} | "
            f"{fmt(stats.get('applied_gates'))} | {fmt(stats.get('max_nodes'))} |"
        )
    print()

    print("## Correctness / Probability Comparison")
    print("| comparison | same labels | support match | max abs prob diff | max label | max state |")
    print("|---|---:|---:|---:|---|---|")
    for key, title in [
        ("multi_vs_sequential", "SliQSim multi vs SliQSim sequential"),
        ("multi_vs_ddsim", "SliQSim multi vs DDSIM"),
        ("sequential_vs_ddsim", "SliQSim sequential vs DDSIM"),
    ]:
        comp = report["comparisons"].get(key, {})
        print(
            f"| {title} | {fmt(comp.get('same_label_set'))} | {fmt(comp.get('support_match'))} | "
            f"{fmt(comp.get('max_abs_probability_diff'))} | {fmt(comp.get('max_diff_label'))} | {fmt(comp.get('max_diff_state'))} |"
        )
    print()

    ddsim = report["ddsim"]
    if ddsim.get("per_label"):
        print("## DDSIM Per-Initial-State Statistics")
        print("| label | status | total wall time (s) | backend time (s) | peak RSS (bytes) | qubits | depth | operations |")
        print("|---|---|---:|---:|---:|---:|---:|---:|")
        for label, item in ddsim["per_label"].items():
            print(
                f"| {label} | {fmt(item.get('status', 'ok'))} | {fmt(item.get('total_wall_time_sec'))} | "
                f"{fmt(item.get('backend_time_sec'))} | {fmt(item.get('observed_peak_rss_bytes'))} | "
                f"{fmt(item.get('num_qubits'))} | {fmt(item.get('circuit_depth'))} | {fmt(item.get('circuit_size'))} |"
            )
        print()
    if ddsim.get("aborted"):
        print("## DDSIM Aborted Labels")
        for label, reason in ddsim["aborted"].items():
            print(f"- {label}: {reason}")
        print()
    if ddsim.get("errors"):
        print("## DDSIM Error Labels")
        for label, reason in ddsim["errors"].items():
            print(f"- {label}: {reason}")
        print()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--sliqsim", type=Path, default=Path("./SliQSim"))
    parser.add_argument("--qasm", type=Path, required=True)
    parser.add_argument("--init-states", type=Path, required=True)
    parser.add_argument("--select-init", default="all")
    parser.add_argument("--type", type=int, choices=[0, 1], default=1)
    parser.add_argument("--shots", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--r", type=int, default=32)
    parser.add_argument("--ddsim-timeout", type=float, default=30.0, help="per-label DDSIM wall-time limit in seconds")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    args.repo = args.repo.resolve()
    args.qasm = args.qasm.resolve()
    args.init_states = args.init_states.resolve()
    if not args.sliqsim.is_absolute():
        args.sliqsim = (args.repo / args.sliqsim).resolve()

    states = parse_init_states(args.init_states, args.select_init)
    selected_labels = [state.label for state in states]

    multi = run_sliqsim(args, sequential=False)
    sequential = run_sliqsim(args, sequential=True)
    ddsim = run_ddsim(args, states)

    comparisons = {
        "multi_vs_sequential": compare_probabilities(multi["probabilities"], sequential["probabilities"]),
        "multi_vs_ddsim": compare_probabilities(multi["probabilities"], ddsim.get("probabilities", {})),
        "sequential_vs_ddsim": compare_probabilities(sequential["probabilities"], ddsim.get("probabilities", {})),
    }

    report = {
        "configuration": {
            "qasm": str(args.qasm),
            "init_states": str(args.init_states),
            "selected_labels": selected_labels,
            "type": args.type,
            "shots": args.shots,
            "seed": args.seed,
            "r": args.r,
            "ddsim_timeout_sec": args.ddsim_timeout,
        },
        "sliqsim_multi": multi,
        "sliqsim_sequential": sequential,
        "ddsim": ddsim,
        "comparisons": comparisons,
    }

    print_report(report)
    if args.json_out:
        args.json_out.write_text(json.dumps(report, indent=2, default=str) + "\n")
        print(f"Raw JSON written to {args.json_out}")

    failed = multi["returncode"] != 0 or sequential["returncode"] != 0
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
