#!/usr/bin/env python3
"""
CAMI Low Complexity Benchmark - Exhaustive Test Suite
Tests all available flags and combinations systematically.
"""

import subprocess
import os
import time
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

# ============================
# CONFIG
# ============================
BIT_POP = r"C:\Users\Daddy\Documents\GitHub\bit-pop\target\release\bit-pop.exe"
INDEX = r"D:\CAMI_low\cami_index.bitpop"
INDEX_HF = r"D:\CAMI_low\cami_index_hf.bitpop"
READS = r"D:\CAMI_low\cami_reads_100k_random.fastq"
GT = r"D:\CAMI_low\gs_read_mapping.binning\gs_read_mapping.tsv"
GENOMES = r"D:\CAMI_low\source_genomes_low\source_genomes"
OUTDIR = r"D:\CAMI_low"
RESULTS_FILE = os.path.join(OUTDIR, "benchmark_results.json")

MAX_WORKERS = 1  # Jedan po jedan
TOTAL_TIMEOUT_HOURS = 5

# ============================
# TEST GENERATOR
# ============================
def gen_tests():
    """Generate comprehensive test matrix."""
    tests = []
    test_id = 0

    # === 1. BASELINE: align mode x top-n ===
    for align in ["xor", "hybrid", "sw"]:
        for tn in [1, 2, 3]:
            test_id += 1
            tests.append({
                "id": test_id,
                "name": f"base_{align}_tn{tn}",
                "desc": f"{align} top-n={tn}",
                "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                        "-o", f"{OUTDIR}/t{test_id:03d}_{align}_tn{tn}.sam",
                        "-a", align, "-t", "16", "--top-n", str(tn)],
                "priority": "high"
            })

    # === 2. SPACED SEEDS ===
    for tn in [1, 2]:
        for pattern in [None, "11111011111111", "11101011101111"]:
            test_id += 1
            name = f"ss_tn{tn}"
            if pattern:
                name += f"_p{len(pattern)}"
                cmd = [BIT_POP, "map", "-i", INDEX, "-r", READS,
                       "-o", f"{OUTDIR}/t{test_id:03d}_ss_{name}.sam",
                       "-a", "xor", "-t", "16", "--top-n", str(tn),
                       "--spaced-seed", "--spaced-seed-pattern", pattern]
            else:
                cmd = [BIT_POP, "map", "-i", INDEX, "-r", READS,
                       "-o", f"{OUTDIR}/t{test_id:03d}_ss_{name}.sam",
                       "-a", "xor", "-t", "16", "--top-n", str(tn),
                       "--spaced-seed"]
            tests.append({
                "id": test_id,
                "name": name,
                "desc": f"spaced-seed top-n={tn}" + (f" pat={len(pattern)}" if pattern else ""),
                "cmd": cmd,
                "priority": "high"
            })

    # === 3. GOLDEN ANCHORS ===
    for tn in [1, 2]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"golden_tn{tn}",
            "desc": f"golden-anchors top-n={tn}",
            "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_golden_tn{tn}.sam",
                    "-a", "xor", "-t", "16", "--top-n", str(tn),
                    "--golden-anchors"],
            "priority": "medium"
        })

    # === 4. SEARCH RADIUS ===
    for radius in [5, 10, 20, 50, 100]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"sr{radius}_tn2",
            "desc": f"search-radius={radius} top-n=2",
            "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_sr{radius}.sam",
                    "-a", "xor", "-t", "16", "--top-n", "2",
                    "--search-radius", str(radius)],
            "priority": "medium"
        })

    # === 5. MIN SCORE ===
    for minscore in ["0.5", "0.6", "0.7", "0.8", "0.9"]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"ms{minscore}_tn2",
            "desc": f"min-score={minscore} top-n=2",
            "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_ms{minscore}.sam",
                    "-a", "xor", "-t", "16", "--top-n", "2",
                    "-m", minscore],
            "priority": "medium"
        })

    # === 6. MIN QUALITY ===
    for minqual in [0, 10, 20, 30]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"mq{minqual}_tn2",
            "desc": f"min-quality={minqual} top-n=2",
            "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_mq{minqual}.sam",
                    "-a", "xor", "-t", "16", "--top-n", "2",
                    "-q", str(minqual)],
            "priority": "low"
        })

    # === 7. CHUNK STRATEGY ===
    for strategy in ["rarest", "golden", "spaced"]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"cs_{strategy}_tn2",
            "desc": f"chunk-strategy={strategy} top-n=2",
            "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_cs_{strategy}.sam",
                    "-a", "xor", "-t", "16", "--top-n", "2",
                    "--chunk-strategy", strategy],
            "priority": "medium"
        })

    # === 8. HF (needs HF index) ===
    for tn in [1, 2]:
        for hfmin in [3, 4, 5]:
            test_id += 1
            tests.append({
                "id": test_id,
                "name": f"hf{hfmin}_tn{tn}",
                "desc": f"hf min={hfmin} top-n={tn}",
                "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                        "-o", f"{OUTDIR}/t{test_id:03d}_hf{hfmin}_tn{tn}.sam",
                        "-a", "xor", "-t", "16", "--top-n", str(tn),
                        "--hf", "--hf-min", str(hfmin)],
                "priority": "high"
            })

    # === 9. HF + SNAPSHOT ===
    for tn in [1, 2]:
        for snp_support in [3, 5]:
            test_id += 1
            tests.append({
                "id": test_id,
                "name": f"hf_snp{snp_support}_tn{tn}",
                "desc": f"hf+snp(support={snp_support}) top-n={tn}",
                "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                        "-o", f"{OUTDIR}/t{test_id:03d}_hf_snp{snp_support}_tn{tn}.sam",
                        "-a", "xor", "-t", "16", "--top-n", str(tn),
                        "--hf", "--hf-min", "3",
                        "--snp-detect", "--snp-min-support", str(snp_support)],
                "priority": "high"
            })

    # === 10. HF + GOLDEN ===
    for tn in [1, 2]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"hf_golden_tn{tn}",
            "desc": f"hf+golden top-n={tn}",
            "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_hf_golden_tn{tn}.sam",
                    "-a", "xor", "-t", "16", "--top-n", str(tn),
                    "--hf", "--hf-min", "3", "--golden-anchors"],
            "priority": "medium"
        })

    # === 11. HF + SPACED SEED ===
    for tn in [1, 2]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"hf_ss_tn{tn}",
            "desc": f"hf+spaced-seed top-n={tn}",
            "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_hf_ss_tn{tn}.sam",
                    "-a", "xor", "-t", "16", "--top-n", str(tn),
                    "--hf", "--hf-min", "3", "--spaced-seed"],
            "priority": "medium"
        })

    # === 12. HF + SEARCH RADIUS ===
    for radius in [5, 10, 20, 50]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"hf_sr{radius}_tn2",
            "desc": f"hf+search-radius={radius} top-n=2",
            "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_hf_sr{radius}.sam",
                    "-a", "xor", "-t", "16", "--top-n", "2",
                    "--hf", "--hf-min", "3", "--search-radius", str(radius)],
            "priority": "medium"
        })

    # === 13. SNAPSHOT ALONE ===
    for snp_support in [3, 5, 10]:
        for tn in [1, 2]:
            test_id += 1
            tests.append({
                "id": test_id,
                "name": f"snp{snp_support}_tn{tn}",
                "desc": f"snp(support={snp_support}) top-n={tn}",
                "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                        "-o", f"{OUTDIR}/t{test_id:03d}_snp{snp_support}_tn{tn}.sam",
                        "-a", "xor", "-t", "16", "--top-n", str(tn),
                        "--snp-detect", "--snp-min-support", str(snp_support)],
                "priority": "medium"
            })

    # === 14. SNAPSHOT + GOLDEN ===
    for tn in [1, 2]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"snp_golden_tn{tn}",
            "desc": f"snp+golden top-n={tn}",
            "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_snp_golden_tn{tn}.sam",
                    "-a", "xor", "-t", "16", "--top-n", str(tn),
                    "--snp-detect", "--snp-min-support", "3",
                    "--golden-anchors"],
            "priority": "low"
        })

    # === 15. SNAPSHOT + SPACED SEED ===
    for tn in [1, 2]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"snp_ss_tn{tn}",
            "desc": f"snp+spaced-seed top-n={tn}",
            "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_snp_ss_tn{tn}.sam",
                    "-a", "xor", "-t", "16", "--top-n", str(tn),
                    "--snp-detect", "--snp-min-support", "3",
                    "--spaced-seed"],
            "priority": "low"
        })

    # === 16. TRIPLE COMBO: HF + SNAPSHOT + GOLDEN ===
    for tn in [1, 2]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"hf_snp_golden_tn{tn}",
            "desc": f"hf+snp+golden top-n={tn}",
            "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_hf_snp_golden_tn{tn}.sam",
                    "-a", "xor", "-t", "16", "--top-n", str(tn),
                    "--hf", "--hf-min", "3",
                    "--snp-detect", "--snp-min-support", "3",
                    "--golden-anchors"],
            "priority": "high"
        })

    # === 17. TRIPLE COMBO: HF + SNAPSHOT + SPACED SEED ===
    for tn in [1, 2]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"hf_snp_ss_tn{tn}",
            "desc": f"hf+snp+spaced-seed top-n={tn}",
            "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_hf_snp_ss_tn{tn}.sam",
                    "-a", "xor", "-t", "16", "--top-n", str(tn),
                    "--hf", "--hf-min", "3",
                    "--snp-detect", "--snp-min-support", "3",
                    "--spaced-seed"],
            "priority": "high"
        })

    # === 18. QUAD COMBO: HF + SNAPSHOT + GOLDEN + SPACED SEED ===
    for tn in [1, 2]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"hf_snp_golden_ss_tn{tn}",
            "desc": f"hf+snp+golden+ss top-n={tn}",
            "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_hf_snp_golden_ss_tn{tn}.sam",
                    "-a", "xor", "-t", "16", "--top-n", str(tn),
                    "--hf", "--hf-min", "3",
                    "--snp-detect", "--snp-min-support", "3",
                    "--golden-anchors", "--spaced-seed"],
            "priority": "high"
        })

    # === 19. ALIGN MODE + HF ===
    for align in ["hybrid", "sw"]:
        for tn in [1, 2]:
            test_id += 1
            tests.append({
                "id": test_id,
                "name": f"{align}_hf_tn{tn}",
                "desc": f"{align}+hf top-n={tn}",
                "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                        "-o", f"{OUTDIR}/t{test_id:03d}_{align}_hf_tn{tn}.sam",
                        "-a", align, "-t", "16", "--top-n", str(tn),
                        "--hf", "--hf-min", "3"],
                "priority": "medium"
            })

    # === 20. ALIGN MODE + HF + SNAPSHOT ===
    for align in ["hybrid", "sw"]:
        for tn in [1, 2]:
            test_id += 1
            tests.append({
                "id": test_id,
                "name": f"{align}_hf_snp_tn{tn}",
                "desc": f"{align}+hf+snp top-n={tn}",
                "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                        "-o", f"{OUTDIR}/t{test_id:03d}_{align}_hf_snp_tn{tn}.sam",
                        "-a", align, "-t", "16", "--top-n", str(tn),
                        "--hf", "--hf-min", "3",
                        "--snp-detect", "--snp-min-support", "3"],
                "priority": "medium"
            })

    # === 21. CHUNK-CONSENSUS ===
    for pcts in ["0.01,0.10,0.50", "0.01,0.05,0.10", "0.05,0.10,0.20", "0.01,0.02,0.05"]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"cc_{pcts.replace('.', '')}",
            "desc": f"chunk-consensus {pcts}",
            "cmd": [BIT_POP, "chunk-consensus", "-i", INDEX,
                    "-c", pcts, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_cc_{pcts.replace('.', '')}.sam",
                    "-t", "16"],
            "priority": "medium"
        })

    # === 22. CHUNK-CONSENSUS + STRATEGY ===
    for strategy in ["majority", "weighted_score"]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"cc_{strategy}",
            "desc": f"chunk-consensus strategy={strategy}",
            "cmd": [BIT_POP, "chunk-consensus", "-i", INDEX,
                    "-c", "0.01,0.10,0.50", "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_cc_{strategy}.sam",
                    "-t", "16", "--strategy", strategy],
            "priority": "low"
        })

    # === 23. CHUNK-CONSENSUS + MIN AGREEMENT ===
    for minagree in [2, 3]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"cc_ma{minagree}",
            "desc": f"chunk-consensus min-agreement={minagree}",
            "cmd": [BIT_POP, "chunk-consensus", "-i", INDEX,
                    "-c", "0.01,0.10,0.50", "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_cc_ma{minagree}.sam",
                    "-t", "16", "--min-agreement", str(minagree)],
            "priority": "low"
        })

    # === 24. AUTO-K ===
    test_id += 1
    tests.append({
        "id": test_id,
        "name": "autok_tn2",
        "desc": "auto-k top-n=2",
        "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                "-o", f"{OUTDIR}/t{test_id:03d}_autok.sam",
                "-a", "xor", "-t", "16", "--top-n", "2", "--auto-k"],
        "priority": "low"
    })

    # === 25. READ TYPE ===
    for rtype in ["short", "long"]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"rt_{rtype}_tn2",
            "desc": f"read-type={rtype} top-n=2",
            "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_rt_{rtype}.sam",
                    "-a", "xor", "-t", "16", "--top-n", "2",
                    "--read-type", rtype],
            "priority": "low"
        })

    # === 26. CHAIN MODE ===
    test_id += 1
    tests.append({
        "id": test_id,
        "name": "chain_tn2",
        "desc": "chain mode top-n=2",
        "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                "-o", f"{OUTDIR}/t{test_id:03d}_chain.sam",
                "-a", "chain", "-t", "16", "--top-n", "2"],
        "priority": "low"
    })

    # === 27. SOFTCLIP MODE ===
    test_id += 1
    tests.append({
        "id": test_id,
        "name": "softclip_tn2",
        "desc": "softclip mode top-n=2",
        "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                "-o", f"{OUTDIR}/t{test_id:03d}_softclip.sam",
                "-a", "softclip", "-t", "16", "--top-n", "2"],
        "priority": "low"
    })

    # === 28. RECONCILE TOP-N (paired-end not available, skip) ===
    # === 29. TOP-N VARIANTS ===
    for tn in [4, 5]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"base_xor_tn{tn}",
            "desc": f"xor top-n={tn}",
            "cmd": [BIT_POP, "map", "-i", INDEX, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_xor_tn{tn}.sam",
                    "-a", "xor", "-t", "16", "--top-n", str(tn)],
            "priority": "low"
        })

    # === 30. HF + TOP-N VARIANTS ===
    for tn in [4, 5]:
        test_id += 1
        tests.append({
            "id": test_id,
            "name": f"hf_tn{tn}",
            "desc": f"hf top-n={tn}",
            "cmd": [BIT_POP, "map", "-i", INDEX_HF, "-r", READS,
                    "-o", f"{OUTDIR}/t{test_id:03d}_hf_tn{tn}.sam",
                    "-a", "xor", "-t", "16", "--top-n", str(tn),
                    "--hf", "--hf-min", "3"],
            "priority": "low"
        })

    return tests


# ============================
# ACCURACY CHECKER
# ============================
def run_accuracy(sam_file):
    """Run accuracy script and parse results."""
    try:
        result = subprocess.run(
            ["python", "C:\\Users\\Daddy\\Documents\\GitHub\\bit-pop\\scripts\\cami_accuracy.py",
             sam_file, GT, GENOMES],
            capture_output=True, text=True, timeout=300
        )
        output = result.stdout

        # Parse accuracy
        import re
        acc_match = re.search(r'Correct predictions:\s+\d+\s+\((\d+\.?\d*)%\s*\)', output)
        mapped_match = re.search(r'Mapping complete:\s+(\d+)\/(\d+) reads mapped', output)
        if not mapped_match:
            mapped_match = re.search(r'Loaded (\d+) unique read mappings', output)

        accuracy = float(acc_match.group(1)) if acc_match else None
        mapped = int(mapped_match.group(1)) if mapped_match else None

        return {"accuracy": accuracy, "mapped": mapped, "raw_output": output[:500]}
    except Exception as e:
        return {"accuracy": None, "mapped": None, "error": str(e)}


# ============================
# TEST RUNNER
# ============================
def run_test(test):
    """Run a single test and return results."""
    start = time.time()
    name = test["name"]
    cmd = test["cmd"]
    sam_file = cmd[cmd.index("-o") + 1]

    print(f"[{name}] Starting: {' '.join(cmd[4:10])}...")

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600
        )
        elapsed = time.time() - start

        if result.returncode != 0:
            print(f"[{name}] FAILED after {elapsed:.0f}s: {result.stderr[:200]}")
            return {"name": name, "desc": test["desc"], "status": "failed",
                    "time": elapsed, "accuracy": None, "mapped": None}

        if not os.path.exists(sam_file):
            print(f"[{name}] No output file after {elapsed:.0f}s")
            return {"name": name, "desc": test["desc"], "status": "no_output",
                    "time": elapsed, "accuracy": None, "mapped": None}

        acc_result = run_accuracy(sam_file)
        elapsed = time.time() - start

        status = "done"
        print(f"[{name}] DONE in {elapsed:.0f}s | "
              f"Mapped: {acc_result['mapped']:,} | "
              f"Accuracy: {acc_result['accuracy']:.2f}%")

        return {
            "name": name,
            "desc": test["desc"],
            "status": status,
            "time": round(elapsed, 1),
            "accuracy": acc_result["accuracy"],
            "mapped": acc_result["mapped"],
            "sam_file": sam_file
        }

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        print(f"[{name}] TIMEOUT after {elapsed:.0f}s")
        return {"name": name, "desc": test["desc"], "status": "timeout",
                "time": elapsed, "accuracy": None, "mapped": None}
    except Exception as e:
        elapsed = time.time() - start
        print(f"[{name}] ERROR after {elapsed:.0f}s: {e}")
        return {"name": name, "desc": test["desc"], "status": "error",
                "time": elapsed, "accuracy": None, "mapped": None}


# ============================
# MAIN
# ============================
def main():
    tests = gen_tests()
    print(f"{'='*60}")
    print(f"  CAMI Low Complexity - Benchmark Suite")
    print(f"  Tests: {len(tests)}")
    print(f"  Workers: {MAX_WORKERS}")
    print(f"  Timeout: {TOTAL_TIMEOUT_HOURS}h")
    print(f"  Started: {datetime.now().strftime('%H:%M:%S')}")
    print(f"{'='*60}\n")

    results = []
    start_time = time.time()

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(run_test, t): t for t in tests}

        for future in as_completed(futures):
            elapsed_total = time.time() - start_time

            # Check total timeout
            if elapsed_total > TOTAL_TIMEOUT_HOURS * 3600:
                print(f"\n{'='*60}")
                print(f"  TIMEOUT: {TOTAL_TIMEOUT_HOURS}h reached!")
                print(f"{'='*60}")
                executor.shutdown(wait=False, cancel_futures=True)
                break

            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"Unhandled error: {e}")

    # ============================
    # RESULTS
    # ============================
    print(f"\n{'='*60}")
    print(f"  RESULTS SUMMARY")
    print(f"{'='*60}\n")

    # Sort by accuracy (descending), then by mapped (descending)
    valid = [r for r in results if r.get("accuracy")]
    valid.sort(key=lambda x: (-x["accuracy"], -x.get("mapped", 0)))

    print(f"{'Name':<35} {'Desc':<30} {'Mapped':>10} {'Accuracy':>10} {'Time':>8}")
    print("-" * 95)

    for r in valid:
        acc_str = f"{r['accuracy']:.2f}%" if r["accuracy"] else "N/A"
        mapped_str = f"{r['mapped']:,}" if r["mapped"] else "N/A"
        time_str = f"{r['time']:.0f}s"
        print(f"{r['name']:<35} {r['desc']:<30} {mapped_str:>10} {acc_str:>10} {time_str:>8}")

    print(f"\n{'='*60}")
    print(f"  TOP 10 BY ACCURACY")
    print(f"{'='*60}\n")

    for i, r in enumerate(valid[:10], 1):
        print(f"  {i}. {r['name']:<35} {r['accuracy']:.2f}%  (mapped: {r['mapped']:,})")

    # Save results
    with open(RESULTS_FILE, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {RESULTS_FILE}")
    print(f"Completed: {len(results)}/{len(tests)} tests")
    print(f"Total time: {(time.time() - start_time)/60:.1f} minutes")


if __name__ == "__main__":
    main()
