#!/usr/bin/env bash
#
# End-to-end demonstration of issues #26 and #27 and the fix on this branch.
#
# Runs the two MWE harnesses twice:
#   1. POST-FIX: current working tree (the branch's fix). Asserts the fix
#      works — PHARMCAT_VCFPREPROCESSOR workdir has no symlinks; per-sample
#      BCFTOOLS_VIEW gets the right target-pass BED.
#   2. PRE-FIX:  same MWEs against `git stash`-reverted versions of the
#      two changed files. Asserts the bugs reproduce — staged inputs are
#      symlinks back to source (=> htslib writes propagate to shared
#      source .fai); BCFTOOLS_VIEW silently mispairs samples and BEDs.
#
# The stash is popped on exit (including on Ctrl-C / errors) so this script
# never leaves the working tree in a half-state.
#
# Requirements: bash, git, pixi (https://pixi.sh).
# All other tooling — nextflow, samtools, htslib, bcftools — is provisioned
# on the fly by `pixi exec`.

set -euo pipefail

cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
INPUTS_DIR="$PWD/inputs"
PATCH_FILES=(
    "conf/modules/pharmcat_vcf_processing.config"
    "subworkflows/local/pharmcat_vcf_processing.nf"
)
REPO_ROOT="$(git rev-parse --show-toplevel)"
WORK_BASE="$PWD/.mwe-work"
stashed=0
ok=0; fail=0

c_red='\033[31m'; c_grn='\033[32m'; c_yel='\033[33m'; c_blu='\033[34m'; c_off='\033[0m'

restore_stash() {
    if [[ "${stashed}" == "1" ]]; then
        echo
        echo -e "${c_blu}Restoring the fix from stash...${c_off}"
        (cd "${REPO_ROOT}" && git stash pop --quiet) \
            || echo -e "${c_red}WARN: stash pop failed; run \`git stash list\` and recover manually${c_off}"
    fi
}
trap restore_stash EXIT

ensure_inputs() {
    mkdir -p "${INPUTS_DIR}"
    for f in sampleA.vcf.gz sampleA.vcf.gz.tbi sampleB.vcf.gz sampleB.vcf.gz.tbi \
             A_target.pass.bed B_target.pass.bed \
             reference.fna.bgz reference.fna.bgz.fai \
             positions.vcf.gz positions.vcf.gz.tbi \
             uni.vcf.gz uni.vcf.gz.tbi; do
        : > "${INPUTS_DIR}/${f}"
    done
}

run_nf() {
    local script="$1" wd="$2"
    rm -rf "${wd}" "${wd}.nextflow"
    pixi exec --spec 'nextflow=25.10' --spec 'htslib' --spec 'bcftools' -- \
        nextflow run "${script}" -c mwe.config -stub \
            -work-dir "${wd}" --outdir "${wd}/out" 2>&1 \
        | tail -n 5
}

find_vcfpreprocessor_workdir() {
    local wd="$1"
    find "${wd}" -mindepth 2 -maxdepth 2 -type d 2>/dev/null | while read -r d; do
        if grep -lq 'pharmcat_vcf_preprocessor\|missing_pgx_var' "$d/.command.sh" 2>/dev/null; then
            echo "$d"; break
        fi
    done
}

list_bcftools_view_pairings() {
    # echoes lines: "sample bed"
    local wd="$1"
    find "${wd}" -mindepth 2 -maxdepth 2 -type d 2>/dev/null | while read -r d; do
        if grep -q 'gzip > sample' "$d/.command.sh" 2>/dev/null; then
            local sample bed
            sample=$(grep -oE 'sample[AB]' "$d/.command.sh" | head -1)
            bed=$(ls "$d" | grep '\.bed$' | head -1)
            echo "${sample} ${bed}"
        fi
    done
}

assert() {
    local desc="$1" want="$2" got="$3"
    if [[ "${want}" == "${got}" ]]; then
        echo -e "  ${c_grn}PASS${c_off}  ${desc}: got '${got}'"
        ok=$((ok+1))
    else
        echo -e "  ${c_red}FAIL${c_off}  ${desc}: want '${want}', got '${got}'"
        fail=$((fail+1))
    fi
}

run_phase() {
    local label="$1" expect_symlinks="$2" expect_pair_a="$3" expect_pair_b="$4"

    echo
    echo -e "${c_blu}=== ${label} ===${c_off}"

    echo "--- #26 (test_26_stage.nf) ---"
    run_nf test_26_stage.nf "${WORK_BASE}/26"
    local vcf_d
    vcf_d=$(find_vcfpreprocessor_workdir "${WORK_BASE}/26")
    if [[ -z "${vcf_d}" ]]; then
        echo -e "  ${c_red}FAIL${c_off}  could not locate PHARMCAT_VCFPREPROCESSOR workdir"
        fail=$((fail+1)); return
    fi
    local symlinks
    symlinks=$(find "${vcf_d}" -maxdepth 1 -type l | wc -l)
    if [[ "${expect_symlinks}" == "zero" ]]; then
        assert "#26 input staging is by copy (0 symlinks)" "0" "${symlinks}"
    else
        if [[ "${symlinks}" -gt 0 ]]; then
            echo -e "  ${c_grn}PASS${c_off}  #26 pre-fix bug reproduces: ${symlinks} symlinks back to source"
            ok=$((ok+1))
            echo "        e.g. $(find "${vcf_d}" -maxdepth 1 -type l -printf '%f -> %l\n' | head -1)"
        else
            echo -e "  ${c_red}FAIL${c_off}  expected symlinks but found 0"
            fail=$((fail+1))
        fi
    fi

    echo "--- #27 (test_27_drift.nf) ---"
    run_nf test_27_drift.nf "${WORK_BASE}/27"
    local pair_a pair_b
    pair_a=$(list_bcftools_view_pairings "${WORK_BASE}/27" | awk '$1=="sampleA"{print $2}')
    pair_b=$(list_bcftools_view_pairings "${WORK_BASE}/27" | awk '$1=="sampleB"{print $2}')
    assert "#27 sampleA's BCFTOOLS_VIEW gets bed" "${expect_pair_a}" "${pair_a}"
    assert "#27 sampleB's BCFTOOLS_VIEW gets bed" "${expect_pair_b}" "${pair_b}"
}

# --------------------------------------------------------------------------
# Sanity checks
# --------------------------------------------------------------------------
if ! command -v pixi >/dev/null 2>&1; then
    echo "ERROR: pixi not found. Install from https://pixi.sh" >&2
    exit 1
fi
if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "ERROR: not inside a git working tree" >&2
    exit 1
fi

ensure_inputs
rm -rf "${WORK_BASE}"

# --------------------------------------------------------------------------
# Phase 1 — POST-FIX (current working tree carries the fix)
# --------------------------------------------------------------------------
run_phase "PHASE 1 — POST-FIX (this branch's working tree)" \
    "zero"             \
    "A_target.pass.bed" \
    "B_target.pass.bed"

# --------------------------------------------------------------------------
# Phase 2 — PRE-FIX (stash the fix in the two patched files)
# --------------------------------------------------------------------------
echo
echo -e "${c_yel}Stashing the fix in:${c_off}"
for f in "${PATCH_FILES[@]}"; do echo "  ${f}"; done

stash_paths=()
for f in "${PATCH_FILES[@]}"; do stash_paths+=( "${REPO_ROOT}/${f}" ); done

if (cd "${REPO_ROOT}" && git stash push --quiet -m "mwe-demo: temporarily revert fix" -- "${PATCH_FILES[@]}"); then
    stashed=1
else
    echo "WARN: stash push reported nothing to stash. Either the working tree was clean," >&2
    echo "      or the fix was already at HEAD. Running pre-fix phase anyway." >&2
fi

run_phase "PHASE 2 — PRE-FIX (working tree reverted)" \
    "many"             \
    "B_target.pass.bed" \
    "A_target.pass.bed"

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
echo
total=$((ok + fail))
if [[ "${fail}" == "0" ]]; then
    echo -e "${c_grn}All ${ok}/${total} assertions passed:${c_off}"
    echo "  • POST-FIX: #26 staged copies (no symlinks); #27 correct per-sample BED pairing."
    echo "  • PRE-FIX:  #26 staged symlinks back to source; #27 silently mispairs samples."
    exit 0
else
    echo -e "${c_red}${fail}/${total} assertion(s) failed.${c_off}"
    exit 1
fi
