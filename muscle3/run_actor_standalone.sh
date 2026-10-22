#!/bin/bash
# Standalone MUSCLE3 test of the linear-model actor with one model: voltage_driver.py
# applies a waveform on one power supply through the actor (actor_standalone.ymmsl),
# the test of utils/run_supply_sequence.m with one supply, through MUSCLE3.
#   run_actor_standalone.sh [model_dir] [key=value ...]
# model_dir: a model folder (default linear_models/plasmaless_model). Keys: channel,
# waveform (step | pulse | bipolar), amplitude [V], t_on [s], width [s] (driver
# settings), t_end [s], dt [s] (dt also sets t_interval); defaults in
# actor_standalone.ymmsl. Example:
#   run_actor_standalone.sh linear_models/niceSimpleLin_105073_t84.4176 channel=VS3 waveform=pulse width=0.05
# Results: <run dir>/channel_response.mat, to overlay on the MATLAB result with
# utils/run_supply_sequence.m (muscle3_results); run dir
# /scratch/users/$USER/linear_model/m3_<date> (also the MUSCLE3 logs).
# Uses 2 cores (driver + one single-threaded MATLAB): run it in Slurm, e.g.
#   srun --partition=all_debug --time=00:20:00 --cpus-per-task=2 --mem=8G muscle3/run_actor_standalone.sh ...

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export LINEAR_MODEL_REPO="$(dirname "$here")"
export OMP_NUM_THREADS=1
linear_model_dir="$LINEAR_MODEL_REPO/linear_models/plasmaless_model"
if [[ $# -gt 0 && "$1" != *=* ]]; then
    linear_model_dir="$(cd "$1" 2>/dev/null && pwd)"
    if [[ -z "$linear_model_dir" ]]; then
        echo "run_actor_standalone.sh: no folder $1" >&2
        exit 2
    fi
    shift
fi
if [[ ! -f "$linear_model_dir/linear_model.json" ]]; then
    echo "run_actor_standalone.sh: $linear_model_dir has no linear_model.json" >&2
    exit 2
fi
run_dir=/scratch/users/$USER/linear_model/m3_$(date +%Y%m%d_%H%M%S)
mkdir -p "$run_dir"
results_file="$run_dir/channel_response.mat"

# a number written as a YAML float (decimal point; YAML 1.1 reads 1e-3 as a string)
yaml_float() {
    if [[ ! "$1" =~ ^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$ ]]; then
        echo "run_actor_standalone.sh: '$1' is not a number" >&2
        return 1
    fi
    awk -v v="$1" 'BEGIN { s = sprintf("%.15g", v + 0)
        if (s ~ /[eE]/ && s !~ /\./) sub(/[eE]/, ".0e", s); else if (s !~ /[.eE]/) s = s ".0"
        print s }'
}

# settings with absolute paths (a relative path would resolve in the instance work dir)
# and the key=value overrides
{
    echo "ymmsl_version: v0.2"
    echo "settings:"
    echo "  driver.md_uri: \"imas:hdf5?path=$LINEAR_MODEL_REPO/data/md_dd4\""
    echo "  driver.results_file: \"$results_file\""
    echo "  driver.model_dir: \"$linear_model_dir\""
    echo "  plant.linear_model_dir: \"$linear_model_dir\""
    for kv in "$@"; do
        key="${kv%%=*}"
        value="${kv#*=}"
        case "$key" in
            channel|waveform)
                echo "  driver.$key: \"$value\"" ;;
            amplitude|t_on|width)
                value="$(yaml_float "$value")" || exit 2
                echo "  driver.$key: $value" ;;
            t_end)
                value="$(yaml_float "$value")" || exit 2
                echo "  t_end: $value" ;;
            dt)
                value="$(yaml_float "$value")" || exit 2
                echo "  dt: $value"
                echo "  t_interval: $value" ;;
            *)
                echo "run_actor_standalone.sh: unknown key '$key' (channel, waveform, amplitude, t_on, width, t_end, dt)" >&2
                exit 2 ;;
        esac
    done
} > "$run_dir/settings.ymmsl"

module use /work/projects/pds/modules/all
module load IMAS-MUSCLE3/1.0.0-intel-2025b-pds

muscle_manager --start-all --run-dir "$run_dir" "$here/actor_standalone.ymmsl" "$run_dir/settings.ymmsl"
rc=$?

echo "muscle_manager exit code: $rc"
echo "Model:   $linear_model_dir"
echo "Logs:    $run_dir (muscle3_manager.log, instances/*/stdout.txt, stderr.txt)"
echo "Results: $results_file"
exit $rc
