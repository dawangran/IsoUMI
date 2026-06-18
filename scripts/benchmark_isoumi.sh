#!/bin/sh
set -eu

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
bin="$repo_dir/src/isoumi"
input=
out_dir=
dataset=
threads=${ISOUMI_THREADS:-4}
buckets=${ISOUMI_BUCKETS:-64}
case_set=standard

usage() {
  cat >&2 <<EOF
Usage: sh scripts/benchmark_isoumi.sh [OPTIONS] <input.bam|input.sam> <out-dir>

Options:
  --dataset <NAME>   Dataset label written to metrics.tsv (default: input basename)
  --threads <INT>    IsoUMI worker threads (default: ${ISOUMI_THREADS:-4})
  --buckets <INT>    IsoUMI bucket count (default: ${ISOUMI_BUCKETS:-64})
  --case-set <NAME>  standard, sensitivity, or all (default: standard)
  --bin <PATH>       IsoUMI binary path (default: $repo_dir/src/isoumi)
  -h, --help         Show this help
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --dataset)
      test $# -ge 2 || { echo "--dataset requires a value" >&2; exit 2; }
      dataset=$2
      shift 2
      ;;
    --threads)
      test $# -ge 2 || { echo "--threads requires a value" >&2; exit 2; }
      threads=$2
      shift 2
      ;;
    --buckets)
      test $# -ge 2 || { echo "--buckets requires a value" >&2; exit 2; }
      buckets=$2
      shift 2
      ;;
    --case-set)
      test $# -ge 2 || { echo "--case-set requires a value" >&2; exit 2; }
      case_set=$2
      shift 2
      ;;
    --bin)
      test $# -ge 2 || { echo "--bin requires a value" >&2; exit 2; }
      bin=$2
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "unknown option: $1" >&2
      usage
      exit 2
      ;;
    *)
      if [ -z "$input" ]; then
        input=$1
      elif [ -z "$out_dir" ]; then
        out_dir=$1
      else
        echo "unexpected argument: $1" >&2
        usage
        exit 2
      fi
      shift
      ;;
  esac
done

if [ -z "$input" ] || [ -z "$out_dir" ]; then
  usage
  exit 2
fi

case "$case_set" in
  standard|sensitivity|all) ;;
  *) echo "--case-set must be one of: standard, sensitivity, all" >&2; exit 2 ;;
esac

if [ ! -x "$bin" ]; then
  echo "IsoUMI binary not found; run make first" >&2
  exit 2
fi

if [ ! -r "$input" ]; then
  echo "input file is not readable: $input" >&2
  exit 2
fi

mkdir -p "$out_dir"
metrics="$out_dir/metrics.tsv"
dataset=${dataset:-$(basename -- "$input")}

case "$input" in
  *.sam|*.SAM) input_format=SAM ;;
  *.bam|*.BAM) input_format=BAM ;;
  *) input_format=unknown ;;
esac

file_size_bytes=$(wc -c < "$input" | tr -d ' ')
input_reads=NA
input_cells=NA
input_umis=NA

count_sam_rows() {
  awk '$0 !~ /^@/ {n++} END {print n+0}' "$1"
}

count_sam_tag_values() {
  tag=$1
  file=$2
  awk -v tag="$tag" '
    BEGIN { FS="\t"; prefix=tag ":Z:" }
    $0 !~ /^@/ {
      for (i = 12; i <= NF; i++) {
        if (index($i, prefix) == 1) {
          seen[substr($i, length(prefix) + 1)] = 1
          break
        }
      }
    }
    END {
      n = 0
      for (v in seen) n++
      print n + 0
    }
  ' "$file"
}

count_stream_tag_values() {
  tag=$1
  awk -v tag="$tag" '
    BEGIN { FS="\t"; prefix=tag ":Z:" }
    {
      for (i = 12; i <= NF; i++) {
        if (index($i, prefix) == 1) {
          seen[substr($i, length(prefix) + 1)] = 1
          break
        }
      }
    }
    END {
      n = 0
      for (v in seen) n++
      print n + 0
    }
  '
}

if [ "$input_format" = "SAM" ]; then
  input_reads=$(count_sam_rows "$input")
  input_cells=$(count_sam_tag_values CB "$input")
  input_umis=$(count_sam_tag_values UR "$input")
elif command -v samtools >/dev/null 2>&1; then
  input_reads=$(samtools view -c "$input")
  input_cells=$(samtools view "$input" | count_stream_tag_values CB)
  input_umis=$(samtools view "$input" | count_stream_tag_values UR)
fi

detect_time_mode() {
  probe="${TMPDIR:-/tmp}/isoumi-time-probe.$$"
  if /usr/bin/time -f '%M' sh -c ':' >"$probe.out" 2>"$probe.err"; then
    rm -f "$probe.out" "$probe.err"
    echo gnu
    return
  fi
  if /usr/bin/time -l sh -c ':' >"$probe.out" 2>"$probe.err"; then
    rm -f "$probe.out" "$probe.err"
    echo bsd
    return
  fi
  rm -f "$probe.out" "$probe.err"
  echo none
}

time_mode=$(detect_time_mode)
version=$("$bin" --version 2>&1)

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  dataset label version input_path input_format input_size_bytes input_reads input_cell_barcodes input_umis \
  threads buckets started_at_utc ended_at_utc elapsed_seconds peak_rss_kb peak_rss_mb exit_status mode_params command \
  > "$metrics"

timed_status=0
timed_elapsed=0
timed_started_at=
timed_ended_at=
timed_peak_rss_kb=NA
timed_peak_rss_mb=NA

run_timed() {
  label=$1
  shift
  stdout_log="$out_dir/$label.stdout.log"
  stderr_log="$out_dir/$label.stderr.log"
  time_log="$out_dir/$label.time.txt"
  timed_started_at=$(date -u "+%Y-%m-%dT%H:%M:%SZ")
  start=$(date +%s)

  set +e
  if [ "$time_mode" = "gnu" ]; then
    /usr/bin/time -f 'peak_rss_kb	%M' -o "$time_log" "$@" >"$stdout_log" 2>"$stderr_log"
    timed_status=$?
    timed_peak_rss_kb=$(awk -F '\t' '$1 == "peak_rss_kb" {print $2; found=1} END {if (!found) print "NA"}' "$time_log")
  elif [ "$time_mode" = "bsd" ]; then
    /usr/bin/time -l "$@" >"$stdout_log" 2>"$stderr_log"
    timed_status=$?
    timed_peak_rss_kb=$(awk '/maximum resident set size/ {printf "%d\n", int(($1 + 1023) / 1024); found=1} END {if (!found) print "NA"}' "$stderr_log")
    printf 'peak_rss_kb\t%s\n' "$timed_peak_rss_kb" > "$time_log"
  else
    "$@" >"$stdout_log" 2>"$stderr_log"
    timed_status=$?
    timed_peak_rss_kb=NA
    printf 'peak_rss_kb\tNA\n' > "$time_log"
  fi
  set -e

  end=$(date +%s)
  timed_ended_at=$(date -u "+%Y-%m-%dT%H:%M:%SZ")
  timed_elapsed=$((end - start))
  if [ "$timed_peak_rss_kb" = "NA" ]; then
    timed_peak_rss_mb=NA
  else
    timed_peak_rss_mb=$(awk -v kb="$timed_peak_rss_kb" 'BEGIN {printf "%.2f", kb / 1024.0}')
  fi
}

run_case() {
  label=$1
  shift
  prefix="$out_dir/$label"
  mode_params=$*
  command="$bin --bam $input --out $prefix --threads $threads --buckets $buckets $mode_params"
  printf '%s\n' "Running $label: $mode_params" >&2
  run_timed "$label" "$bin" --bam "$input" --out "$prefix" --threads "$threads" --buckets "$buckets" "$@"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$dataset" "$label" "$version" "$input" "$input_format" "$file_size_bytes" "$input_reads" "$input_cells" "$input_umis" \
    "$threads" "$buckets" "$timed_started_at" "$timed_ended_at" "$timed_elapsed" "$timed_peak_rss_kb" "$timed_peak_rss_mb" \
    "$timed_status" "$mode_params" "$command" \
    >> "$metrics"
  if [ "$timed_status" -ne 0 ]; then
    echo "benchmark case failed: $label (see $out_dir/$label.stderr.log)" >&2
    exit "$timed_status"
  fi
}

run_standard_cases() {
  run_case default --emit-tsv --emit-explain
  run_case endbin50 --end-bin 50 --emit-tsv --emit-explain
  run_case conservative --ham 1 --ratio 0.05 --min-merge-confidence 0.60 --emit-tsv --emit-explain
  run_case baseline_no_structure --no-structure --emit-tsv --emit-explain
}

run_sensitivity_cases() {
  run_case ratio_0_05 --ratio 0.05 --emit-tsv --emit-explain
  run_case ratio_0_10 --ratio 0.10 --emit-tsv --emit-explain
  run_case ratio_0_30 --ratio 0.30 --emit-tsv --emit-explain
  run_case ratio_0_60 --ratio 0.60 --emit-tsv --emit-explain
  run_case endbin25 --end-bin 25 --emit-tsv --emit-explain
  run_case endbin100 --end-bin 100 --emit-tsv --emit-explain
  run_case confidence_0_00 --min-merge-confidence 0.00 --emit-tsv --emit-explain
  run_case confidence_0_45 --min-merge-confidence 0.45 --emit-tsv --emit-explain
  run_case confidence_0_60 --min-merge-confidence 0.60 --emit-tsv --emit-explain
}

case "$case_set" in
  standard)
    run_standard_cases
    ;;
  sensitivity)
    run_sensitivity_cases
    ;;
  all)
    run_standard_cases
    run_sensitivity_cases
    ;;
esac

printf '%s\n' "Wrote benchmark metrics to $metrics"
