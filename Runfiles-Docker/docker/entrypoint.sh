#!/usr/bin/env bash
set -Eeuo pipefail

runfiles_root="${MPEM_RUNFILES_ROOT:-/work/Runfiles}"
checkpoint_root="${MPEM_CHECKPOINT_ROOT:-/checkpoints}"
cpu_profile="${MPEM_CPU_PROFILE:-/opt/multipem/docker/cpu-profile.R}"

usage() {
  cat <<'EOF'
Usage:
  multipem-container list
  multipem-container run ANALYSIS [full|calibration|event]
EOF
}

list_analyses() {
  find "$runfiles_root" -type f \( -name runMPEM.r -o -name runMPEM_0.r \) \
    -printf '%h\n' |
    sed "s#^${runfiles_root}/##" |
    sort -u
}

run_deck() {
  local script="$1"
  local output="${script%.r}.out"
  local cpu_log="${script%.r}.cpu-adjustments.log"
  local status_file=".mpem-${script%.r}.status"
  local -a pipeline_status
  local rc

  if [[ ! -f "$script" ]]; then
    printf 'Missing input deck: %s/%s\n' "$PWD" "$script" >&2
    return 66
  fi

  printf 'state=running\nstarted=%s\nscript=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$script" > "$status_file"
  printf 'Starting %s in %s\n' "$script" "$PWD"

  # --restore provides the workspace handoff required by rapid-event decks.
  # The input decks retain complete control of seeds, algorithms, parallel
  # plans and saves. The profile limits only the number of concurrent workers.
  : > "$cpu_log"
  set +e
  MPEM_CPU_LOG="$PWD/$cpu_log" \
    R_PROFILE_USER="$cpu_profile" \
    Rscript --restore "$script" 2>&1 | tee "$output"
  pipeline_status=("${PIPESTATUS[@]}")
  rc=${pipeline_status[0]}
  if [[ $rc -eq 0 && ${pipeline_status[1]} -ne 0 ]]; then
    rc=${pipeline_status[1]}
  fi
  set -e

  printf 'state=%s\nfinished=%s\nscript=%s\nexit_code=%s\n' \
    "$([[ $rc -eq 0 ]] && printf complete || printf failed)" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$script" "$rc" > "$status_file"
  printf 'Finished %s with exit code %s\n' "$script" "$rc"
  return "$rc"
}

save_calibration_checkpoint() {
  local checkpoint_file="$checkpoint_root/calibration.RData"
  local temporary_file="$checkpoint_file.tmp.$$"
  local checksum_file="$checkpoint_root/calibration.RData.sha256"
  local checksum_temporary="$checksum_file.tmp.$$"

  if [[ ! -f .RData ]]; then
    if [[ -f runMPEM_0.r ]]; then
      printf 'Calibration completed without creating %s/.RData.\n' "$PWD" >&2
      return 66
    fi
    return 0
  fi

  mkdir -p "$checkpoint_root"
  cp .RData "$temporary_file"
  chmod 0444 "$temporary_file"
  mv -f "$temporary_file" "$checkpoint_file"
  (
    cd "$checkpoint_root"
    if command -v sha256sum >/dev/null 2>&1; then
      sha256sum calibration.RData > "${checksum_temporary##*/}"
    elif command -v shasum >/dev/null 2>&1; then
      shasum -a 256 calibration.RData > "${checksum_temporary##*/}"
    else
      printf 'No SHA-256 utility is available.\n' >&2
      exit 69
    fi
    mv -f "${checksum_temporary##*/}" "${checksum_file##*/}"
  ) || return $?
  printf 'Saved pristine calibration checkpoint to %s\n' "$checkpoint_file"
}

restore_calibration_checkpoint() {
  local checkpoint_file="$checkpoint_root/calibration.RData"
  local checksum_file="$checkpoint_root/calibration.RData.sha256"
  local temporary_file=".RData.checkpoint.$$"

  if [[ ! -f "$checkpoint_file" ]]; then
    printf 'Event stage requires pristine checkpoint %s.\n' "$checkpoint_file" >&2
    return 66
  fi

  if [[ ! -f "$checksum_file" ]]; then
    printf 'Event stage requires checkpoint checksum %s.\n' "$checksum_file" >&2
    return 66
  fi
  if ! (
    cd "$checkpoint_root"
    if command -v sha256sum >/dev/null 2>&1; then
      sha256sum -c calibration.RData.sha256 >/dev/null
    elif command -v shasum >/dev/null 2>&1; then
      shasum -a 256 -c calibration.RData.sha256 >/dev/null
    else
      printf 'No SHA-256 utility is available.\n' >&2
      exit 69
    fi
  ); then
    printf 'Calibration checkpoint checksum verification failed.\n' >&2
    return 65
  fi

  cp "$checkpoint_file" "$temporary_file"
  chmod 0600 "$temporary_file"
  mv -f "$temporary_file" .RData
  printf 'Restored working .RData from pristine calibration checkpoint.\n'
}

run_stage() {
  local analysis="$1"
  local stage="$2"

  printf '\nMultiPEM job invocation: started=%s analysis=%s stage=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$analysis" "$stage"

  case "$stage" in
    calibration)
      run_deck runMPEM.r || return $?
      save_calibration_checkpoint
      ;;
    event)
      restore_calibration_checkpoint || return $?
      run_deck runMPEM_0.r
      ;;
    full)
      run_deck runMPEM.r || return $?
      save_calibration_checkpoint || return $?
      if [[ -f runMPEM_0.r ]]; then
        restore_calibration_checkpoint || return $?
        run_deck runMPEM_0.r
      fi
      ;;
  esac
}

run_analysis() {
  local analysis="${1:-}"
  local stage="${2:-full}"
  local analysis_dir
  local -a pipeline_status
  local rc

  if [[ -z "$analysis" || "$analysis" == /* || "$analysis" == *..* ]]; then
    printf 'Invalid analysis path: %s\n' "$analysis" >&2
    return 64
  fi
  case "$stage" in
    full|calibration|event) ;;
    *) printf 'Invalid stage: %s\n' "$stage" >&2; return 64 ;;
  esac

  analysis_dir="$runfiles_root/$analysis"
  if [[ ! -d "$analysis_dir" ]]; then
    printf 'Unknown analysis: %s\n' "$analysis" >&2
    return 66
  fi
  cd "$analysis_dir"

  set +e
  run_stage "$analysis" "$stage" 2>&1 | tee -a multipem-job.out
  pipeline_status=("${PIPESTATUS[@]}")
  set -e
  rc=${pipeline_status[0]}
  if [[ $rc -eq 0 && ${pipeline_status[1]} -ne 0 ]]; then
    rc=${pipeline_status[1]}
  fi
  return "$rc"
}

command="${1:-}"
case "$command" in
  list)
    list_analyses
    ;;
  run)
    shift
    run_analysis "$@"
    ;;
  help|-h|--help|'')
    usage
    ;;
  *)
    printf 'Unknown command: %s\n' "$command" >&2
    usage >&2
    exit 64
    ;;
esac
