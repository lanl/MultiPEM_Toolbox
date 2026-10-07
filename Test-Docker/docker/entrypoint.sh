#!/usr/bin/env bash
set -Eeuo pipefail

test_root="${TEST_MPEM_TEST_ROOT:-/work/Test}"
cpu_profile="${TEST_MPEM_CPU_PROFILE:-/opt/testmpem/docker/cpu-profile.R}"

usage() {
  cat <<'EOF'
Usage:
  testmpem-container run TEST
  testmpem-container help

TEST is global or a path relative to Test, such as IYDT/Seismic.
EOF
}

valid_test_name() {
  local value="$1"
  [[ -n "$value" && "$value" != /* && "$value" != *\\* ]] || return 1
  [[ "/$value/" != *"/../"* && "/$value/" != *"/./"* && "$value" != *"//"* ]]
}

run_test() {
  local test_name="${1:-}"
  local test_dir
  local status_file=".testmpem.status"
  local -a pipeline_status
  local rc

  if [[ "$test_name" == "global" ]]; then
    test_dir="$test_root"
  elif valid_test_name "$test_name"; then
    test_dir="$test_root/$test_name"
  else
    printf 'Invalid test name: %s\n' "$test_name" >&2
    return 64
  fi
  if [[ ! -f "$test_dir/tests.r" ]]; then
    printf 'Unknown test suite or missing tests.r: %s\n' "$test_name" >&2
    return 66
  fi

  cd "$test_dir"
  printf 'state=running\nstarted=%s\ntest=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$test_name" > "$status_file"
  : > tests.cpu-adjustments.log

  printf '\nMultiPEM test invocation: started=%s test=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$test_name" | tee -a test-job.out

  set +e
  TEST_MPEM_CPU_LOG="$PWD/tests.cpu-adjustments.log" \
    R_PROFILE_USER="$cpu_profile" \
    Rscript --no-restore --no-save \
      /opt/testmpem/docker/run-test.R tests.r .RData \
      2>&1 | tee tests.out | tee -a test-job.out
  pipeline_status=("${PIPESTATUS[@]}")
  set -e
  rc=${pipeline_status[0]}
  if [[ $rc -eq 0 && ${pipeline_status[1]} -ne 0 ]]; then rc=${pipeline_status[1]}; fi
  if [[ $rc -eq 0 && ${pipeline_status[2]} -ne 0 ]]; then rc=${pipeline_status[2]}; fi

  printf 'state=%s\nfinished=%s\ntest=%s\nexit_code=%s\n' \
    "$([[ $rc -eq 0 ]] && printf complete || printf failed)" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$test_name" "$rc" > "$status_file"
  printf 'Finished %s with exit code %s\n' "$test_name" "$rc" | tee -a test-job.out
  return "$rc"
}

command="${1:-}"
case "$command" in
  run)
    shift
    [[ $# -eq 1 ]] || { usage >&2; exit 64; }
    run_test "$1"
    ;;
  help|-h|--help|'') usage ;;
  *)
    printf 'Unknown command: %s\n' "$command" >&2
    usage >&2
    exit 64
    ;;
esac
