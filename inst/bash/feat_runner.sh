# Shared by immediate execution and scheduler jobs. Each model owns a sibling lock.
feat_runner() {
  fsf="$1"; odir="$2"; token="$3"; action="$4"; njobs="${5:-}"
  lock="${fsf}.lock"
  if [ ! -f "$lock/owner" ] || [ "$(cat "$lock/owner")" != "$token" ]; then
    echo "FEAT reservation is missing or belongs to another job: $fsf" >&2
    exit 1
  fi
  started=0; interrupted=0; feat_pid=
  start_time=$(date "+%Y-%m-%d %H:%M:%S")
  finish_feat() {
    exit_code=$?
    trap - EXIT
    if [ "$started" = 1 ]; then
      mkdir -p "$odir" || exit_code=1
      if [ "$exit_code" -eq 0 ]; then
        status_file="$odir/.feat_complete"
        artifact_status=COMPLETED
      else
        status_file="$odir/.feat_fail"
        artifact_status=FAILED
        rm -f "$odir/.feat_complete" "$odir/.feat_validated.rds"
      fi
      if ! printf '%s\n%s\n' "$start_time" "$(date "+%Y-%m-%d %H:%M:%S")" > "${status_file}.tmp" ||
          ! mv "${status_file}.tmp" "$status_file"; then
        rm -f "$odir/.feat_complete"
        exit_code=1
        artifact_status=FAILED
      fi
      write_job_manifest "$odir" "fsl_l${feat_level}_feat" "$artifact_status" "$fsf"
    fi
    # Interrupted FSL descendants may still be active: retain the lock conservatively.
    if [ "$interrupted" = 0 ] && [ "$(cat "$lock/owner" 2>/dev/null)" = "$token" ]; then
      rm -f "$lock/owner"
      rmdir "$lock"
    fi
    exit "$exit_code"
  }
  interrupt_feat() {
    interrupted=1
    if [ -n "$feat_pid" ]; then kill -TERM "$feat_pid" 2>/dev/null; wait "$feat_pid"; fi
    exit 143
  }
  trap finish_feat EXIT
  trap interrupt_feat TERM INT HUP
  Rscript "$feat_integrity_script" prepare "$fsf" "$action" || exit $?
  started=1
  if [ "${feat_use_retry:-0}" = 1 ]; then
    "$feat_retry_binary" "$feat_binary" "$fsf" "$njobs" &
  elif [ -n "$njobs" ]; then
    "$feat_binary" "$fsf" -P "$njobs" &
  else
    "$feat_binary" "$fsf" &
  fi
  feat_pid=$!
  wait "$feat_pid"
  exit_code=$?
  feat_pid=
  if [ -n "${feat_fallback_logger:-}" ]; then
    Rscript "$feat_fallback_logger" --feat_dir "$odir" --level "$feat_level" \
      --log_txt "${feat_log_txt:-NULL}" --log_json "${feat_log_json:-NULL}" \
      --threshold "${feat_log_threshold:-warn}" || true
  fi
  [ "$exit_code" -eq 0 ] || exit "$exit_code"
  if [ "$feat_level" = 1 ]; then bash "$feat_reg_script" "$odir" || exit $?; fi
  Rscript "$feat_integrity_script" validate "$fsf" || exit $?
}
