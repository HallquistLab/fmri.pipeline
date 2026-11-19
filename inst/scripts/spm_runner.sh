# Embedded in the generated scheduler script; all stages share one failure path.
function spm_runner() {
  odir="$1"
  rm -f "${odir}/.spm_fail" "${odir}/.spm_complete" "${odir}/.spm_stage_complete" \
    "${odir}/spm_validated_token.txt"
  start_time=$(date)

  run_matlab() {
    local script_path="$1"
    if [ ! -s "$script_path" ]; then
      echo "Missing SPM script: $script_path" >&2
      return 1
    fi
    # MATLAB literals escape an apostrophe by doubling it.
    local escaped_path="${script_path//\'/\'\'}"
    local cmd_str="try; run('${escaped_path}'); catch ME; disp(getReport(ME,'extended')); exit(1); end; exit(0);"
    "$matlab_cmd" $matlab_args "$cmd_str"
  }

  spm_stages() {
    if [ -f "${odir}/gunzip_commands.sh" ]; then
      bash "${odir}/gunzip_commands.sh" || return $?
    fi
    if [ "$run_setup" -eq 1 ]; then
      run_matlab "${odir}/${spm_setup_script}" || return $?
    fi
    if [ "$run_glm" -eq 1 ]; then
      run_matlab "${odir}/${spm_glm_script}" || return $?
    fi
    if [ "$run_contrasts" -eq 1 ]; then
      if [ "$spm_level" -eq 1 ]; then
        if [ ! -s "${odir}/setup_spm_contrasts.sh" ]; then
          echo "Missing SPM contrast mapping script in $odir" >&2
          return 1
        fi
        bash "${odir}/setup_spm_contrasts.sh" || return $?
      fi
      run_matlab "${odir}/${spm_contrast_script}" || return $?
    fi
    if [ "$run_setup" -eq 1 ] && [ "$run_glm" -eq 1 ] && [ "$run_contrasts" -eq 1 ]; then
      run_matlab "${odir}/validate_spm_outputs.m" || return $?
      [ -s "${odir}/spm_expected_outputs.txt" ] || return 1
      cmp -s "${odir}/spm_setup_token.txt" "${odir}/spm_validated_token.txt" || return 1
      while IFS= read -r expected; do
        [ -s "$expected" ] || { echo "Missing SPM output: $expected" >&2; return 1; }
      done < "${odir}/spm_expected_outputs.txt"
    fi
    return 0
  }

  spm_stages
  exit_code=$?
  end_time=$(date)
  if [ "$exit_code" -ne 0 ]; then
    status_file="${odir}/.spm_fail"
    artifact_status="FAILED"
  elif [ "$run_setup" -eq 1 ] && [ "$run_glm" -eq 1 ] && [ "$run_contrasts" -eq 1 ]; then
    status_file="${odir}/.spm_complete"
    artifact_status="COMPLETED"
    if [ "$cleanup_tmp" -eq 1 ] && [ -f "${odir}/.nifti_tmpdir" ]; then
      tmpdir=$(cat "${odir}/.nifti_tmpdir")
      rm -f "${tmpdir}"/run_*.nii
      rmdir "$tmpdir" 2>/dev/null || true
    fi
  else
    status_file="${odir}/.spm_stage_complete"
    artifact_status="INCOMPLETE"
  fi
  printf '%s\n%s\n' "$start_time" "$end_time" > "$status_file"
  write_job_manifest "$odir" "spm_l${spm_level}" "$artifact_status" "$odir"
  return "$exit_code"
}
