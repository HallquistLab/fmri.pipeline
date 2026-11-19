function spm_write_output_manifest(outdir, write_residuals)
% Validate the fitted model and export actual SPM.xCon names and filenames.
loaded = load(fullfile(outdir, 'SPM.mat'), 'SPM');
SPM = loaded.SPM;
assert(isfield(SPM, 'Vbeta') && ~isempty(SPM.Vbeta), 'SPM estimation has no beta images.');
assert(numel(SPM.Vbeta) == size(SPM.xX.X, 2), 'SPM beta count disagrees with the design.');
assert(isfield(SPM, 'xCon') && ~isempty(SPM.xCon), 'SPM has no estimated contrasts.');
expected_names = regexp(strtrim(fileread(fullfile(outdir, 'spm_expected_contrasts.txt'))), '\r?\n', 'split');
names = {SPM.xCon.name};
assert(numel(unique(names)) == numel(names), 'Duplicate SPM contrast names.');
assert(isequal(sort(expected_names(:)), sort(names(:))), 'SPM contrasts do not match the requested specification.');
files = {fullfile(outdir, 'SPM.mat')};
for k = 1:numel(SPM.Vbeta)
    files{end+1} = image_path(outdir, SPM.Vbeta(k).fname);
end
for field = {'VM', 'VResMS'}
    if isfield(SPM, field{1}) && ~isempty(SPM.(field{1}))
        files{end+1} = image_path(outdir, SPM.(field{1}).fname);
    end
end
rows = cell(numel(SPM.xCon), 1);
for k = 1:numel(SPM.xCon)
    con = SPM.xCon(k);
    assert(~isempty(con.Vcon) && ~isempty(con.Vspm), 'Contrast %s was not estimated.', con.name);
    coef = image_path(outdir, con.Vcon.fname);
    stat = image_path(outdir, con.Vspm.fname);
    files = [files, {coef, stat}];
    rows{k} = sprintf('%d\t"%s"\t%s\t"%s"\t"%s"', k, ...
        strrep(con.name, '"', '""'), con.STAT, strrep(coef, '"', '""'), strrep(stat, '"', '""'));
end
if write_residuals
    assert(isfield(SPM, 'xY') && isfield(SPM.xY, 'VY'), 'Cannot determine expected residual images.');
    for k = 1:numel(SPM.xY.VY)
        files{end+1} = fullfile(outdir, sprintf('Res_%04d.nii', k));
    end
end
for k = 1:numel(files)
    info = dir(files{k});
    assert(~isempty(info) && info.bytes > 0, 'Missing or empty SPM output: %s', files{k});
    if ~strcmp(files{k}, fullfile(outdir, 'SPM.mat'))
        spm_vol(files{k}); % reject invalid image headers
    end
end
manifest = fullfile(outdir, 'spm_contrasts.tsv');
write_lines(manifest, [{'contrast_number	contrast_name	stat	coef_file	stat_file'}; rows]);
files = [files, {manifest, fullfile(outdir, 'spm_expected_contrasts.txt')}];
write_lines(fullfile(outdir, 'spm_expected_outputs.txt'), unique(files, 'stable'));
copyfile(fullfile(outdir, 'spm_setup_token.txt'), fullfile(outdir, 'spm_validated_token.txt'));
end

function path = image_path(outdir, path)
if isempty(regexp(path, '^(/|[A-Za-z]:[\\/])', 'once'))
    path = fullfile(outdir, path);
end
end

function write_lines(path, lines)
fid = fopen(path, 'w');
assert(fid ~= -1, 'Cannot write SPM manifest: %s', path);
closer = onCleanup(@() fclose(fid));
for k = 1:numel(lines)
    fprintf(fid, '%s\n', lines{k});
end
end
