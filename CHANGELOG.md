# Changelog

## v0.0.28 - 2026-06-06

### Changed

- Added `number_of_bins` support to `bw_track` so bigWig summaries can be
  controlled by bin count instead of only `binsize`.
- Added `fill` style support to `bw_track` and improved handling of constant
  ranges and single-bin small regions.
- Added optional `nans_to_zeros` and signal `transform` controls for bigWig
  plotting.
- Added on-the-fly bigWig `operation` support with `second_bw`, including
  expressions such as `file - second_file` and
  `log2((1 + file) / (1 + second_file))`.

## v0.0.27 - 2026-06-06

### Changed

- Improved `gene_track` row assignment so genes and labels are laid out using
  occupied genomic spans instead of simple round-robin rows.
- Made automatic gene labels more conservative in dense regions. Labels are
  shown only when they fit the available rows; dense views keep gene structure
  readable instead of forcing overlapping names.
- Added `max_labels` and `all_labels_inside` controls for gene label behavior.
- Preserved explicit `show_label` names as priority labels while suppressing
  crowded automatic labels.
- Included genes that overlap the plotted region even when their start or end
  falls outside the view.
