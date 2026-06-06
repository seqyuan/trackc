# Changelog

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
