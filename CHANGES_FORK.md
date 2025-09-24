# STARsolo Fork – Changes Since Upstream

This fork extends STARsolo with a two-pass workflow that allows CB/UB tags to be added to unsorted BAM output without requiring coordinate sorting. Key changes:

- **New CLI flag – `--soloAddTagsToUnsorted`** (default `no`): when set to `yes` alongside `--outSAMtype BAM Unsorted`, STAR captures alignments in a temporary container, performs Solo correction, then injects the corrected CB/UB tags into the final unsorted BAM.
- **Pass-1 capture**: unsorted BAM output is redirected to a new `BAMoutputSoloTmp` writer that stores `[size][payload][trailer]` triplets containing the alignment record and `iReadAll` identifier.
- **Pass-2 replay**: a new module, `BAMunsortedAddSoloTags`, replays the temporary file, invokes `SoloFeature::addBAMtags`, and writes the final unsorted BAM while preserving tag semantics identical to the coordinate-sorted pipeline.
- **Shared plumbing updates**: parameter parsing, SAM header emission, and STAR’s orchestration flow were updated to open/close the appropriate streams based on the new flag. Buffer guardrails were added to prevent oversized records from overflowing fixed scratch buffers.
- **Testing harness**: the `testing/` directory now contains scripts that reproduce the sorted baseline, run the new mode, and compare BAM tags, alignment fields, Solo matrices, and logs to guarantee behavioural parity.

Refer to `docs/two_pass_unsorted_usage.md` for usage instructions and to `docs/two_pass_unsorted_tech.md` for implementation details.
