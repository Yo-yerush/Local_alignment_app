# Sequence Alignment Studio

A local JavaScript app for exploratory alignment of DNA, RNA, and protein FASTA sequences.

Use the web app directly at: <https://yo-yerush.github.io/Local_alignment_app/>

Open `index.html` in any modern browser. No installation, Node.js, R, Python, server, or internet connection is required.

## What it does

- Reads pasted FASTA or uploaded FASTA files.
- Auto-detects DNA, RNA, or protein sequence type, with manual override.
- Performs:
  - Pairwise alignment of the first two sequences.
  - Multiple sequence alignment using a simple progressive center-star method.
- Supports alignment algorithms:
  - Global / Needleman-Wunsch.
  - Local / Smith-Waterman.
  - Semi-global / overlap alignment.
- Supports affine gap penalties with gap-open and gap-extension values.
- Includes scoring models:
  - DNA/RNA simple identity.
  - DNA/RNA IUPAC ambiguity-aware matching.
  - DNA/RNA transition/transversion scoring.
  - Protein simple identity.
  - Protein biochemical similarity groups.
  - BLOSUM62.
  - PAM250.
- Shows a colored alignment plot with:
  - Nucleotide coloring.
  - Protein Clustal-like coloring.
  - Protein Zappo-style groups.
  - Protein Taylor-like groups.
  - Conservation-only coloring.
  - Conservation black/white and monochrome modes.
  - Optional sequence logo row, consensus row, similarity marks, and column numbers.
  - Customizable colors, box spacing/roundness, output font, and output font size.
- Exports:
  - Aligned FASTA.
  - Alignment plot as SVG, including wrapping, optional logo/consensus/similarity rows, row-end positions, and legend.
- Supports custom structural annotations & highlights:
  - **Range Tracks**: Add empty spacer rows above or below the alignment to draw structural shapes (Helix Wave `〰️`, Loop Arch `⅏`, Cylinder Box `▱`, Helix-Loop-Helix combo `🧬`, Curly Brackets, Straight Lines, Square Brackets, or Arrow Lines) with centered text labels and optional translucent column highlights.
  - **Region Overlays**: Highlight specific columns, residues, or rows with background blocks, borders, or divider lines without drawing text labels on the plot.

## Suggested use

Use this app for small-to-medium exploratory alignments, teaching, checking PCR fragments, comparing protein family members, checking conserved domains, or preparing simple alignment figures.

For many sequences, very long sequences, or publication-grade multiple sequence alignment, validate with tools such as MAFFT, MUSCLE, Clustal Omega, EMBOSS needle/water, DECIPHER, Biostrings, or msa in R.

## Algorithm/model guide

- **Global / Needleman-Wunsch**: best when sequences should align from end to end and are similar in length.
- **Local / Smith-Waterman**: best when only part of the sequence is shared, such as domains, motifs, or fragments.
- **Semi-global / overlap**: best for partial sequences, reads, contigs, or incomplete ends.
- **DNA/RNA simple identity**: best for closely related nucleotide sequences.
- **IUPAC ambiguity-aware**: best when bases include N, R, Y, W, S, K, M, B, D, H, V.
- **Transition/transversion**: useful for evolutionary nucleotide comparisons.
- **Protein biochemical similarity**: useful for teaching or rapid exploratory protein comparisons.
- **BLOSUM62**: good default protein model for moderate similarity.
- **PAM250**: more permissive; useful for more divergent proteins.

## Files

- `index.html` — app layout.
- `src/style.css` — interface and alignment plot styling.
- `src/app.js` — parser, alignment algorithms, scoring models, viewer, and export logic.
- `examples_dna.fasta` — small DNA demo.
- `examples_protein.fasta` — small protein demo.

## Notes and limitations

The multiple sequence alignment is a simple progressive center-star implementation. It is useful and transparent, but it is not intended to replace mature MSA programs for difficult biological inference.

The SVG export is best for short to medium alignments. Very long alignments can produce large SVG files, although wrapped export keeps the plot from becoming one extremely wide row.

## Authors

Developed by **Yonatan Yerushalmy** (Rachel Amir's group).
GitHub Repository: [Local_alignment_app](https://github.com/Yo-yerush/Local_alignment_app)
