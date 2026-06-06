/*
  Sequence Alignment Studio
  Self-contained JavaScript implementation for local browser use.
  Algorithms: Needleman-Wunsch, Smith-Waterman, semi-global with affine gaps,
  plus progressive center-star multiple sequence alignment.
*/

const AA = "ARNDCQEGHILKMFPSTWYV";
const DNA_ALPHABET = "ACGTRYSWKMBDHVN";
const RNA_ALPHABET = "ACGURYSWKMBDHVN";

const IUPAC_DNA = {
  A: new Set(["A"]), C: new Set(["C"]), G: new Set(["G"]), T: new Set(["T"]), U: new Set(["T"]),
  R: new Set(["A", "G"]), Y: new Set(["C", "T"]), S: new Set(["G", "C"]), W: new Set(["A", "T"]),
  K: new Set(["G", "T"]), M: new Set(["A", "C"]), B: new Set(["C", "G", "T"]), D: new Set(["A", "G", "T"]),
  H: new Set(["A", "C", "T"]), V: new Set(["A", "C", "G"]), N: new Set(["A", "C", "G", "T"]),
  "-": new Set([])
};

const PROTEIN_GROUPS = {
  hydrophobic: new Set("VILM".split("")),
  aromatic: new Set("FWY".split("")),
  positive: new Set("KRH".split("")),
  negative: new Set("DE".split("")),
  polar: new Set("STNQ".split("")),
  small: new Set("GAP".split("")),
  sulfur: new Set("CM".split("")),
  special: new Set("GP".split(""))
};

const WEAK_PROTEIN_GROUPS = [
  new Set("AVLIM".split("")),
  new Set("FYW".split("")),
  new Set("STNQ".split("")),
  new Set("KRH".split("")),
  new Set("DE".split("")),
  new Set("AGST".split("")),
  new Set("NDEQ".split("")),
  new Set("ILVF".split(""))
];

const SCORING_MODELS = [
  {
    id: "simple_nt",
    types: ["dna", "rna"],
    name: "DNA/RNA simple identity",
    short: "Fast default for closely related nucleotide sequences.",
    explanation: "Scores exact base matches with your Match score and all mismatches with your Mismatch score. Use for PCR products, variants, short genes, or sequences expected to be highly similar."
  },
  {
    id: "iupac_nt",
    types: ["dna", "rna"],
    name: "DNA/RNA IUPAC ambiguity-aware",
    short: "Handles N, R, Y, W, S, K, M, B, D, H, V partial matches.",
    explanation: "Useful when FASTA contains ambiguous bases. Exact matches get Match score; partially compatible ambiguity codes receive an intermediate score; incompatible bases receive Mismatch score."
  },
  {
    id: "transition_nt",
    types: ["dna", "rna"],
    name: "DNA/RNA transition/transversion",
    short: "Penalizes transversions more than transitions.",
    explanation: "Good for evolutionary comparison of nucleotide sequences, because A↔G and C↔T/U transitions are often more frequent than transversions."
  },
  {
    id: "protein_identity",
    types: ["protein"],
    name: "Protein simple identity",
    short: "Simple exact-match protein scoring.",
    explanation: "A quick model for very similar proteins or checking whether two protein sequences are nearly identical. It does not reward conservative amino-acid substitutions."
  },
  {
    id: "protein_biochemical",
    types: ["protein"],
    name: "Protein biochemical similarity",
    short: "Rewards conservative substitutions by amino-acid property groups.",
    explanation: "Useful for teaching and exploratory comparisons. Exact matches are best; amino acids with similar charge, polarity, aromaticity, or hydrophobicity receive intermediate scores."
  },
  {
    id: "blosum62",
    types: ["protein"],
    name: "BLOSUM62",
    short: "General-purpose protein substitution matrix.",
    explanation: "A strong default for protein alignment when sequences have moderate similarity. It is widely used for local and global protein alignment."
  },
  {
    id: "pam250",
    types: ["protein"],
    name: "PAM250",
    short: "Protein matrix for more divergent sequences.",
    explanation: "More permissive than BLOSUM62. Consider it for more distantly related proteins, but inspect alignments carefully."
  }
];

const BLOSUM62 = parseMatrix(`
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
`);

const PAM250 = parseMatrix(`
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  2 -2  0  0 -2  0  0  1 -1 -1 -2 -1 -1 -3  1  1  1 -6 -3  0
R -2  6  0 -1 -4  1 -1 -3  2 -2 -3  3  0 -4  0  0 -1  2 -4 -2
N  0  0  2  2 -4  1  1  0  2 -2 -3  1 -2 -3  0  1  0 -4 -2 -2
D  0 -1  2  4 -5  2  3  1  1 -2 -4  0 -3 -6 -1  0  0 -7 -4 -2
C -2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3  0 -2 -8  0 -2
Q  0  1  1  2 -5  4  2 -1  3 -2 -2  1 -1 -5  0 -1 -1 -5 -4 -2
E  0 -1  1  3 -5  2  4  0  1 -2 -3  0 -2 -5 -1  0  0 -7 -4 -2
G  1 -3  0  1 -3 -1  0  5 -2 -3 -4 -2 -3 -5  0  1  0 -7 -5 -1
H -1  2  2  1 -3  3  1 -2  6 -2 -2  0 -2 -2  0 -1 -1 -3  0 -2
I -1 -2 -2 -2 -2 -2 -2 -3 -2  5  2 -2  2  1 -2 -1  0 -5 -1  4
L -2 -3 -3 -4 -6 -2 -3 -4 -2  2  6 -3  4  2 -3 -3 -2 -2 -1  2
K -1  3  1  0 -5  1  0 -2  0 -2 -3  5  0 -5 -1  0  0 -3 -4 -2
M -1  0 -2 -3 -5 -1 -2 -3 -2  2  4  0  6  0 -2 -2 -1 -4 -2  2
F -3 -4 -3 -6 -4 -5 -5 -5 -2  1  2 -5  0  9 -5 -3 -3  0  7 -1
P  1  0  0 -1 -3  0 -1  0  0 -2 -3 -1 -2 -5  6  1  0 -6 -5 -1
S  1  0  1  0  0 -1  0  1 -1 -1 -3  0 -2 -3  1  2  1 -2 -3 -1
T  1 -1  0  0 -2 -1  0  0 -1  0 -2  0 -1 -3  0  1  3 -5 -3  0
W -6  2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4  0 -6 -2 -5 17  0 -6
Y -3 -4 -2 -4  0 -4 -4 -5  0 -1 -1 -4 -2  7 -5 -3 -3  0 10 -2
V  0 -2 -2 -2 -2 -2 -2 -1 -2  4  2 -2  2 -1 -1 -1  0 -6 -2  4
`);

const DOM = {
  fastaInput: document.getElementById("fastaInput"),
  fileInput: document.getElementById("fileInput"),
  clearInput: document.getElementById("clearInput"),
  inputSummary: document.getElementById("inputSummary"),
  seqType: document.getElementById("seqType"),
  taskMode: document.getElementById("taskMode"),
  algorithm: document.getElementById("algorithm"),
  scoringModel: document.getElementById("scoringModel"),
  gapOpen: document.getElementById("gapOpen"),
  gapExtend: document.getElementById("gapExtend"),
  matchScore: document.getElementById("matchScore"),
  mismatchScore: document.getElementById("mismatchScore"),
  colorScheme: document.getElementById("colorScheme"),
  wrapCols: document.getElementById("wrapCols"),
  modelExplanation: document.getElementById("modelExplanation"),
  runBtn: document.getElementById("runBtn"),
  copyFastaBtn: document.getElementById("copyFastaBtn"),
  downloadFastaBtn: document.getElementById("downloadFastaBtn"),
  downloadSvgBtn: document.getElementById("downloadSvgBtn"),
  showConsensus: document.getElementById("showConsensus"),
  showSimilarity: document.getElementById("showSimilarity"),
  resultStats: document.getElementById("resultStats"),
  legend: document.getElementById("legend"),
  viewer: document.getElementById("viewer"),
  alignedOutput: document.getElementById("alignedOutput"),
  loadDemoDna: document.getElementById("loadDemoDna"),
  loadDemoProtein: document.getElementById("loadDemoProtein"),
  colorPickerGrid: document.getElementById("colorPickerGrid"),
  resetColorsBtn: document.getElementById("resetColorsBtn"),
  boxRoundness: document.getElementById("boxRoundness"),
  boxRoundnessVal: document.getElementById("boxRoundnessVal"),
  makeWhiteBtn: document.getElementById("makeWhiteBtn"),
  boxGap: document.getElementById("boxGap"),
  boxGapVal: document.getElementById("boxGapVal"),
  annotationsList: document.getElementById("annotationsList"),
  annTypeSelect: document.getElementById("annTypeSelect"),
  annStartCol: document.getElementById("annStartCol"),
  annEndCol: document.getElementById("annEndCol"),
  annStartCol2: document.getElementById("annStartCol2"),
  annEndCol2: document.getElementById("annEndCol2"),
  annStartCol3: document.getElementById("annStartCol3"),
  annEndCol3: document.getElementById("annEndCol3"),
  annStartCol4: document.getElementById("annStartCol4"),
  annEndCol4: document.getElementById("annEndCol4"),
  annHasSeg2: document.getElementById("annHasSeg2"),
  annHasSeg3: document.getElementById("annHasSeg3"),
  annHasSeg4: document.getElementById("annHasSeg4"),
  seg2Row: document.getElementById("seg2Row"),
  seg3Row: document.getElementById("seg3Row"),
  seg4Row: document.getElementById("seg4Row"),
  annSeqSelect: document.getElementById("annSeqSelect"),
  annResPos: document.getElementById("annResPos"),
  annNameInput: document.getElementById("annNameInput"),
  annLabelInput: document.getElementById("annLabelInput"),
  annTrackLabelLabel: document.getElementById("annTrackLabelLabel"),
  annColorInput: document.getElementById("annColorInput"),
  annTextColorInput: document.getElementById("annTextColorInput"),
  addAnnBtn: document.getElementById("addAnnBtn"),
  cancelAnnBtn: document.getElementById("cancelAnnBtn"),
  rangeInputsRow: document.getElementById("rangeInputsRow"),
  residueInputsRow: document.getElementById("residueInputsRow"),
  showColumnNumbers: document.getElementById("showColumnNumbers")
};

let currentResult = null;
let annotationsArray = [];
let editingIndex = null;

const DEFAULT_COLORS = {
  // Nucleotide
  "nt-A-bg": "#bbf7d0", "nt-A-fg": "#14532d",
  "nt-C-bg": "#bfdbfe", "nt-C-fg": "#1e3a8a",
  "nt-G-bg": "#fed7aa", "nt-G-fg": "#7c2d12",
  "nt-T-bg": "#fecaca", "nt-T-fg": "#7f1d1d",
  "nt-U-bg": "#fecaca", "nt-U-fg": "#7f1d1d",
  "nt-N-bg": "#e5e7eb", "nt-N-fg": "#374151",
  
  // Protein groups
  "pro-hydrophobic-bg": "#fde68a", "pro-hydrophobic-fg": "#713f12",
  "pro-aromatic-bg": "#fbcfe8", "pro-aromatic-fg": "#831843",
  "pro-positive-bg": "#bfdbfe", "pro-positive-fg": "#1e3a8a",
  "pro-negative-bg": "#fecaca", "pro-negative-fg": "#7f1d1d",
  "pro-polar-bg": "#bbf7d0", "pro-polar-fg": "#14532d",
  "pro-small-bg": "#e9d5ff", "pro-small-fg": "#581c87",
  "pro-special-bg": "#ddd6fe", "pro-special-fg": "#4c1d95",
  "pro-sulfur-bg": "#fef08a", "pro-sulfur-fg": "#713f12",

  // General & Conservation
  "conserved-bg": "#dcfce7", "conserved-fg": "#166534",
  "simconserved-bg": "#fef9c3", "simconserved-fg": "#854d0e",
  "mono-bg": "#f8fafc", "mono-fg": "#0f172a",
  "gap-bg": "#e7e9ee", "gap-fg": "#8a94a6",

  // Consensus row
  "consensus-bg": "transparent", "consensus-fg": "#0f172a",

  // Outlines
  "identical-outline-border": "transparent",
  "similar-outline-border": "transparent"
};

let currentColors = { ...DEFAULT_COLORS };
let currentRoundness = 0;
let currentGap = 1;

function applyColors() {
  for (const [key, val] of Object.entries(currentColors)) {
    document.documentElement.style.setProperty(`--color-${key}`, val);
  }
}

function applyRoundness() {
  document.documentElement.style.setProperty(`--cell-radius`, `${currentRoundness}px`);
  DOM.boxRoundnessVal.textContent = `${currentRoundness}px`;
  DOM.boxRoundness.value = currentRoundness;
}

function applyGap() {
  document.documentElement.style.setProperty(`--cell-gap`, `${currentGap}px`);
  document.documentElement.style.setProperty(`--cell-border-color`, currentGap === 0 ? "transparent" : "rgba(0,0,0,0.04)");
  DOM.boxGapVal.textContent = `${currentGap}px`;
  DOM.boxGap.value = currentGap;
}

function parseResiduePositions(inputStr) {
  const positions = [];
  const parts = inputStr.split(",").map(p => p.trim());
  for (const part of parts) {
    if (!part) continue;
    if (part.includes("-")) {
      const [startStr, endStr] = part.split("-");
      const start = parseInt(startStr, 10);
      const end = parseInt(endStr, 10);
      if (!isNaN(start) && !isNaN(end) && start > 0 && end > 0) {
        const min = Math.min(start, end);
        const max = Math.max(start, end);
        for (let p = min; p <= max; p++) {
          positions.push(p);
        }
      }
    } else {
      const pos = parseInt(part, 10);
      if (!isNaN(pos) && pos > 0) {
        positions.push(pos);
      }
    }
  }
  return [...new Set(positions)].sort((a, b) => a - b);
}

function getAnnotationsData() {
  const rangeAnnotations = annotationsArray.filter(a => a.type === "range");
  const residueAnnotations = new Map();
  for (const ann of annotationsArray.filter(a => a.type === "residue")) {
    const positions = ann.positions || (ann.position ? [ann.position] : []);
    for (const pos of positions) {
      residueAnnotations.set(`${ann.sequence}:${pos}`, { label: ann.name || ann.label, color: ann.color });
    }
  }
  return { rangeAnnotations, residueAnnotations };
}

function updateAnnotationSeqSelect(records) {
  if (!DOM.annSeqSelect) return;
  try {
    const recs = records || parseFasta(DOM.fastaInput.value);
    DOM.annSeqSelect.innerHTML = "";
    if (recs.length === 0) {
      const opt = document.createElement("option");
      opt.value = "";
      opt.textContent = "(Load sequences first)";
      DOM.annSeqSelect.appendChild(opt);
      return;
    }
    for (const rec of recs) {
      const opt = document.createElement("option");
      opt.value = rec.name;
      opt.textContent = rec.name;
      DOM.annSeqSelect.appendChild(opt);
    }
  } catch (e) {
    DOM.annSeqSelect.innerHTML = `<option value="">(Error parsing sequences)</option>`;
  }
}

function renderAnnotationsList() {
  if (!DOM.annotationsList) return;
  DOM.annotationsList.innerHTML = "";
  
  if (annotationsArray.length === 0) {
    DOM.annotationsList.innerHTML = `<div class="summary muted" style="text-align: center; padding: 12px; font-style: italic;">No features added yet.</div>`;
    return;
  }
  
  annotationsArray.forEach((ann, idx) => {
    const item = document.createElement("div");
    item.className = "annotations-list-item";
    
    const badge = document.createElement("span");
    badge.className = `ann-badge badge-${ann.type}`;
    badge.textContent = ann.type;
    
    const desc = document.createElement("span");
    desc.className = "ann-desc";
    if (ann.type === "range") {
      let descText = "";
      let rangesText = "";
      if (ann.ranges && ann.ranges.length > 0) {
        rangesText = "Col " + ann.ranges.map(r => `${r.start}–${r.end}`).join(", ");
      } else {
        rangesText = `Col ${ann.start}–${ann.end}`;
      }
      const labelPart = ann.label ? ` [${ann.label}]` : "";
      descText = `${rangesText}: ${ann.name || ""}${labelPart}`;
      desc.textContent = descText;
    } else {
      desc.textContent = `${ann.sequence} Pos ${ann.positionInput || ann.position || ""}: ${ann.name || ann.label || ""}`;
    }
    desc.title = desc.textContent;
    
    const swatch = document.createElement("span");
    swatch.className = "ann-color-swatch";
    swatch.style.backgroundColor = ann.color;
    
    const actions = document.createElement("div");
    actions.style.display = "flex";
    actions.style.gap = "6px";
    
    // Edit Button
    const editBtn = document.createElement("button");
    editBtn.type = "button";
    editBtn.innerHTML = "✏️";
    editBtn.title = "Edit annotation";
    editBtn.style.background = "transparent";
    editBtn.style.border = "none";
    editBtn.style.cursor = "pointer";
    editBtn.style.padding = "0 4px";
    editBtn.style.fontSize = "12px";
    editBtn.addEventListener("click", () => {
      editingIndex = idx;
      DOM.addAnnBtn.textContent = "Save Feature";
      DOM.cancelAnnBtn.style.display = "inline-block";
      
      DOM.annTypeSelect.value = ann.type;
      DOM.annNameInput.value = ann.name || "";
      DOM.annLabelInput.value = ann.label || "";
      DOM.annColorInput.value = ann.color;
      
      if (ann.type === "range") {
        DOM.rangeInputsRow.style.display = "flex";
        DOM.residueInputsRow.style.display = "none";
        DOM.annTrackLabelLabel.style.display = "flex";
        DOM.annTextColorLabel.style.display = "flex";
        DOM.annTextColorInput.value = ann.textColor || "#ffffff";
        
        DOM.annHasSeg2.checked = false;
        DOM.annHasSeg3.checked = false;
        DOM.annHasSeg4.checked = false;
        DOM.seg2Row.style.display = "none";
        DOM.seg3Row.style.display = "none";
        DOM.seg4Row.style.display = "none";
        
        DOM.annStartCol.value = "";
        DOM.annEndCol.value = "";
        DOM.annStartCol2.value = "";
        DOM.annEndCol2.value = "";
        DOM.annStartCol3.value = "";
        DOM.annEndCol3.value = "";
        DOM.annStartCol4.value = "";
        DOM.annEndCol4.value = "";
        
        if (ann.ranges) {
          DOM.annStartCol.value = ann.ranges[0] ? ann.ranges[0].start : "";
          DOM.annEndCol.value = ann.ranges[0] ? ann.ranges[0].end : "";
          
          if (ann.ranges[1]) {
            DOM.annHasSeg2.checked = true;
            DOM.seg2Row.style.display = "flex";
            DOM.annStartCol2.value = ann.ranges[1].start;
            DOM.annEndCol2.value = ann.ranges[1].end;
            
            if (ann.ranges[2]) {
              DOM.annHasSeg3.checked = true;
              DOM.seg3Row.style.display = "flex";
              DOM.annStartCol3.value = ann.ranges[2].start;
              DOM.annEndCol3.value = ann.ranges[2].end;
              
              if (ann.ranges[3]) {
                DOM.annHasSeg4.checked = true;
                DOM.seg4Row.style.display = "flex";
                DOM.annStartCol4.value = ann.ranges[3].start;
                DOM.annEndCol4.value = ann.ranges[3].end;
              }
            }
          }
        } else {
          DOM.annStartCol.value = ann.start || "";
          DOM.annEndCol.value = ann.end || "";
        }
      } else {
        DOM.rangeInputsRow.style.display = "none";
        DOM.residueInputsRow.style.display = "flex";
        DOM.annTrackLabelLabel.style.display = "none";
        DOM.annTextColorLabel.style.display = "none";
        DOM.annSeqSelect.value = ann.sequence;
        DOM.annResPos.value = ann.positionInput || ann.position || "";
      }
    });
    
    // Delete Button
    const deleteBtn = document.createElement("button");
    deleteBtn.type = "button";
    deleteBtn.className = "ann-delete-btn";
    deleteBtn.innerHTML = "&times;";
    deleteBtn.title = "Delete annotation";
    deleteBtn.addEventListener("click", () => {
      if (editingIndex === idx) {
        cancelEditMode();
      } else if (editingIndex !== null && editingIndex > idx) {
        editingIndex--;
      }
      annotationsArray.splice(idx, 1);
      onAnnotationsUpdated();
    });
    
    actions.append(editBtn, deleteBtn);
    item.append(badge, desc, swatch, actions);
    DOM.annotationsList.appendChild(item);
  });
}

function onAnnotationsUpdated() {
  renderAnnotationsList();
  if (currentResult) {
    const { rangeAnnotations, residueAnnotations } = getAnnotationsData();
    currentResult.rangeAnnotations = rangeAnnotations;
    currentResult.residueAnnotations = residueAnnotations;
    renderAlignmentViewer(currentResult);
  }
}

function cancelEditMode() {
  editingIndex = null;
  DOM.addAnnBtn.textContent = "Add Feature";
  DOM.cancelAnnBtn.style.display = "none";
  
  DOM.annNameInput.value = "";
  DOM.annLabelInput.value = "";
  DOM.annStartCol.value = "";
  DOM.annEndCol.value = "";
  DOM.annStartCol2.value = "";
  DOM.annEndCol2.value = "";
  DOM.annStartCol3.value = "";
  DOM.annEndCol3.value = "";
  DOM.annStartCol4.value = "";
  DOM.annEndCol4.value = "";
  DOM.annResPos.value = "";
  DOM.annColorInput.value = "#315efb";
  DOM.annTextColorInput.value = "#ffffff";
  
  DOM.annHasSeg2.checked = false;
  DOM.annHasSeg3.checked = false;
  DOM.annHasSeg4.checked = false;
  DOM.seg2Row.style.display = "none";
  DOM.seg3Row.style.display = "none";
  DOM.seg4Row.style.display = "none";
}

function formatSimpleLabel(W, label) {
  if (W <= 0) return "";
  let text = label || "";
  if (text.length > W) {
    text = text.slice(0, W);
  }
  const padTotal = W - text.length;
  const padLeft = Math.floor(padTotal / 2);
  const padRight = padTotal - padLeft;
  return " ".repeat(padLeft) + text + " ".repeat(padRight);
}

const SCHEME_CATEGORIES = {
  nucleotide: [
    { key: "nt-A", label: "A (Adenine)", preview: "A", defaults: { bg: "nt-A-bg", fg: "nt-A-fg" } },
    { key: "nt-C", label: "C (Cytosine)", preview: "C", defaults: { bg: "nt-C-bg", fg: "nt-C-fg" } },
    { key: "nt-G", label: "G (Guanine)", preview: "G", defaults: { bg: "nt-G-bg", fg: "nt-G-fg" } },
    { key: "nt-T", label: "T (Thymine)", preview: "T", defaults: { bg: "nt-T-bg", fg: "nt-T-fg" } },
    { key: "nt-U", label: "U (Uracil)", preview: "U", defaults: { bg: "nt-U-bg", fg: "nt-U-fg" } },
    { key: "nt-N", label: "N (Ambiguous)", preview: "N", defaults: { bg: "nt-N-bg", fg: "nt-N-fg" } },
    { key: "gap", label: "Gap (-)", preview: "-", defaults: { bg: "gap-bg", fg: "gap-fg" } },
    { key: "consensus", label: "Consensus Row", preview: "C", defaults: { bg: "consensus-bg", fg: "consensus-fg" } },
    { key: "identical-outline", label: "Identical Outline", preview: "*", borderOnly: true, defaults: { border: "identical-outline-border" } },
    { key: "similar-outline", label: "Similar Outline", preview: ":", borderOnly: true, defaults: { border: "similar-outline-border" } }
  ],
  identity: [
    { key: "conserved", label: "Identical Col", preview: "C", defaults: { bg: "conserved-bg", fg: "conserved-fg" } },
    { key: "simconserved", label: "Similar Col", preview: "S", defaults: { bg: "simconserved-bg", fg: "simconserved-fg" } },
    { key: "mono", label: "Variable", preview: "V", defaults: { bg: "mono-bg", fg: "mono-fg" } },
    { key: "gap", label: "Gap (-)", preview: "-", defaults: { bg: "gap-bg", fg: "gap-fg" } },
    { key: "consensus", label: "Consensus Row", preview: "C", defaults: { bg: "consensus-bg", fg: "consensus-fg" } },
    { key: "identical-outline", label: "Identical Outline", preview: "*", borderOnly: true, defaults: { border: "identical-outline-border" } },
    { key: "similar-outline", label: "Similar Outline", preview: ":", borderOnly: true, defaults: { border: "similar-outline-border" } }
  ],
  mono: [
    { key: "mono", label: "Residue/Base", preview: "X", defaults: { bg: "mono-bg", fg: "mono-fg" } },
    { key: "gap", label: "Gap (-)", preview: "-", defaults: { bg: "gap-bg", fg: "gap-fg" } },
    { key: "consensus", label: "Consensus Row", preview: "C", defaults: { bg: "consensus-bg", fg: "consensus-fg" } },
    { key: "identical-outline", label: "Identical Outline", preview: "*", borderOnly: true, defaults: { border: "identical-outline-border" } },
    { key: "similar-outline", label: "Similar Outline", preview: ":", borderOnly: true, defaults: { border: "similar-outline-border" } }
  ],
  protein: [
    { key: "pro-hydrophobic", label: "Hydrophobic (VILM)", preview: "H", defaults: { bg: "pro-hydrophobic-bg", fg: "pro-hydrophobic-fg" } },
    { key: "pro-aromatic", label: "Aromatic (FWY)", preview: "A", defaults: { bg: "pro-aromatic-bg", fg: "pro-aromatic-fg" } },
    { key: "pro-positive", label: "Positive (KRH)", preview: "+", defaults: { bg: "pro-positive-bg", fg: "pro-positive-fg" } },
    { key: "pro-negative", label: "Negative (DE)", preview: "-", defaults: { bg: "pro-negative-bg", fg: "pro-negative-fg" } },
    { key: "pro-polar", label: "Polar (STNQ)", preview: "P", defaults: { bg: "pro-polar-bg", fg: "pro-polar-fg" } },
    { key: "pro-small", label: "Small/Special (GP)", preview: "S", defaults: { bg: "pro-small-bg", fg: "pro-small-fg" } },
    { key: "pro-sulfur", label: "Sulfur/Cys (CM)", preview: "C", defaults: { bg: "pro-sulfur-bg", fg: "pro-sulfur-fg" } },
    { key: "gap", label: "Gap (-)", preview: "-", defaults: { bg: "gap-bg", fg: "gap-fg" } },
    { key: "consensus", label: "Consensus Row", preview: "C", defaults: { bg: "consensus-bg", fg: "consensus-fg" } },
    { key: "identical-outline", label: "Identical Outline", preview: "*", borderOnly: true, defaults: { border: "identical-outline-border" } },
    { key: "similar-outline", label: "Similar Outline", preview: ":", borderOnly: true, defaults: { border: "similar-outline-border" } }
  ]
};

function renderColorCustomizer(settings) {
  if (!settings) {
    DOM.colorPickerGrid.innerHTML = `<div class="summary muted">Run an alignment to customize colors.</div>`;
    return;
  }
  
  let scheme = settings.colorScheme;
  if (scheme === "auto") scheme = settings.type === "protein" ? "clustal" : "nucleotide";
  
  let categoryKey = "mono";
  if (scheme === "nucleotide" || (settings.type !== "protein" && scheme !== "identity" && scheme !== "mono")) {
    categoryKey = "nucleotide";
  } else if (scheme === "identity") {
    categoryKey = "identity";
  } else if (scheme === "mono") {
    categoryKey = "mono";
  } else {
    categoryKey = "protein";
  }
  
  const categories = SCHEME_CATEGORIES[categoryKey];
  DOM.colorPickerGrid.innerHTML = "";
  
  for (const cat of categories) {
    const item = document.createElement("div");
    item.className = "color-picker-item";
    
    const info = document.createElement("div");
    info.className = "color-picker-info";
    
    const swatch = document.createElement("span");
    swatch.className = "color-picker-swatch";
    swatch.textContent = cat.preview;
    
    if (cat.borderOnly) {
      const borderVar = cat.defaults.border;
      swatch.style.border = `2px solid var(--color-${borderVar})`;
      swatch.style.background = `color-mix(in srgb, var(--color-${borderVar}) 22%, transparent)`;
      swatch.style.color = "var(--text)";
    } else {
      swatch.style.background = `var(--color-${cat.defaults.bg})`;
      swatch.style.color = `var(--color-${cat.defaults.fg})`;
    }
    
    const label = document.createElement("span");
    label.className = "color-picker-label";
    label.textContent = cat.label;
    label.title = cat.label;
    
    info.append(swatch, label);
    
    const inputs = document.createElement("div");
    inputs.className = "color-picker-inputs";
    
    if (cat.borderOnly) {
      const borderKey = cat.defaults.border;
      const borderWrapper = document.createElement("div");
      borderWrapper.className = "color-input-wrapper";
      
      const input = document.createElement("input");
      input.type = "color";
      let borderVal = currentColors[borderKey] || "#000000";
      if (borderVal === "transparent") borderVal = "#ffffff";
      input.value = borderVal;
      input.title = "Outline Color";
      input.addEventListener("input", e => {
        currentColors[borderKey] = e.target.value;
        applyColors();
        swatch.style.border = `2px solid ${e.target.value}`;
        swatch.style.background = `color-mix(in srgb, ${e.target.value} 22%, transparent)`;
      });
      
      const inputLabel = document.createElement("span");
      inputLabel.textContent = "Line";
      
      borderWrapper.append(input, inputLabel);
      inputs.append(borderWrapper);
    } else {
      const bgKey = cat.defaults.bg;
      const fgKey = cat.defaults.fg;
      
      const bgWrapper = document.createElement("div");
      bgWrapper.className = "color-input-wrapper";
      const bgInput = document.createElement("input");
      bgInput.type = "color";
      let bgVal = currentColors[bgKey] || "#ffffff";
      if (bgVal === "transparent") bgVal = "#ffffff";
      bgInput.value = bgVal;
      bgInput.title = "Background Color";
      bgInput.addEventListener("input", e => {
        currentColors[bgKey] = e.target.value;
        applyColors();
        swatch.style.background = e.target.value;
      });
      const bgLabel = document.createElement("span");
      bgLabel.textContent = "Bg";
      bgWrapper.append(bgInput, bgLabel);
      
      const fgWrapper = document.createElement("div");
      fgWrapper.className = "color-input-wrapper";
      const fgInput = document.createElement("input");
      fgInput.type = "color";
      fgInput.value = currentColors[fgKey] || "#000000";
      fgInput.title = "Text Color";
      fgInput.addEventListener("input", e => {
        currentColors[fgKey] = e.target.value;
        applyColors();
        swatch.style.color = e.target.value;
      });
      const fgLabel = document.createElement("span");
      fgLabel.textContent = "Text";
      fgWrapper.append(fgInput, fgLabel);
      
      inputs.append(bgWrapper, fgWrapper);
    }
    
    item.append(info, inputs);
    DOM.colorPickerGrid.appendChild(item);
  }
}

const DEMO_DNA = `>Arabidopsis_CGS_like_fragment
ATGGCTGCTACGTTGATGACTGCTGAGGAGAATCTTCTCAGGATGGGATTCGAGAACTTGCTT
>mto1_like_fragment
ATGGCTGCTACGTTGATGACTGCCGAGGAGAATCTTCTCAGGATGGGATTTGAGAACTTGCTT
>Brassica_CGS_like_fragment
ATGGCTGCTACGCTGATGACCGCTGAGGAGAACCTTCTCAGGATGGGATTCAAGAACTTGCTT`;

const DEMO_PROTEIN = `>RAB11A_human_fragment
MGTRDDEYDYLFKVVLIGDSGVGKSNLLSRFTRNEFNLESKSTIGVEFATRSIQVDGKTIKAQIWDTAGQERYRAITSAYYRGAVGALLVYDIAKHLTYENVERWLKELRDHADSNIVIMLVGNKSDLRHLRAVPTDEARAFAEKNN
>RAB11A_mouse_fragment
MGTRDDEYDYLFKVVLIGDSGVGKSNLLSRFTRNEFNLESKSTIGVEFATRSIQVDGKTIKAQIWDTAGQERYRAITSAYYRGAVGALLVYDIAKHLTYENVERWLKELRDHADSNIVIMLVGNKSDLRHLRAVPTDEARAFAEKNN
>RAB11B_fragment
MGTRDDEYDYLFKVVLIGDSGVGKSNLLSRFTRNEFNLESKSTIGVEFATRTIQVDGKTIKAQIWDTAGQERYRAITSAYYRGAVGALLVYDIAKHLTYENVERWLKELRDHADSNIVIMLVGNKSDLKHLRAVPTDEARAFAEKNN`;

function parseMatrix(text) {
  const lines = text.trim().split(/\n+/).map(l => l.trim()).filter(Boolean);
  const header = lines[0].split(/\s+/);
  const matrix = {};
  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split(/\s+/);
    const row = parts[0];
    matrix[row] = {};
    for (let j = 1; j < parts.length; j++) matrix[row][header[j - 1]] = Number(parts[j]);
  }
  return matrix;
}

function parseFasta(text) {
  const clean = text.replace(/\r/g, "").trim();
  if (!clean) return [];
  const records = [];
  let current = null;
  for (const rawLine of clean.split("\n")) {
    const line = rawLine.trim();
    if (!line || line.startsWith(";")) continue;
    if (line.startsWith(">")) {
      if (current) records.push(current);
      const fullName = line.slice(1).trim() || `sequence_${records.length + 1}`;
      current = { name: fullName, sequence: "" };
    } else {
      if (!current) current = { name: `sequence_${records.length + 1}`, sequence: "" };
      current.sequence += line.replace(/\s+/g, "").toUpperCase();
    }
  }
  if (current) records.push(current);
  return records.map((r, index) => ({ ...r, index, sequence: r.sequence.replace(/\./g, "-") }));
}

function validateSequences(records) {
  if (records.length < 2) throw new Error("Please provide at least two FASTA sequences.");
  const bad = records.filter(r => !r.sequence || !/[A-Z*-]/.test(r.sequence));
  if (bad.length) throw new Error("Every FASTA record must contain a sequence.");
  const duplicate = records.find((r, i) => records.findIndex(x => x.name === r.name) !== i);
  if (duplicate) throw new Error(`Duplicate sequence name: ${duplicate.name}. Please make names unique.`);
}

function detectType(records) {
  const letters = records.map(r => r.sequence.replace(/[-*]/g, "")).join("").toUpperCase();
  if (!letters) return "dna";
  const unique = new Set(letters.split(""));
  const dnaOk = [...unique].every(ch => DNA_ALPHABET.includes(ch));
  const rnaOk = [...unique].every(ch => RNA_ALPHABET.includes(ch));
  const hasU = unique.has("U");
  const hasT = unique.has("T");
  if (rnaOk && hasU && !hasT) return "rna";
  if (dnaOk && !hasU) return "dna";
  return "protein";
}

function availableModels(type) {
  return SCORING_MODELS.filter(m => m.types.includes(type));
}

function populateModels() {
  const records = parseFasta(DOM.fastaInput.value);
  const selectedType = DOM.seqType.value === "auto" ? detectType(records) : DOM.seqType.value;
  const models = availableModels(selectedType);
  const current = DOM.scoringModel.value;
  DOM.scoringModel.innerHTML = "";
  for (const model of models) {
    const option = document.createElement("option");
    option.value = model.id;
    option.textContent = model.name;
    DOM.scoringModel.appendChild(option);
  }
  if (models.some(m => m.id === current)) DOM.scoringModel.value = current;
  else if (selectedType === "protein") DOM.scoringModel.value = "blosum62";
  else DOM.scoringModel.value = "simple_nt";
  updateExplanation();
}

function updateExplanation() {
  const model = SCORING_MODELS.find(m => m.id === DOM.scoringModel.value) || SCORING_MODELS[0];
  DOM.modelExplanation.innerHTML = `
    <ul>
      <li><b>Global / Needleman-Wunsch:</b> align full-length sequences end-to-end. Use for sequences of similar length that should align across their whole span.</li>
      <li><b>Local / Smith-Waterman:</b> find the best matching region. Use when only a domain, motif, exon, or fragment is shared.</li>
      <li><b>Semi-global / overlap:</b> do not penalize terminal gaps strongly. Use for reads/fragments, partial contigs, primer products, or sequences with incomplete ends.</li>
      <li><b>${escapeHtml(model.name)}:</b> ${escapeHtml(model.explanation)}</li>
      <li><b>Pairwise:</b> best for two sequences. <b>Multiple:</b> uses a simple progressive center-star MSA; it is convenient and local, but for many divergent sequences validate with MAFFT/MUSCLE/Clustal Omega.</li>
    </ul>`;
}

function updateInputSummary() {
  try {
    const records = parseFasta(DOM.fastaInput.value);
    updateAnnotationSeqSelect(records);
    if (!records.length) {
      DOM.inputSummary.textContent = "No sequences loaded.";
      populateModels();
      return;
    }
    const type = DOM.seqType.value === "auto" ? detectType(records) : DOM.seqType.value;
    const lengths = records.map(r => r.sequence.replace(/-/g, "").length);
    DOM.inputSummary.textContent = `${records.length} sequence(s), ${type.toUpperCase()}, length range ${Math.min(...lengths)}–${Math.max(...lengths)} residues/bases.`;
    populateModels();
  } catch (err) {
    DOM.inputSummary.textContent = err.message;
    updateAnnotationSeqSelect([]);
  }
}

function getSettings(type) {
  return {
    type,
    taskMode: DOM.taskMode.value,
    algorithm: DOM.algorithm.value,
    scoringModel: DOM.scoringModel.value,
    gapOpen: Number(DOM.gapOpen.value),
    gapExtend: Number(DOM.gapExtend.value),
    match: Number(DOM.matchScore.value),
    mismatch: Number(DOM.mismatchScore.value),
    colorScheme: DOM.colorScheme.value,
    wrapCols: clamp(Number(DOM.wrapCols.value) || 80, 20, 200),
    showConsensus: DOM.showConsensus.checked,
    showSimilarity: DOM.showSimilarity.checked,
    showColumnNumbers: DOM.showColumnNumbers.checked
  };
}

function clamp(x, lo, hi) { return Math.max(lo, Math.min(hi, x)); }

function scoreChars(aRaw, bRaw, settings) {
  const a = normalizeResidue(aRaw, settings.type);
  const b = normalizeResidue(bRaw, settings.type);
  if (a === "-" || b === "-") return Number.NEGATIVE_INFINITY / 4;
  const model = settings.scoringModel;

  if (settings.type === "dna" || settings.type === "rna") {
    const aa = a === "U" ? "T" : a;
    const bb = b === "U" ? "T" : b;
    if (model === "transition_nt") {
      if (aa === bb) return settings.match;
      if ((aa === "A" && bb === "G") || (aa === "G" && bb === "A") || (aa === "C" && bb === "T") || (aa === "T" && bb === "C")) return Math.round(settings.mismatch / 2);
      return settings.mismatch * 2;
    }
    if (model === "iupac_nt") {
      const setA = IUPAC_DNA[aa] || new Set([aa]);
      const setB = IUPAC_DNA[bb] || new Set([bb]);
      if (aa === bb) return settings.match;
      if (intersects(setA, setB)) return Math.max(settings.mismatch + 1, Math.floor((settings.match + settings.mismatch) / 2));
      return settings.mismatch;
    }
    return aa === bb ? settings.match : settings.mismatch;
  }

  if (model === "blosum62") return matrixScore(BLOSUM62, a, b, settings);
  if (model === "pam250") return matrixScore(PAM250, a, b, settings);
  if (model === "protein_biochemical") {
    if (a === b) return 5;
    if (sameNamedGroup(a, b)) return 2;
    if (WEAK_PROTEIN_GROUPS.some(g => g.has(a) && g.has(b))) return 1;
    return -3;
  }
  return a === b ? 5 : -4;
}

function matrixScore(matrix, a, b, settings) {
  if (matrix[a] && matrix[a][b] !== undefined) return matrix[a][b];
  if (a === b) return settings.match;
  if (sameNamedGroup(a, b)) return 1;
  return settings.mismatch;
}

function normalizeResidue(ch, type) {
  const c = (ch || "-").toUpperCase();
  if (type === "rna" && c === "T") return "U";
  if (type === "dna" && c === "U") return "T";
  return c;
}

function intersects(a, b) {
  for (const x of a) if (b.has(x)) return true;
  return false;
}

function sameNamedGroup(a, b) {
  return Object.values(PROTEIN_GROUPS).some(g => g.has(a) && g.has(b));
}

function alignPairAffine(seqA, seqB, settings) {
  const a = seqA.replace(/-/g, "").toUpperCase();
  const b = seqB.replace(/-/g, "").toUpperCase();
  const n = a.length, m = b.length;
  const NEG = -1e12;
  const gapOpen = settings.gapOpen;
  const gapExtend = settings.gapExtend;
  const mode = settings.algorithm;

  const M = make2D(n + 1, m + 1, NEG);
  const X = make2D(n + 1, m + 1, NEG); // gap in B: consume A
  const Y = make2D(n + 1, m + 1, NEG); // gap in A: consume B
  const traceM = make2D(n + 1, m + 1, null);
  const traceX = make2D(n + 1, m + 1, null);
  const traceY = make2D(n + 1, m + 1, null);

  M[0][0] = 0;
  X[0][0] = Y[0][0] = NEG;

  if (mode === "global") {
    for (let i = 1; i <= n; i++) {
      X[i][0] = gapOpen + (i - 1) * gapExtend;
      traceX[i][0] = i === 1 ? "M" : "X";
    }
    for (let j = 1; j <= m; j++) {
      Y[0][j] = gapOpen + (j - 1) * gapExtend;
      traceY[0][j] = j === 1 ? "M" : "Y";
    }
  } else if (mode === "semiglobal") {
    for (let i = 1; i <= n; i++) { M[i][0] = 0; X[i][0] = 0; traceM[i][0] = "M"; traceX[i][0] = "X"; }
    for (let j = 1; j <= m; j++) { M[0][j] = 0; Y[0][j] = 0; traceM[0][j] = "M"; traceY[0][j] = "Y"; }
  } else { // local
    for (let i = 0; i <= n; i++) { M[i][0] = 0; X[i][0] = 0; Y[i][0] = 0; }
    for (let j = 0; j <= m; j++) { M[0][j] = 0; X[0][j] = 0; Y[0][j] = 0; }
  }

  let best = { score: 0, i: n, j: m, state: "M" };

  for (let i = 1; i <= n; i++) {
    for (let j = 1; j <= m; j++) {
      const s = scoreChars(a[i - 1], b[j - 1], settings);
      const prevForM = maxState([
        [M[i - 1][j - 1], "M"],
        [X[i - 1][j - 1], "X"],
        [Y[i - 1][j - 1], "Y"]
      ]);
      M[i][j] = prevForM[0] + s;
      traceM[i][j] = prevForM[1];

      const prevForX = maxState([
        [M[i - 1][j] + gapOpen, "M"],
        [X[i - 1][j] + gapExtend, "X"],
        [Y[i - 1][j] + gapOpen, "Y"]
      ]);
      X[i][j] = prevForX[0];
      traceX[i][j] = prevForX[1];

      const prevForY = maxState([
        [M[i][j - 1] + gapOpen, "M"],
        [Y[i][j - 1] + gapExtend, "Y"],
        [X[i][j - 1] + gapOpen, "X"]
      ]);
      Y[i][j] = prevForY[0];
      traceY[i][j] = prevForY[1];

      if (mode === "local") {
        if (M[i][j] < 0) { M[i][j] = 0; traceM[i][j] = null; }
        if (X[i][j] < 0) { X[i][j] = 0; traceX[i][j] = null; }
        if (Y[i][j] < 0) { Y[i][j] = 0; traceY[i][j] = null; }
      }

      if (mode === "local") {
        const candidate = maxState([[M[i][j], "M"], [X[i][j], "X"], [Y[i][j], "Y"]]);
        if (candidate[0] > best.score) best = { score: candidate[0], i, j, state: candidate[1] };
      }
    }
  }

  if (mode === "global") {
    const end = maxState([[M[n][m], "M"], [X[n][m], "X"], [Y[n][m], "Y"]]);
    best = { score: end[0], i: n, j: m, state: end[1] };
  } else if (mode === "semiglobal") {
    best = { score: NEG, i: n, j: m, state: "M" };
    for (let i = 0; i <= n; i++) {
      const cand = maxState([[M[i][m], "M"], [X[i][m], "X"], [Y[i][m], "Y"]]);
      if (cand[0] > best.score) best = { score: cand[0], i, j: m, state: cand[1] };
    }
    for (let j = 0; j <= m; j++) {
      const cand = maxState([[M[n][j], "M"], [X[n][j], "X"], [Y[n][j], "Y"]]);
      if (cand[0] > best.score) best = { score: cand[0], i: n, j, state: cand[1] };
    }
  }

  let i = best.i, j = best.j, state = best.state;
  const outA = [], outB = [];

  // Add unpenalized terminal overhangs for semiglobal without moving the traceback start.
  if (mode === "semiglobal") {
    const suffixA = [];
    const suffixB = [];
    for (let ti = i; ti < n; ti++) { suffixA.push(a[ti]); suffixB.push("-"); }
    for (let tj = j; tj < m; tj++) { suffixA.push("-"); suffixB.push(b[tj]); }
    tracebackCore();
    outA.reverse();
    outB.reverse();
    return finalizePair(outA.concat(suffixA).join(""), outB.concat(suffixB).join(""), best.score, settings);
  }

  tracebackCore();

  if (mode === "global") {
    while (i > 0) { outA.push(a[i - 1]); outB.push("-"); i--; }
    while (j > 0) { outA.push("-"); outB.push(b[j - 1]); j--; }
  }

  outA.reverse();
  outB.reverse();
  return finalizePair(outA.join(""), outB.join(""), best.score, settings);

  function tracebackCore() {
    while (i > 0 || j > 0) {
      if (mode === "local") {
        const val = state === "M" ? M[i][j] : state === "X" ? X[i][j] : Y[i][j];
        if (val <= 0) break;
      }
      if (mode === "semiglobal" && (i === 0 || j === 0)) break;
      if (state === "M") {
        const prev = traceM[i][j];
        outA.push(a[i - 1]); outB.push(b[j - 1]);
        i--; j--; state = prev || "M";
      } else if (state === "X") {
        const prev = traceX[i][j];
        outA.push(a[i - 1]); outB.push("-");
        i--; state = prev || "M";
      } else {
        const prev = traceY[i][j];
        outA.push("-"); outB.push(b[j - 1]);
        j--; state = prev || "M";
      }
    }
    if (mode === "semiglobal") {
      while (i > 0) { outA.push(a[i - 1]); outB.push("-"); i--; }
      while (j > 0) { outA.push("-"); outB.push(b[j - 1]); j--; }
    }
  }
}

function finalizePair(alignedA, alignedB, score, settings) {
  const stats = computePairStats(alignedA, alignedB, settings);
  return { alignedA, alignedB, score, stats };
}

function make2D(rows, cols, fill) {
  return Array.from({ length: rows }, () => Array(cols).fill(fill));
}

function maxState(items) {
  let best = items[0];
  for (const item of items) if (item[0] > best[0]) best = item;
  return best;
}

function computePairStats(a, b, settings) {
  let aligned = 0, identical = 0, similar = 0, gaps = 0;
  for (let i = 0; i < a.length; i++) {
    if (a[i] === "-" || b[i] === "-") { gaps++; continue; }
    aligned++;
    if (normalizeResidue(a[i], settings.type) === normalizeResidue(b[i], settings.type)) identical++;
    else if (areSimilar(a[i], b[i], settings)) similar++;
  }
  return {
    length: a.length,
    aligned,
    gaps,
    identical,
    similar,
    identityPct: aligned ? 100 * identical / aligned : 0,
    similarityPct: aligned ? 100 * (identical + similar) / aligned : 0
  };
}

function areSimilar(aRaw, bRaw, settings) {
  const a = normalizeResidue(aRaw, settings.type);
  const b = normalizeResidue(bRaw, settings.type);
  if (a === "-" || b === "-") return false;
  if (a === b) return true;
  if (settings.type === "dna" || settings.type === "rna") {
    const aa = a === "U" ? "T" : a;
    const bb = b === "U" ? "T" : b;
    const setA = IUPAC_DNA[aa] || new Set([aa]);
    const setB = IUPAC_DNA[bb] || new Set([bb]);
    return intersects(setA, setB);
  }
  if (settings.scoringModel === "blosum62") return matrixScore(BLOSUM62, a, b, settings) > 0;
  if (settings.scoringModel === "pam250") return matrixScore(PAM250, a, b, settings) > 0;
  return sameNamedGroup(a, b) || WEAK_PROTEIN_GROUPS.some(g => g.has(a) && g.has(b));
}

function runAlignment(records, settings) {
  const task = settings.taskMode === "auto" ? (records.length === 2 ? "pairwise" : "msa") : settings.taskMode;
  if (task === "pairwise") {
    const pair = alignPairAffine(records[0].sequence, records[1].sequence, settings);
    return {
      type: "pairwise",
      names: [records[0].name, records[1].name],
      aligned: [pair.alignedA, pair.alignedB],
      score: pair.score,
      pairStats: pair.stats,
      settings
    };
  }
  return progressiveMSA(records, settings);
}

function progressiveMSA(records, settings) {
  const centerIndex = chooseCenter(records, settings);
  const center = records[centerIndex];
  const order = records.map((_, i) => i).filter(i => i !== centerIndex)
    .sort((i, j) => records[j].sequence.length - records[i].sequence.length);

  let msa = [{ name: center.name, originalIndex: center.index, aligned: center.sequence.replace(/-/g, "").toUpperCase() }];
  let centerAligned = msa[0].aligned;
  const pairScores = [];

  for (const idx of order) {
    const rec = records[idx];
    const pair = alignPairAffine(center.sequence, rec.sequence, settings);
    pairScores.push(pair.score);
    const merged = mergeCenterAlignment(msa, centerAligned, pair.alignedA, { name: rec.name, originalIndex: rec.index, aligned: pair.alignedB });
    msa = merged.msa;
    centerAligned = merged.centerAligned;
  }

  msa.sort((a, b) => a.originalIndex - b.originalIndex);
  const aligned = msa.map(r => r.aligned);
  const names = msa.map(r => r.name);
  return {
    type: "msa",
    names,
    aligned,
    score: pairScores.reduce((a, b) => a + b, 0),
    centerName: center.name,
    centerIndex: center.index,
    settings,
    msaStats: computeMSAStats(aligned, settings)
  };
}

function chooseCenter(records, settings) {
  if (records.length <= 2) return 0;
  let bestIdx = 0, bestScore = -Infinity;
  for (let i = 0; i < records.length; i++) {
    let total = 0;
    for (let j = 0; j < records.length; j++) {
      if (i === j) continue;
      const pair = alignPairAffine(records[i].sequence, records[j].sequence, { ...settings, algorithm: "global" });
      total += pair.stats.identityPct - pair.stats.gaps * 0.05;
    }
    if (total > bestScore) { bestScore = total; bestIdx = i; }
  }
  return bestIdx;
}

function mergeCenterAlignment(currentMsa, currentCenter, pairCenter, newRecord) {
  const newRows = currentMsa.map(r => ({ ...r, aligned: "" }));
  const appended = { name: newRecord.name, originalIndex: newRecord.originalIndex, aligned: "" };
  let i = 0, j = 0;

  while (i < currentCenter.length || j < pairCenter.length) {
    const c = i < currentCenter.length ? currentCenter[i] : null;
    const p = j < pairCenter.length ? pairCenter[j] : null;

    if (c !== null && p !== null && c !== "-" && p !== "-") {
      copyExistingColumn(i, newRecord.aligned[j]);
      i++; j++;
    } else if (c === "-" && (p === null || p !== "-")) {
      copyExistingColumn(i, "-");
      i++;
    } else if ((c === null || c !== "-") && p === "-") {
      insertGapColumn(newRecord.aligned[j]);
      j++;
    } else if (c === "-" && p === "-") {
      copyExistingColumn(i, newRecord.aligned[j]);
      i++; j++;
    } else if (c !== null && p === null) {
      copyExistingColumn(i, "-");
      i++;
    } else if (c === null && p !== null) {
      insertGapColumn(newRecord.aligned[j]);
      j++;
    } else {
      break;
    }
  }

  newRows.push(appended);
  const centerAligned = newRows[0].aligned;
  return { msa: newRows, centerAligned };

  function copyExistingColumn(colIndex, newChar) {
    for (let r = 0; r < currentMsa.length; r++) newRows[r].aligned += currentMsa[r].aligned[colIndex] || "-";
    appended.aligned += newChar || "-";
  }

  function insertGapColumn(newChar) {
    for (let r = 0; r < currentMsa.length; r++) newRows[r].aligned += "-";
    appended.aligned += newChar || "-";
  }
}

function computeMSAStats(aligned, settings) {
  const len = aligned[0]?.length || 0;
  let identicalCols = 0, similarCols = 0, variableCols = 0, gapCols = 0;
  for (let col = 0; col < len; col++) {
    const chars = aligned.map(s => s[col]);
    const nonGap = chars.filter(c => c !== "-");
    if (!nonGap.length) { gapCols++; continue; }
    if (new Set(nonGap.map(c => normalizeResidue(c, settings.type))).size === 1 && nonGap.length === chars.length) identicalCols++;
    else if (columnMostlySimilar(nonGap, settings)) similarCols++;
    else variableCols++;
  }
  return { length: len, identicalCols, similarCols, variableCols, gapCols };
}

function columnMostlySimilar(chars, settings) {
  if (chars.length < 2) return false;
  let pairs = 0, sim = 0;
  for (let i = 0; i < chars.length; i++) {
    for (let j = i + 1; j < chars.length; j++) {
      pairs++;
      if (areSimilar(chars[i], chars[j], settings)) sim++;
    }
  }
  return pairs > 0 && sim / pairs >= 0.6;
}

function consensusAt(aligned, col, settings) {
  const counts = new Map();
  for (const row of aligned) {
    const c = row[col];
    if (!c || c === "-") continue;
    const n = normalizeResidue(c, settings.type);
    counts.set(n, (counts.get(n) || 0) + 1);
  }
  if (!counts.size) return "-";
  const sorted = [...counts.entries()].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
  const top = sorted[0];
  const total = [...counts.values()].reduce((a, b) => a + b, 0);
  if (top[1] / total >= 0.5) return top[0];
  return settings.type === "protein" ? "X" : "N";
}

function similarityMark(aligned, col, settings) {
  const chars = aligned.map(s => s[col]).filter(c => c !== "-");
  if (!chars.length) return " ";
  if (chars.length === aligned.length && new Set(chars.map(c => normalizeResidue(c, settings.type))).size === 1) return "*";
  return columnMostlySimilar(chars, settings) ? ":" : " ";
}

function renderResult(result) {
  currentResult = result;
  const settings = result.settings;
  const fasta = toFasta(result.names, result.aligned);
  DOM.alignedOutput.value = fasta;
  DOM.copyFastaBtn.disabled = false;
  DOM.downloadFastaBtn.disabled = false;
  DOM.downloadSvgBtn.disabled = false;
  renderStats(result);
  renderLegend(settings);
  renderAlignmentViewer(result);
  renderColorCustomizer(settings);
}

function renderStats(result) {
  if (result.type === "pairwise") {
    const s = result.pairStats;
    DOM.resultStats.textContent = `Pairwise score ${round(result.score, 2)} | aligned length ${s.length} | identity ${s.identityPct.toFixed(1)}% | similarity ${s.similarityPct.toFixed(1)}% | gaps ${s.gaps}.`;
  } else {
    const s = result.msaStats;
    DOM.resultStats.textContent = `MSA with ${result.names.length} sequences | length ${s.length} | identical columns ${s.identicalCols} | similar columns ${s.similarCols} | center sequence: ${result.centerName}.`;
  }
}

function renderAlignmentViewer(result) {
  const wrap = result.settings.wrapCols;
  const len = result.aligned[0]?.length || 0;
  DOM.viewer.innerHTML = "";
  for (let start = 0; start < len; start += wrap) {
    const end = Math.min(start + wrap, len);
    const block = document.createElement("div");
    block.className = "aln-block";
    
    // Render Position Ruler Row
    if (result.settings.showColumnNumbers) {
      const L = end - start;
      const rulerChars = Array(L).fill(" ");
      for (let col = start; col < end; col++) {
        const colNum = col + 1;
        if (colNum === 1 || colNum % 10 === 0) {
          const numStr = String(colNum);
          for (let i = 0; i < numStr.length; i++) {
            const idx = col - start + i;
            if (idx < L) rulerChars[idx] = numStr[i];
          }
        } else if (colNum % 5 === 0) {
          const idx = col - start;
          if (idx < L && rulerChars[idx] === " ") rulerChars[idx] = "·";
        }
      }
      const rulerRow = renderSpecialRow("Position", rulerChars.join(""), start, end, result, "ruler-row");
      block.appendChild(rulerRow);
    }
    
    // Render sequence rows
    for (let row = 0; row < result.aligned.length; row++) {
      block.appendChild(renderAlignmentRow(result.names[row], result.aligned[row], start, end, result, row));
    }
    
    // Render Approach 1 - Range Annotations
    if (result.rangeAnnotations && result.rangeAnnotations.length > 0) {
      result.rangeAnnotations.forEach((ann, annIndex) => {
        const annotSeq = Array(len).fill(" ");
        if (ann.ranges && ann.ranges.length > 0) {
          for (const r of ann.ranges) {
            const W = r.end - r.start + 1;
            const segmentStr = formatSimpleLabel(W, ann.label);
            for (let i = 0; i < W; i++) {
              annotSeq[r.start - 1 + i] = segmentStr[i];
            }
          }
        } else if (ann.start && ann.end) {
          const W = ann.end - ann.start + 1;
          const segmentStr = formatSimpleLabel(W, ann.label);
          for (let i = 0; i < W; i++) {
            annotSeq[ann.start - 1 + i] = segmentStr[i];
          }
        }
        const rowSeq = annotSeq.slice(start, end).join("");
        const row = renderSpecialRow(ann.name || ann.label || "", rowSeq, start, end, result, "annotation-row", ann);
        row.setAttribute("data-ann-index", annIndex);
        block.appendChild(row);
      });
    }
    
    if (result.settings.showConsensus) {
      let cons = "";
      for (let col = start; col < end; col++) cons += consensusAt(result.aligned, col, result.settings);
      const row = renderSpecialRow("Consensus", cons, start, end, result, "consensus-row");
      block.appendChild(row);
    }
    if (result.settings.showSimilarity) {
      let marks = "";
      for (let col = start; col < end; col++) marks += similarityMark(result.aligned, col, result.settings);
      const row = renderSpecialRow("Similarity", marks, start, end, result, "similarity-row");
      block.appendChild(row);
    }
    DOM.viewer.appendChild(block);
  }
  redrawSVGOverlays(result);
}

function redrawSVGOverlays(result) {
  if (!result || !DOM.viewer) return;
  const blocks = DOM.viewer.querySelectorAll(".aln-block");
  blocks.forEach(block => {
    let svg = block.querySelector(".aln-block-svg-overlay");
    if (!svg) {
      svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
      svg.setAttribute("class", "aln-block-svg-overlay");
      block.appendChild(svg);
    }
    svg.innerHTML = "";
    
    const blockRect = block.getBoundingClientRect();
    if (blockRect.width === 0 || blockRect.height === 0) return;
    
    const annRows = block.querySelectorAll(".aln-row.annotation-row");
    annRows.forEach(rowEl => {
      const annIndex = parseInt(rowEl.getAttribute("data-ann-index"), 10);
      if (isNaN(annIndex)) return;
      const ann = result.rangeAnnotations && result.rangeAnnotations[annIndex];
      if (!ann) return;
      
      const cellElements = Array.from(rowEl.querySelectorAll(".cell"));
      const segments = [];
      let currentSegment = [];
      for (const cell of cellElements) {
        if (cell.classList.contains("ann-block-cell")) {
          currentSegment.push(cell);
        } else {
          if (currentSegment.length > 0) {
            segments.push(currentSegment);
            currentSegment = [];
          }
        }
      }
      if (currentSegment.length > 0) {
        segments.push(currentSegment);
      }
      
      segments.forEach(segment => {
        const startCell = segment[0];
        const endCell = segment[segment.length - 1];
        
        const startRect = startCell.getBoundingClientRect();
        const endRect = endCell.getBoundingClientRect();
        
        const x1 = startRect.left - blockRect.left;
        const x2 = endRect.right - blockRect.left;
        const y = startRect.bottom - blockRect.top;
        
        const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
        const h = 5;
        const r = 3;
        const xMid = (x1 + x2) / 2;
        
        let d = "";
        if (x2 - x1 < 8) {
          d = `M ${x1} ${y} L ${xMid} ${y + h} L ${x2} ${y}`;
        } else {
          d = `M ${x1} ${y} 
               Q ${x1} ${y + h} ${x1 + r} ${y + h} 
               L ${xMid - r} ${y + h} 
               Q ${xMid} ${y + h} ${xMid} ${y + h + r} 
               Q ${xMid} ${y + h} ${xMid + r} ${y + h} 
               L ${x2 - r} ${y + h} 
               Q ${x2} ${y + h} ${x2} ${y}`;
        }
        
        path.setAttribute("d", d);
        path.setAttribute("fill", "none");
        path.setAttribute("stroke", ann.color || "var(--primary)");
        path.setAttribute("stroke-width", "1.5");
        svg.appendChild(path);
      });
    });
  });
}


function renderAlignmentRow(name, sequence, start, end, result, rowIndex) {
  const row = document.createElement("div");
  row.className = "aln-row";
  const nameEl = document.createElement("div");
  nameEl.className = "aln-name";
  nameEl.title = name;
  nameEl.textContent = name;
  const cells = document.createElement("div");
  cells.className = "aln-cells";
  
  const residueAnnotations = result.residueAnnotations;
  let unalignedPos = sequence.slice(0, start).replace(/-/g, "").length;
  
  for (let col = start; col < end; col++) {
    const char = sequence[col];
    const cellEl = makeCell(char, result, col);
    if (char !== "-") {
      unalignedPos++;
      if (residueAnnotations) {
        const ann = residueAnnotations.get(`${name}:${unalignedPos}`);
        if (ann) {
          cellEl.classList.add("annotated");
          cellEl.style.setProperty("--annotation-color", ann.color);
          cellEl.setAttribute("data-tooltip", `${ann.label} (${name} Pos ${unalignedPos})`);
        }
      }
    }
    cells.appendChild(cellEl);
  }
  
  const pos = document.createElement("div");
  pos.className = "aln-pos";
  pos.textContent = residuePosition(sequence.slice(0, end));
  row.append(nameEl, cells, pos);
  return row;
}

function renderSpecialRow(name, sequence, start, end, result, className, annotationObj = null) {
  const row = document.createElement("div");
  row.className = `aln-row ${className}`;
  const nameEl = document.createElement("div");
  nameEl.className = "aln-name";
  nameEl.textContent = name;
  const cells = document.createElement("div");
  cells.className = "aln-cells";
  for (let i = 0; i < sequence.length; i++) cells.appendChild(makeCell(sequence[i], result, start + i, className, annotationObj));
  const pos = document.createElement("div");
  pos.className = "aln-pos";
  pos.textContent = "";
  row.append(nameEl, cells, pos);
  return row;
}

function makeCell(char, result, col, special = "", annotationObj = null) {
  const span = document.createElement("span");
  const c = char || " ";
  if (special.includes("annotation") || special.includes("ruler")) {
    span.textContent = c === " " ? "\u00A0" : c;
  } else {
    span.textContent = c === " " ? "·" : c;
  }
  
  span.className = `cell ${cellClass(c, result.settings, result.aligned, col)} ${conservationClass(result.aligned, col, c, result.settings)}`;
  
  if (special.includes("similarity") || special.includes("ruler")) span.className = "cell mono";
  
  if (special.includes("annotation") && annotationObj) {
    span.className = "cell";
    let inRange = false;
    if (annotationObj.ranges && annotationObj.ranges.length > 0) {
      for (const r of annotationObj.ranges) {
        if (col >= r.start - 1 && col <= r.end - 1) {
          inRange = true;
          break;
        }
      }
    } else if (annotationObj.start && annotationObj.end) {
      if (col >= annotationObj.start - 1 && col <= annotationObj.end - 1) {
        inRange = true;
      }
    }
    
    if (inRange) {
      span.classList.add("ann-block-cell");
      span.style.backgroundColor = annotationObj.color;
      span.style.color = annotationObj.textColor || "#ffffff";
    }
  }
  
  return span;
}

function residuePosition(prefix) {
  return String(prefix.replace(/-/g, "").length || 0);
}

function conservationClass(aligned, col, char, settings) {
  if (char === "-" || char === " ") return "";
  const mark = similarityMark(aligned, col, settings);
  if (mark === "*") return "identical";
  if (mark === ":") return "similar";
  return "variable";
}

function cellClass(charRaw, settings, aligned, col) {
  const char = normalizeResidue(charRaw, settings.type);
  if (char === "-") return "gap";
  let scheme = settings.colorScheme;
  if (scheme === "auto") scheme = settings.type === "protein" ? "clustal" : "nucleotide";
  if (scheme === "mono") return "mono";
  if (scheme === "identity") {
    const mark = similarityMark(aligned, col, settings);
    if (mark === "*") return "conserved";
    if (mark === ":") return "simconserved";
    return "mono";
  }
  if (scheme === "nucleotide" || settings.type !== "protein") return `nt-${char}`;
  if (scheme === "zappo") return proteinZappoClass(char);
  if (scheme === "taylor") return proteinTaylorClass(char);
  return proteinClustalClass(char);
}

function proteinClustalClass(c) {
  if ("AVLIM".includes(c)) return "pro-hydrophobic";
  if ("FWY".includes(c)) return "pro-aromatic";
  if ("KRH".includes(c)) return "pro-positive";
  if ("DE".includes(c)) return "pro-negative";
  if ("STNQ".includes(c)) return "pro-polar";
  if ("GP".includes(c)) return "pro-special";
  if ("C".includes(c)) return "pro-sulfur";
  return "mono";
}

function proteinZappoClass(c) {
  if ("AILMFPWV".includes(c)) return "pro-hydrophobic";
  if ("KRH".includes(c)) return "pro-positive";
  if ("DE".includes(c)) return "pro-negative";
  if ("STNQ".includes(c)) return "pro-polar";
  if ("CY".includes(c)) return "pro-sulfur";
  if ("G".includes(c)) return "pro-small";
  return "mono";
}

function proteinTaylorClass(c) {
  if ("VILM".includes(c)) return "pro-hydrophobic";
  if ("FYWH".includes(c)) return "pro-aromatic";
  if ("KR".includes(c)) return "pro-positive";
  if ("DE".includes(c)) return "pro-negative";
  if ("STNQ".includes(c)) return "pro-polar";
  if ("AGP".includes(c)) return "pro-small";
  if ("C".includes(c)) return "pro-sulfur";
  return "mono";
}

function renderLegend(settings) {
  let scheme = settings.colorScheme === "auto" ? (settings.type === "protein" ? "clustal" : "nucleotide") : settings.colorScheme;
  const items = [];
  if (scheme === "nucleotide" || settings.type !== "protein") {
    items.push(["A", "nt-A"], ["C", "nt-C"], ["G", "nt-G"], [settings.type === "rna" ? "U" : "T", settings.type === "rna" ? "nt-U" : "nt-T"], ["Ambiguous", "nt-N"], ["Gap", "gap"]);
  } else if (scheme === "identity") {
    items.push(["Identical column", "conserved"], ["Similar column", "simconserved"], ["Variable", "mono"], ["Gap", "gap"]);
  } else if (scheme === "mono") {
    items.push(["Residue/base", "mono"], ["Gap", "gap"]);
  } else {
    items.push(["Hydrophobic", "pro-hydrophobic"], ["Aromatic", "pro-aromatic"], ["Positive", "pro-positive"], ["Negative", "pro-negative"], ["Polar", "pro-polar"], ["Small/special", "pro-small"], ["Sulfur", "pro-sulfur"], ["Gap", "gap"]);
  }
  items.push(["Outline: identical", "identical"], ["Outline: similar", "similar"]);
  DOM.legend.innerHTML = items.map(([label, cls]) => `<span class="legend-item"><span class="legend-swatch cell ${cls}"></span>${escapeHtml(label)}</span>`).join("");
}

function toFasta(names, aligned) {
  return names.map((name, i) => `>${name}\n${wrapString(aligned[i], 80)}`).join("\n");
}

function wrapString(text, width) {
  const out = [];
  for (let i = 0; i < text.length; i += width) out.push(text.slice(i, i + width));
  return out.join("\n");
}

function downloadText(filename, text, mime = "text/plain") {
  const blob = new Blob([text], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

function makeAlignmentSVG(result) {
  const settings = result.settings;
  const cellW = 14, cellH = 18, nameW = 190, top = 36, left = 12;
  const stepW = cellW + currentGap;
  const stepH = cellH + currentGap;
  const len = result.aligned[0]?.length || 0;
  const annCount = result.rangeAnnotations ? result.rangeAnnotations.length : 0;
  const hasRuler = settings.showColumnNumbers ? 1 : 0;
  const rows = result.aligned.length + annCount + (settings.showConsensus ? 1 : 0) + (settings.showSimilarity ? 1 : 0) + hasRuler;
  const width = Math.max(900, left + nameW + len * stepW + 80);
  const height = top + rows * stepH + 70;
  const esc = escapeHtml;
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`;
  svg += `<rect width="100%" height="100%" fill="#ffffff"/><text x="${left}" y="24" font-family="Arial" font-size="18" font-weight="700">Sequence alignment</text>`;
  let y = top;
  
  // Render Position Ruler Row
  if (settings.showColumnNumbers) {
    svg += `<text x="${left}" y="${y + 13}" font-family="Arial" font-size="11" font-weight="700" fill="#64748b">Position</text>`;
    const rulerChars = Array(len).fill(" ");
    for (let col = 0; col < len; col++) {
      const colNum = col + 1;
      if (colNum === 1 || colNum % 10 === 0) {
        const numStr = String(colNum);
        for (let i = 0; i < numStr.length; i++) {
          const idx = col + i;
          if (idx < len) rulerChars[idx] = numStr[i];
        }
      } else if (colNum % 5 === 0) {
        if (rulerChars[col] === " ") rulerChars[col] = "·";
      }
    }
    for (let c = 0; c < len; c++) {
      svg += svgRulerCell(rulerChars[c], left + nameW + c * stepW, y, cellW, cellH);
    }
    y += stepH;
  }
  
  // Render sequence rows with Approach 2 residue annotations
  for (let r = 0; r < result.aligned.length; r++) {
    const name = result.names[r];
    const sequence = result.aligned[r];
    svg += `<text x="${left}" y="${y + 13}" font-family="monospace" font-size="12" fill="#334155">${esc(name.slice(0, 26))}</text>`;
    let unalignedPos = 0;
    for (let c = 0; c < len; c++) {
      const char = sequence[c];
      if (char !== "-") unalignedPos++;
      const ann = result.residueAnnotations && result.residueAnnotations.get(`${name}:${unalignedPos}`);
      svg += svgCell(char, left + nameW + c * stepW, y, cellW, cellH, result, c, false, false, ann);
    }
    y += stepH;
  }
  
  // Render Approach 1 range annotations
  if (result.rangeAnnotations && result.rangeAnnotations.length > 0) {
    for (const ann of result.rangeAnnotations) {
      svg += `<text x="${left}" y="${y + 13}" font-family="Arial" font-size="11" font-weight="700" fill="#64748b">${esc((ann.name || ann.label || "").slice(0, 26))}</text>`;
      
      const annotSeq = Array(len).fill(" ");
      if (ann.ranges && ann.ranges.length > 0) {
        for (const r of ann.ranges) {
          const W = r.end - r.start + 1;
          const segmentStr = formatSimpleLabel(W, ann.label);
          for (let i = 0; i < W; i++) {
            annotSeq[r.start - 1 + i] = segmentStr[i];
          }
        }
      } else if (ann.start && ann.end) {
        const W = ann.end - ann.start + 1;
        const segmentStr = formatSimpleLabel(W, ann.label);
        for (let i = 0; i < W; i++) {
          annotSeq[ann.start - 1 + i] = segmentStr[i];
        }
      }

      for (let c = 0; c < len; c++) {
        let isAnnotated = false;
        if (ann.ranges && ann.ranges.length > 0) {
          for (const r of ann.ranges) {
            if (c >= r.start - 1 && c <= r.end - 1) {
              isAnnotated = true;
              break;
            }
          }
        } else if (ann.start && ann.end) {
          if (c >= ann.start - 1 && c <= ann.end - 1) {
            isAnnotated = true;
          }
        }
        svg += svgAnnotationCell(annotSeq[c], left + nameW + c * stepW, y, cellW, cellH, isAnnotated, ann.color, ann.textColor || "#ffffff");
      }
      y += stepH;
    }
  }
  
  if (settings.showConsensus) {
    svg += `<text x="${left}" y="${y + 13}" font-family="monospace" font-size="12" font-weight="700" fill="#0f766e">Consensus</text>`;
    for (let c = 0; c < len; c++) svg += svgCell(consensusAt(result.aligned, c, settings), left + nameW + c * stepW, y, cellW, cellH, result, c, true);
    y += stepH;
  }
  if (settings.showSimilarity) {
    svg += `<text x="${left}" y="${y + 13}" font-family="monospace" font-size="12" fill="#64748b">Similarity</text>`;
    for (let c = 0; c < len; c++) svg += svgCell(similarityMark(result.aligned, c, settings), left + nameW + c * stepW, y, cellW, cellH, result, c, false, true);
    y += stepH;
  }
  svg += `<text x="${left}" y="${height - 18}" font-family="Arial" font-size="12" fill="#64748b">Generated locally by Sequence Alignment Studio. * identical column, : similar column.</text>`;
  svg += `</svg>`;
  return svg;
}

function svgRulerCell(ch, x, y, w, h) {
  const label = ch === " " ? "" : ch;
  return `<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="none" stroke="none"/>` +
    `<text x="${x + w / 2}" y="${y + 13}" text-anchor="middle" font-family="monospace" font-size="10" font-weight="600" fill="#64748b">${escapeHtml(label)}</text>`;
}

function svgAnnotationCell(ch, x, y, w, h, isAnnotated, color, textColor) {
  const bg = isAnnotated ? color : "#ffffff";
  const fg = isAnnotated ? textColor : "#ffffff00";
  const stroke = isAnnotated ? color : "transparent";
  const label = ch === " " ? "" : ch;
  const sw = isAnnotated ? 0.5 : 0;
  return `<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${currentRoundness}" fill="${bg}" stroke="${stroke}" stroke-width="${sw}"/>` +
    `<text x="${x + w / 2}" y="${y + 13}" text-anchor="middle" font-family="monospace" font-size="11" font-weight="700" fill="${fg}">${escapeHtml(label)}</text>`;
}

function svgCell(ch, x, y, w, h, result, col, consensus = false, similarity = false, residueAnn = null) {
  const colors = svgColors(ch, result, col, consensus, similarity);
  const label = ch === " " ? "·" : ch;
  let sw = colors.sw;
  let stroke = colors.stroke;
  if (!consensus && !similarity && currentGap === 0 && sw === 0.5) {
    sw = 0;
  }
  
  if (residueAnn) {
    stroke = residueAnn.color;
    sw = 1.5;
  }
  
  return `<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${currentRoundness}" fill="${colors.bg}" stroke="${stroke}" stroke-width="${sw}"/>` +
    `<text x="${x + w / 2}" y="${y + 13}" text-anchor="middle" font-family="monospace" font-size="11" font-weight="700" fill="${colors.fg}">${escapeHtml(label)}</text>`;
}

function svgColors(ch, result, col, consensus, similarity) {
  if (similarity) return { bg: "#ffffff", fg: "#475569", stroke: "#ffffff", sw: 0 };
  if (consensus) {
    const bg = currentColors["consensus-bg"] || "#ccfbf1";
    const fg = currentColors["consensus-fg"] || "#0f766e";
    const isTrans = bg === "transparent";
    return { bg: isTrans ? "#ffffff" : bg, fg, stroke: isTrans ? "transparent" : fg, sw: isTrans ? 0 : 1 };
  }
  const cls = cellClass(ch, result.settings, result.aligned, col);
  const bg = currentColors[`${cls}-bg`] || currentColors["mono-bg"];
  const fg = currentColors[`${cls}-fg`] || currentColors["mono-fg"];
  const cons = conservationClass(result.aligned, col, ch, result.settings);
  if (cons === "identical") {
    const stroke = currentColors["identical-outline-border"] || "#16a34a";
    return { bg, fg, stroke, sw: stroke === "transparent" ? 0 : 1.5 };
  }
  if (cons === "similar") {
    const stroke = currentColors["similar-outline-border"] || "#ca8a04";
    return { bg, fg, stroke, sw: stroke === "transparent" ? 0 : 1.2 };
  }
  return { bg, fg, stroke: "#ffffff", sw: 0.5 };
}

function escapeHtml(text) {
  return String(text).replace(/[&<>"]/g, ch => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[ch]));
}

function round(x, digits = 2) {
  const p = 10 ** digits;
  return Math.round(x * p) / p;
}

function showError(message) {
  DOM.viewer.innerHTML = `<div class="error-box">${escapeHtml(message)}</div>`;
  DOM.resultStats.textContent = "Alignment failed.";
  DOM.copyFastaBtn.disabled = true;
  DOM.downloadFastaBtn.disabled = true;
  DOM.downloadSvgBtn.disabled = true;
}

function runFromUI() {
  try {
    const records = parseFasta(DOM.fastaInput.value);
    validateSequences(records);
    const type = DOM.seqType.value === "auto" ? detectType(records) : DOM.seqType.value;
    const settings = getSettings(type);
    if (settings.gapOpen > 0 || settings.gapExtend > 0) {
      throw new Error("Gap open and gap extend should usually be negative penalties, for example -10 and -1.");
    }
    if (records.length > 60) throw new Error("This local demo app is intended for small-to-medium alignments. Please use fewer than 60 sequences or a command-line MSA tool.");
    const maxLen = Math.max(...records.map(r => r.sequence.replace(/-/g, "").length));
    if (maxLen > 2500) throw new Error("For browser performance, keep individual sequences below ~2500 residues/bases. Use MAFFT/MUSCLE/EMBOSS for longer alignments.");
    const result = runAlignment(records, settings);
    
    // Retrieve and attach annotations
    const { rangeAnnotations, residueAnnotations } = getAnnotationsData();
    result.rangeAnnotations = rangeAnnotations;
    result.residueAnnotations = residueAnnotations;
    
    renderResult(result);
  } catch (err) {
    showError(err.message || String(err));
  }
}

DOM.fileInput.addEventListener("change", async event => {
  const file = event.target.files[0];
  if (!file) return;
  DOM.fastaInput.value = await file.text();
  updateInputSummary();
});
DOM.clearInput.addEventListener("click", () => {
  DOM.fastaInput.value = "";
  DOM.alignedOutput.value = "";
  DOM.viewer.innerHTML = "";
  DOM.legend.innerHTML = "";
  DOM.resultStats.textContent = "Run an alignment to see results.";
  currentResult = null;
  annotationsArray = [];
  onAnnotationsUpdated();
  renderColorCustomizer(null);
  updateInputSummary();
});
DOM.fastaInput.addEventListener("input", updateInputSummary);
DOM.seqType.addEventListener("change", updateInputSummary);
DOM.scoringModel.addEventListener("change", updateExplanation);
DOM.algorithm.addEventListener("change", updateExplanation);
DOM.colorScheme.addEventListener("change", () => { 
  if (currentResult) { 
    currentResult.settings.colorScheme = DOM.colorScheme.value; 
    renderLegend(currentResult.settings); 
    renderAlignmentViewer(currentResult); 
    renderColorCustomizer(currentResult.settings);
  } 
});
DOM.showConsensus.addEventListener("change", () => { if (currentResult) { currentResult.settings.showConsensus = DOM.showConsensus.checked; renderAlignmentViewer(currentResult); } });
DOM.showSimilarity.addEventListener("change", () => { if (currentResult) { currentResult.settings.showSimilarity = DOM.showSimilarity.checked; renderAlignmentViewer(currentResult); } });
DOM.wrapCols.addEventListener("change", () => { if (currentResult) { currentResult.settings.wrapCols = clamp(Number(DOM.wrapCols.value) || 80, 20, 200); renderAlignmentViewer(currentResult); } });
DOM.resetColorsBtn.addEventListener("click", () => {
  currentColors = { ...DEFAULT_COLORS };
  currentRoundness = 0;
  currentGap = 1;
  applyColors();
  applyRoundness();
  applyGap();
  if (currentResult) {
    renderColorCustomizer(currentResult.settings);
  }
});
DOM.boxRoundness.addEventListener("input", e => {
  currentRoundness = Number(e.target.value);
  applyRoundness();
});
DOM.boxGap.addEventListener("input", e => {
  currentGap = Number(e.target.value);
  applyGap();
});
DOM.makeWhiteBtn.addEventListener("click", () => {
  for (const key of Object.keys(currentColors)) {
    if (key.endsWith("-bg")) {
      currentColors[key] = "#ffffff";
    } else if (key.endsWith("-fg")) {
      currentColors[key] = "#0f172a";
    }
  }
  currentColors["consensus-bg"] = "transparent";
  currentColors["consensus-fg"] = "#0f172a";
  currentColors["identical-outline-border"] = "transparent";
  currentColors["similar-outline-border"] = "transparent";
  applyColors();
  if (currentResult) {
    renderColorCustomizer(currentResult.settings);
  }
});
DOM.runBtn.addEventListener("click", runFromUI);
DOM.copyFastaBtn.addEventListener("click", async () => {
  await navigator.clipboard.writeText(DOM.alignedOutput.value);
  DOM.copyFastaBtn.textContent = "Copied";
  setTimeout(() => DOM.copyFastaBtn.textContent = "Copy aligned FASTA", 1200);
});
DOM.downloadFastaBtn.addEventListener("click", () => downloadText("aligned_sequences.fasta", DOM.alignedOutput.value, "text/x-fasta"));
DOM.downloadSvgBtn.addEventListener("click", () => {
  if (!currentResult) return;
  downloadText("alignment_plot.svg", makeAlignmentSVG(currentResult), "image/svg+xml");
});

// Annotations Form Event Listeners
DOM.annTypeSelect.addEventListener("change", () => {
  if (DOM.annTypeSelect.value === "range") {
    DOM.rangeInputsRow.style.display = "flex";
    DOM.residueInputsRow.style.display = "none";
    DOM.annTextColorLabel.style.display = "flex";
    DOM.annTrackLabelLabel.style.display = "flex";
  } else {
    DOM.rangeInputsRow.style.display = "none";
    DOM.residueInputsRow.style.display = "flex";
    DOM.annTextColorLabel.style.display = "none";
    DOM.annTrackLabelLabel.style.display = "none";
  }
});

DOM.annHasSeg2.addEventListener("change", () => {
  if (DOM.annHasSeg2.checked) {
    DOM.seg2Row.style.display = "flex";
  } else {
    DOM.seg2Row.style.display = "none";
    DOM.annHasSeg3.checked = false;
    DOM.seg3Row.style.display = "none";
    DOM.annHasSeg4.checked = false;
    DOM.seg4Row.style.display = "none";
    DOM.annStartCol2.value = "";
    DOM.annEndCol2.value = "";
    DOM.annStartCol3.value = "";
    DOM.annEndCol3.value = "";
    DOM.annStartCol4.value = "";
    DOM.annEndCol4.value = "";
  }
});

DOM.annHasSeg3.addEventListener("change", () => {
  if (DOM.annHasSeg3.checked) {
    DOM.seg3Row.style.display = "flex";
  } else {
    DOM.seg3Row.style.display = "none";
    DOM.annHasSeg4.checked = false;
    DOM.seg4Row.style.display = "none";
    DOM.annStartCol3.value = "";
    DOM.annEndCol3.value = "";
    DOM.annStartCol4.value = "";
    DOM.annEndCol4.value = "";
  }
});

DOM.annHasSeg4.addEventListener("change", () => {
  if (DOM.annHasSeg4.checked) {
    DOM.seg4Row.style.display = "flex";
  } else {
    DOM.seg4Row.style.display = "none";
    DOM.annStartCol4.value = "";
    DOM.annEndCol4.value = "";
  }
});

DOM.cancelAnnBtn.addEventListener("click", cancelEditMode);

DOM.addAnnBtn.addEventListener("click", () => {
  const name = DOM.annNameInput.value.trim();
  const label = DOM.annLabelInput.value.trim();
  const color = DOM.annColorInput.value;
  const type = DOM.annTypeSelect.value;
  
  if (!name) {
    alert("Please enter a Feature Name.");
    return;
  }
  
  let newAnn = { type, name, label, color };
  
  if (type === "range") {
    const start = parseInt(DOM.annStartCol.value, 10);
    const end = parseInt(DOM.annEndCol.value, 10);
    if (isNaN(start) || start <= 0 || isNaN(end) || end <= 0) {
      alert("Please enter valid positive numbers for both start and end columns of Range 1.");
      return;
    }
    if (start > end) {
      alert("Start column 1 must be less than or equal to end column 1.");
      return;
    }
    
    const ranges = [{ start, end }];
    
    if (DOM.annHasSeg2.checked) {
      const start2 = parseInt(DOM.annStartCol2.value, 10);
      const end2 = parseInt(DOM.annEndCol2.value, 10);
      if (isNaN(start2) || start2 <= 0 || isNaN(end2) || end2 <= 0) {
        alert("Please enter valid positive numbers for both start and end columns of Range 2.");
        return;
      }
      if (start2 > end2) {
        alert("Start column 2 must be less than or equal to end column 2.");
        return;
      }
      ranges.push({ start: start2, end: end2 });
      
      if (DOM.annHasSeg3.checked) {
        const start3 = parseInt(DOM.annStartCol3.value, 10);
        const end3 = parseInt(DOM.annEndCol3.value, 10);
        if (isNaN(start3) || start3 <= 0 || isNaN(end3) || end3 <= 0) {
          alert("Please enter valid positive numbers for both start and end columns of Range 3.");
          return;
        }
        if (start3 > end3) {
          alert("Start column 3 must be less than or equal to end column 3.");
          return;
        }
        ranges.push({ start: start3, end: end3 });
        
        if (DOM.annHasSeg4.checked) {
          const start4 = parseInt(DOM.annStartCol4.value, 10);
          const end4 = parseInt(DOM.annEndCol4.value, 10);
          if (isNaN(start4) || start4 <= 0 || isNaN(end4) || end4 <= 0) {
            alert("Please enter valid positive numbers for both start and end columns of Range 4.");
            return;
          }
          if (start4 > end4) {
            alert("Start column 4 must be less than or equal to end column 4.");
            return;
          }
          ranges.push({ start: start4, end: end4 });
        }
      }
    }
    
    newAnn.ranges = ranges;
    newAnn.textColor = DOM.annTextColorInput.value;
    newAnn.start = start;
    newAnn.end = end;
  } else {
    const sequence = DOM.annSeqSelect.value;
    const positionInput = DOM.annResPos.value.trim();
    if (!sequence) {
      alert("Please select a sequence. If no sequences are available, load FASTA first.");
      return;
    }
    if (!positionInput) {
      alert("Please enter unaligned position(s) (e.g. 24 or 24-26, 30).");
      return;
    }
    
    const positions = parseResiduePositions(positionInput);
    if (positions.length === 0) {
      alert("Please enter valid positive unaligned position(s) (e.g. 24-26, 30).");
      return;
    }
    
    newAnn.sequence = sequence;
    newAnn.positions = positions;
    newAnn.positionInput = positionInput;
    newAnn.position = positions[0];
  }
  
  if (editingIndex !== null) {
    annotationsArray[editingIndex] = newAnn;
    cancelEditMode();
  } else {
    annotationsArray.push(newAnn);
    cancelEditMode();
  }
  
  onAnnotationsUpdated();
});

DOM.showColumnNumbers.addEventListener("change", () => {
  if (currentResult) {
    currentResult.settings.showColumnNumbers = DOM.showColumnNumbers.checked;
    renderAlignmentViewer(currentResult);
  }
});

DOM.loadDemoDna.addEventListener("click", () => { 
  DOM.fastaInput.value = DEMO_DNA; 
  annotationsArray = [];
  DOM.seqType.value = "auto"; 
  updateInputSummary(); 
  onAnnotationsUpdated();
});
DOM.loadDemoProtein.addEventListener("click", () => { 
  DOM.fastaInput.value = DEMO_PROTEIN; 
  annotationsArray = [];
  DOM.seqType.value = "auto"; 
  updateInputSummary(); 
  onAnnotationsUpdated();
});

applyColors();
applyRoundness();
applyGap();
renderColorCustomizer(null);
updateInputSummary();
onAnnotationsUpdated();

window.addEventListener("resize", () => {
  if (currentResult) {
    redrawSVGOverlays(currentResult);
  }
});
