# 100 Days of Computational Chemistry

One small, working, documented tool per day. Chemistry students and researchers can reuse them; each one comes from real lab and research workflows (N-CQDs, DFT/xTB, docking/MD, chromene synthesis, spectroscopy).

## Daily rules
- Folder per day: `dayNN-name/` with `README.md` (what, why, how to run), a script or notebook, and a tiny example input and output.
- Minimum viable commit on busy days: a test, a docstring, or a README improvement.
- Commit by evening IST (GitHub counts days in UTC). Never backdate.
- Every Sunday: tag a release, update the index table, post a short summary (LinkedIn/X/Instagram).
- Stack: Python, RDKit, xTB, ASE, MDAnalysis, pandas, matplotlib, Streamlit.

---

## Phase 1 (Days 1-10): Foundations and template
1. Repo setup: license, README, `requirements.txt`, folder template
2. SMILES to 3D geometry (RDKit)
3. Molecular descriptors table from a SMILES list
4. Batch xTB optimization wrapper
5. xTB output parser: energies, gaps, charges
6. XYZ file toolkit: convert, merge, split, center
7. Unit converter: Hartree, eV, kcal/mol, kJ/mol, nm, cm⁻¹
8. Periodic-table property lookup CLI
9. Molecule grid image generator for figures
10. Index page with day links and a GitHub Actions lint check

## Phase 2 (Days 11-20): Nanomaterial and spectroscopy data (N-CQD core)
11. UV-vis plotter with peak detection
12. Tauc plot and optical band gap calculator
13. Fluorescence quantum yield (relative method) calculator
14. PL spectrum plotter with excitation-dependent emission overlay
15. Stokes shift calculator
16. FTIR peak picker with functional-group suggestions
17. XRD Scherrer crystallite size calculator
18. TEM particle size distribution from measurements (histogram plus lognormal fit)
19. Baseline correction and smoothing toolkit for spectra
20. **Mini release v0.1: `spectra-kit`** combining Days 11-19 as a package

## Phase 3 (Days 21-30): DFT and xTB workflows
21. ORCA input generator
22. Psi4 input generator
23. Gaussian-style input to ORCA converter
24. HOMO-LUMO and gap extractor with orbital energy diagram
25. Conformer search with CREST: setup and result ranking
26. Boltzmann-weighted energy averaging
27. Simulated UV-vis from TDDFT output (Gaussian broadening)
28. Simulated IR spectrum from frequency output
29. Doped carbon dot model builder (N-doped graphene fragments)
30. **Mini release v0.2: `qc-workflows`**, a documented DFT/xTB pipeline

## Phase 4 (Days 31-40): Docking and MD
31. Ligand preparation pipeline (protonation, 3D, minimization)
32. Protein cleanup: remove waters, select chain, fix missing atoms
33. AutoDock Vina batch docking script
34. Docking results ranker and table exporter
35. Protein-ligand contact summary
36. Binding pose visualization snapshots
37. MDAnalysis RMSD and RMSF plots
38. Hydrogen-bond occupancy analysis
39. Free-energy result parser and summary table
40. **Mini release v0.3: `dock-md-lite`**, a small docking-to-MD pipeline

## Phase 5 (Days 41-50): Synthesis and green chemistry (chromene work)
41. Reaction yield tracker and statistics notebook
42. Atom economy calculator
43. E-factor and green metrics calculator
44. Reaction condition screening: catalyst loading, time, temperature heatmaps
45. Recyclability plot for catalysts (yield vs cycle)
46. Melting point and TLC Rf logger
47. NMR peak table to formatted experimental-section text
48. Reaction mechanism drawer from SMILES
49. Multicomponent reaction library: scope table generator
50. **Mini release v0.4: `synth-lab-notebook`**

## Phase 6 (Days 51-60): ML for chemistry
51. Molecular fingerprints and similarity search
52. Dataset cleaning for chemical data (duplicates, salts, invalid SMILES)
53. Solubility prediction baseline
54. Toxicity or activity QSAR baseline
55. Train/test split by scaffold
56. Feature importance for descriptors
57. Simple GNN for property prediction
58. Model comparison report generator
59. Applicability domain check
60. **Mini release v0.5: `chem-ml-starter`**

## Phase 7 (Days 61-70): Tools for students (largest audience)
61. CSIR-NET periodic trends quiz
62. Named reactions flashcards with mechanisms
63. Spectroscopy problem generator (IR/NMR/MS)
64. Organic mechanism arrow-pushing practice cards
65. Quantum chemistry formula cheat sheet with code examples
66. Thermodynamics and kinetics solver
67. Point group and symmetry identifier
68. Stoichiometry and molarity calculators
69. Lab report figure and table templates
70. **Mini release v0.6: `netjrf-chem-toolkit`**

## Phase 8 (Days 71-80): Apps and visualization
71. Streamlit app: Tauc plot calculator
72. Streamlit app: SMILES to properties and 3D view
73. Streamlit app: quantum yield calculator
74. Streamlit app: spectra viewer and overlay
75. Streamlit app: docking result browser
76. Publication-quality matplotlib style for chemistry figures
77. Color-blind-safe palettes and plotting presets
78. Graphical abstract template helper
79. Dockerfile for the whole toolkit
80. **Mini release v0.7: `chem-apps`** with live demo links

## Phase 9 (Days 81-90): Quality and reproducibility
81. Unit tests for core functions
82. GitHub Actions CI
83. Type hints and docstrings pass
84. Sample datasets with citations and licenses
85. Packaging for pip (`pyproject.toml`)
86. Documentation site (MkDocs)
87. Benchmark: your tools vs literature values
88. Zenodo DOI for the main release
89. `CITATION.cff` file
90. **Release v1.0** of the unified toolkit

## Phase 10 (Days 91-100): Capstone
91. Full worked example: N-CQD from raw data to band gap to DFT comparison
92. Full worked example: ligand to docking to MD to analysis
93. Tutorial notebook for beginners
94. Video or GIF demo of the best tool
95. Contribution guide and good-first-issue labels
96. Profile README update: pinned repos, preprint link, skills
97. Blog post: what I learned in 100 days
98. Collect and fix feedback from issues
99. Roadmap for the next 100 days
100. Final summary, stats, and thank-you post

---

## Progress
| Day | Topic | Link | Done |
|-----|-------|------|------|
| 1 | Repo setup | | ☐ |
