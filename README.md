# cysloc
The Cys-LoC and CysLOx methods, enable quantitative subcellular cysteine chemoproteomics, including identification of mitochondrial cysteines sensitive to LPS+IFNg-induced oxidative stress.

Available at: https://www.sciencedirect.com/science/article/pii/S2451945623001861

The CysLoc Database and mass spectrometry processing were created and developed by Lisa Boatner.

About: Proteinaceous cysteines function as essential sensors of cellular redox state. Consequently, defining the
cysteine redoxome is a key challenge for functional proteomic studies. While proteome-wide inventories
of cysteine oxidation state are readily achieved using established, widely adopted proteomic methods
such as OxICAT, Biotin Switch, and SP3-Rox, these methods typically assay bulk proteomes and therefore
fail to capture protein localization-dependent oxidative modifications. Here we establish the local cysteine
capture (Cys-LoC) and local cysteine oxidation (Cys-LOx) methods, which together yield compartment-specific cysteine capture and quantitation of cysteine oxidation state. Benchmarking of the Cys-LoC method
across a panel of subcellular compartments revealed more than 3,500 cysteines not previously captured
by whole-cell proteomic analysis. Application of the Cys-LOx method to LPS-stimulated immortalized murine
bone marrow-derived macrophages (iBMDM), revealed previously unidentified, mitochondrially localized
cysteine oxidative modifications upon pro-inflammatory activation, including those associated with oxidative
mitochondrial metabolism.

Yan, T., Julio, A.R., Villanueva, M., Jones, A.E., Ball, A.B., Boatner, L.M., Turmon, A.C., Nguyễn, K.B., Yen, S.L., Desai, H.S. and Divakaruni, A.S., 2023. Proximity-labeling chemoproteomics defines the subcellular cysteinome and inflammation-responsive mitochondrial redoxome. Cell Chemical Biology, 30(7), 811-827.

To run (specifically for post-processing FragPipe MS label free quantitation (LFQ) data with a mouse fasta file): python3 230119_process.py -exp 'lfq' -lpm '463.2366;467.2529' -hpm '463.2366;467.2529' -proref '2022-10-04-decoys-contam-uniprot_mouse_swissprot_1042022.fasta.fas'
