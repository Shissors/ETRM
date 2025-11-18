# Project Overview
## Main.py

A comprehensive Python script designed to automate the process of generating realistic evolutionary protein snapshots and formatting them for high-throughput structural prediction using the AlphaFold Server.

The workflow takes an ancestral protein sequence (or a set of related sequences) and simulates evolution over a configurable number of generations ("snapshots"), incorporating biological accuracy features like conservation and selection pressure. It finally outputs a single JSON file ready for batch submission to the AlphaFold Server.

### Features

Sequence Source Flexibility: Accepts a multi-sequence FASTA file and uses the first sequence as the evolutionary ancestor.

Optional Multiple Sequence Alignment (MSA): Can run MAFFT on the input file to align sequences and determine the true ancestor (requires MAFFT installation).

Biologically Accurate Evolution: Simulation includes three layers of realism:

Conservation Weighting: Sites are less likely to mutate if they are in the protein core.

Weighted Substitution: Mutations favor chemically similar amino acids (approximating a Substitution Matrix).

Lethal Selection: Destructive mutations (e.g., losing a Cysteine or changing core hydrophobicity) are highly likely to be rejected.

Batch Job Generation: Creates individual FASTA snapshots (protein_snap_XX.fasta) and compiles all sequences into a single alphafold_batch_submission.json file.

### Prerequisites

To run this workflow, you need the following installed:

Python 3.6+

MAFFT (Optional but Recommended): Required if you use the --run_msa flag. Must be callable from your system's command line (PATH).

Required Python Libraries: Biopython, pandas, argparse, json, subprocess.

and ChimeraX

### Command

RUN THIS COMMAND FIRST (in the directory of the .py file using cmd)
```
python main.py sequences.fasta --mafft_path "<path to MAFFT>" --num_snapshots 10 --mutation_rate 0.1 --lethal_rejection 0.9 --snapshot_dir my_snapshots
```

## Main2.py 

This script automates the creation of a structural morph animation video from a set of predicted protein structures (e.g., outputs from an AlphaFold batch prediction run). It uses Python to generate a command script for UCSF ChimeraX, aligning the structures and encoding the morph transition into a video file.

### Project Features

Batch Structure Processing: Automatically searches a parent directory for subfolders, treating each subfolder as an evolutionary snapshot or model.

Randomized Input: If a subfolder contains multiple Predicted Structures (.cif or .cif.gz), the script randomly selects one to represent that state.

ChimeraX Script Generation: Creates a temporary ChimeraX command script (.cxc) to handle the entire visualization workflow.

Automated Alignment: Uses the highly robust match command in ChimeraX to align all protein structures based on sequence and C-alpha coordinates.

MP4 Video Export: Encodes the final aligned morph animation as an MP4 video file.

### Prerequisites

To run this workflow, you need the following installed on your system:

Python 3.6+

UCSF ChimeraX: Must be installed and configured.

Required Python Libraries: pathlib (standard), shutil (standard), and subprocess (standard)



### Command 
After excecuting the previous command, use the output fasta file to generate Alphafold Structures using google Deepmind Alphafold 3 Server

After obtaining the structures, extract it into the same folder and use this

```
python main2.py "<save location for the morph.mp4>"  
```

