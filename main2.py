import sys
import random
import subprocess
import argparse
import tempfile
from pathlib import Path
from shutil import which # Used to find executables in PATH
import os

# --- Configuration ---
# 1. Set the path to your ChimeraX executable here.
#    If ChimeraX is installed and added to your system's PATH, you can set this to 'chimerax'
#    or simply let the script try to find it automatically.
CHIMERAX_EXECUTABLE = which('chimerax') or r"C:\Program Files\ChimeraX 1.10.1\bin\ChimeraX.exe"
# 2. Desired output filename
OUTPUT_VIDEO_NAME = "morph_animation.mp4"
# ---------------------

def find_cif_files(folder: Path):
    """Return a list of all CIF files (including .cif.gz) inside a folder."""
    cifs = list(folder.glob("*.cif")) + list(folder.glob("*.cif.gz"))
    return cifs

def pick_random_cif(folder: Path):
    """Return a random CIF file (Path object) inside a folder."""
    cifs = find_cif_files(folder)
    return random.choice(cifs) if cifs else None


def generate_chimerax_script(structure_paths: list[Path]):
    """Return a string containing the ChimeraX commands to perform the morph."""
    lines = []
    lines.append("close all\n")
    lines.append("set bgColor white")
    lines.append("set silhouette true") # Enhance visual depth

    # 1. Open all structures
    for path in structure_paths:
        # Use str(path) for full OS-compatible path string
        lines.append(f"open '{path}'")

    # 2. Visualization and Alignment
    lines.append("cartoon") # Color by Secondary Structure (SSN)
    lines.append("view all") # Center and zoom to fit all structures
    
    # NEW ROBUST ALIGNMENT LOGIC:
    # Use the simpler 'match' command, which performs sequence alignment and superposition 
    # based on C-alpha atoms by default, resolving the previous syntax issues.
    
    # Loop through models 2 through N and align them to the reference model (#1)
    for i in range(2, len(structure_paths) + 1):
        # Match all atoms in model #i to all atoms in model #1.
        lines.append(f"match #{i} to #1")

    N = len(structure_paths)    
    

    # 3. Morph setup and recording
    all_models = ",".join([f"#{i}" if i == 1 else str(i) for i in range(1, N + 1)])
    
    # Start the morph, which automatically begins playing
    lines.append(f"morph {all_models} frames 60")
    lines.append(f"")
    
    # Record the 120 frames of the morph animation
    lines.append("movie record")
    lines.append("wait 541") # Wait for the 120 frames of the morph to complete
    lines.append("movie stop")

    # 4. Save output video
    OUTPUT_VIDEO_NAME = "morph.mp4"
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, OUTPUT_VIDEO_NAME)
    lines.append(f"movie encode output {output_path}")
    lines.append("exit") # Ensure ChimeraX closes cleanly

    return "\n".join(lines)


def main(base_folder_str):
    # Convert input string to a Path object
    base_folder = Path(base_folder_str).resolve()
    
    # Find ChimeraX executable path dynamically if possible
    global CHIMERAX_EXECUTABLE
    if not which('chimerax') and not Path(CHIMERAX_EXECUTABLE).exists():
        print(f"❌ ChimeraX executable not found at '{CHIMERAX_EXECUTABLE}'.")
        print("Please ensure ChimeraX is in your PATH or update the 'CHIMERAX_EXECUTABLE' variable.")
        return
    elif which('chimerax'):
        CHIMERAX_EXECUTABLE = which('chimerax')

    # 1. Collect CIFs
    structure_paths = []
    
    # Iterate over immediate subdirectories
    for sub in base_folder.iterdir():
        if sub.is_dir():
            cif = pick_random_cif(sub)
            if cif:
                structure_paths.append(cif)

    if len(structure_paths) < 2:
        print(f"❌ Need at least 2 folders in '{base_folder}' with CIF files to morph.")
        return

    print(f"✔ Found {len(structure_paths)} structures to morph:")
    for i, p in enumerate(structure_paths):
        print(f"   Model #{i+1}: {p.name}")

    # 2. Create temporary ChimeraX script
    cx_script_text = generate_chimerax_script(structure_paths)

    # Use Path to handle the temporary file creation
    temp_dir = Path(tempfile.gettempdir())
    script_path = temp_dir / "morph_script.cxc"
    
    try:
        script_path.write_text(cx_script_text, encoding='utf-8')
    except IOError as e:
        print(f"❌ Error writing script to temp file: {e}")
        return

    print(f"\n✔ Generated ChimeraX script at: {script_path}")

    # 3. Run ChimeraX
    command = [CHIMERAX_EXECUTABLE, "--nogui", str(script_path)]

    print(f"\n▶ Running ChimeraX (Command: {' '.join(command)})...")
    
    try:
        # Run the command and capture output/errors
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print("--- ChimeraX Output (last lines) ---")
        print('\n'.join(result.stdout.splitlines()[-5:]))
        print("-----------------------------------")
    
    except subprocess.CalledProcessError as e:
        print(f"❌ ChimeraX failed with return code {e.returncode}")
        print("--- Stderr ---")
        print(e.stderr)
        print("--------------")
        return
    except FileNotFoundError:
        # Should be caught by the initial check, but good for safety
        print(f"❌ Error: ChimeraX executable not found at '{CHIMERAX_EXECUTABLE}'")
        return

    # 4. Clean up temporary script
    try:
        script_path.unlink()
    except OSError:
        pass # Ignore if cleanup fails

    print(f"\n✨ Success! Output video saved as {OUTPUT_VIDEO_NAME} in the ChimeraX working directory.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a ChimeraX morph animation video from random CIF files found in subfolders."
    )
    parser.add_argument("folder", help="Parent folder containing subfolders, each with one or more CIF files.")
    args = parser.parse_args()

    main(args.folder)