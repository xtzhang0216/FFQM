#!/usr/bin/env python
import os
import re
import argparse
import glob
from pathlib import Path

# Dictionary of atomic numbers to element symbols
ATOMIC_NUMBERS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca'
}

def check_optimization_convergence(log_content):
    """
    Check if the optimization converged successfully

    Parameters:
    -----------
    log_content : str
        Content of the Gaussian log file

    Returns:
    --------
    tuple
        (bool, str) - (converged, message)
    """
    # Check for normal termination
    if "Normal termination" not in log_content:
        return False, "Abnormal termination detected"

    # Method 1: Check for optimization completion phrases
    if "Stationary point found" in log_content or "Optimization completed" in log_content:
        return True, "Optimization completed successfully"

    # Method 2: Check for optimization convergence criteria
    convergence_pattern = r'Maximum Force\s+(\S+)\s+(\S+)\s+(\S+)\s*\nRMS     Force\s+(\S+)\s+(\S+)\s+(\S+)\s*\nMaximum Displacement\s+(\S+)\s+(\S+)\s+(\S+)\s*\nRMS     Displacement\s+(\S+)\s+(\S+)\s+(\S+)\s*'
    convergence_match = re.search(convergence_pattern, log_content)

    if convergence_match:
        # Extract convergence values and check if all are YES
        try:
            converged = all(convergence_match.group(i) == "YES" for i in [3, 6, 9, 12])

            if converged:
                return True, "Optimization converged successfully"
            else:
                # Create detailed message about which criteria didn't converge
                criteria = ["Maximum Force", "RMS Force", "Maximum Displacement", "RMS Displacement"]
                values = [convergence_match.group(i) for i in [1, 4, 7, 10]]
                thresholds = [convergence_match.group(i) for i in [2, 5, 8, 11]]
                converged_criteria = [convergence_match.group(i) for i in [3, 6, 9, 12]]

                message = "Optimization did not fully converge:\n"
                for i in range(4):
                    message += f"  {criteria[i]}: {values[i]} (Threshold: {thresholds[i]}) - {converged_criteria[i]}\n"

                return False, message
        except (IndexError, AttributeError):
            pass  # If there's an error parsing the convergence info, try other methods

    # Method 3: Check for SCF convergence
    if "SCF Done" in log_content:
        return True, "SCF calculation completed"

    # If we can't determine convergence but the file has normal termination, assume it's OK
    if "Normal termination" in log_content:
        return True, "Normal termination detected"

    return False, "Convergence information not found"

def extract_optimized_geometry(log_file):
    """
    Extract the final optimized geometry from a Gaussian log file

    Parameters:
    -----------
    log_file : str
        Path to the Gaussian log file

    Returns:
    --------
    tuple
        (list of tuples, bool, str) - (atoms, converged, message)
        atoms: List of (atomic_number, x, y, z) tuples representing the optimized geometry
        converged: Whether the optimization converged
        message: Convergence message
    """
    with open(log_file, 'r') as f:
        content = f.read()

    # Check convergence
    converged, message = check_optimization_convergence(content)

    # Try different patterns to find the optimized geometry

    # Pattern 1: Standard orientation
    orientation_sections = re.findall(
        r'Standard orientation:\s*----+\s*\n(?:.*)\s*----+\s*\n(.*?)\n\s*----+',
        content,
        re.DOTALL
    )

    # # Pattern 2: Input orientation (if standard orientation not found)
    # if not orientation_sections:
    #     orientation_sections = re.findall(r'Input orientation:.*?---------------------------------------------------------------------\n(.*?)----', content, re.DOTALL)

    # # Pattern 3: Z-Matrix orientation (if neither standard nor input orientation found)
    # if not orientation_sections:
    #     orientation_sections = re.findall(r'Z-Matrix orientation:.*?---------------------------------------------------------------------\n(.*?)----', content, re.DOTALL)

    # # Pattern 4: Look for final structure in optimization
    # if not orientation_sections:
    #     orientation_sections = re.findall(r'Optimization completed.*?Standard orientation:.*?---------------------------------------------------------------------\n(.*?)----', content, re.DOTALL)

    if not orientation_sections:
        print(f"Warning: No molecular orientation found in {log_file}")
        return None, converged, message

    # Get the last one (final optimized geometry)
    last_orientation = orientation_sections[-1]

    # Parse the coordinates
    atoms = []
    for line in last_orientation.strip().split('\n'):
        parts = line.split()
        if len(parts) >= 6:  # Ensure the line has enough columns
            try:
                center_num = int(parts[0])
                atomic_num = int(parts[1])
                atom_type = int(parts[2])
                x = float(parts[3])
                y = float(parts[4])
                z = float(parts[5])
                atoms.append((atomic_num, x, y, z))
            except (ValueError, IndexError):
                continue

    return atoms, converged, message

def write_pdb(atoms, output_file, molecule_name):
    """
    Write atoms to a PDB file

    Parameters:
    -----------
    atoms : list of tuples
        List of (atomic_number, x, y, z) tuples
    output_file : str
        Path to the output PDB file
    molecule_name : str
        Name of the molecule (used for residue name)
    """
    with open(output_file, 'w') as f:
        f.write(f"TITLE     {molecule_name} optimized with Gaussian\n")
        f.write("MODEL     1\n")

        for i, (atomic_num, x, y, z) in enumerate(atoms, 1):
            element = ATOMIC_NUMBERS.get(atomic_num, 'X')
            atom_name = element
            if len(element) == 1:
                atom_name = f" {element}"

            # Format according to PDB standard
            f.write(f"ATOM  {i:5d} {atom_name:<4s} {molecule_name:3s}    1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n")

        f.write("TER\n")
        f.write("ENDMDL\n")
        f.write("END\n")

def main():
    parser = argparse.ArgumentParser(description="Extract optimized geometry from Gaussian log files and convert to PDB")
    parser.add_argument("--input_dir", default="/pubhome/xtzhang/interaction/FF/geometry_opt/gaussian_results",
                        help="Directory containing Gaussian log files")
    parser.add_argument("--output_dir", default="/pubhome/xtzhang/interaction/FF/geometry_opt/optimized_pdbs",
                        help="Directory to save output PDB files")
    parser.add_argument("--molecule", default=None,
                        help="Specific molecule to process (without extension). If not provided, processes all molecules")

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Get input files
    if args.molecule:
        input_files = [os.path.join(args.input_dir, f"{args.molecule}.log")]
        if not os.path.exists(input_files[0]):
            print(f"Error: Log file for {args.molecule} not found")
            return
    else:
        input_files = glob.glob(os.path.join(args.input_dir, "*.log"))

    # Sort input files to process them in a consistent order
    input_files.sort()

    # Process each log file
    for log_file in input_files:
        base_name = os.path.basename(log_file).split('.')[0]
        output_file = os.path.join(args.output_dir, f"{base_name}.pdb")

        print(f"Processing {base_name}...")

        # Extract optimized geometry
        result = extract_optimized_geometry(log_file)

        # Check if result is a tuple (new format) or just atoms (old format)
        if isinstance(result, tuple) and len(result) == 3:
            atoms, converged, message = result

            # Print convergence information
            if converged:
                print(f"  Optimization converged: {message}")
            else:
                print(f"  Warning: {message}")
        else:
            # Handle old format (just atoms)
            atoms = result

        if atoms:
            # Write to PDB file
            write_pdb(atoms, output_file, base_name)
            print(f"  Optimized geometry saved to {output_file}")
        else:
            print(f"  Failed to extract optimized geometry from {log_file}")

    print("\nAll files processed.")

if __name__ == "__main__":
    main()
