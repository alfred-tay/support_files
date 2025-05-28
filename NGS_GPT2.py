import os
import subprocess
import shutil
import glob
import json
import re
import csv
from collections import defaultdict
import textwrap # Import textwrap for dedenting awk script

# --- Configuration ---
# Define the base directory where all processing will occur.
# It's recommended to run this script from the directory containing your raw FastQ files.
BASE_DIR = os.getcwd()

# Define subdirectories for organizing outputs
ORIGINAL_DIR = os.path.join(BASE_DIR, "Original_dir")
WORKING_DIR = os.path.join(BASE_DIR, "working_dir") # New working directory for active processing
TRIMMED_DIR = os.path.join(BASE_DIR, "Trimmed_dir")
FASTQC_REPORT_DIR = os.path.join(BASE_DIR, "FastQC_reports")
CONTIGS_SPADES_DIR = os.path.join(BASE_DIR, "Contigs_spades")
CONTIGS_SHOVILL_DIR = os.path.join(BASE_DIR, "Contigs_shovill")
# New Contigs_dir for all assembled genomes at the same level as Trimmed_dir and working_dir
CONTIGS_DIR = os.path.join(BASE_DIR, "Contigs_dir") 

# Minimum contig length for assembly filtering (used after assembly for SPAdes, and by Shovill directly)
MIN_CONTIG_LENGTH = 1000

# --- Helper Functions ---

def run_command(command, description="", stdout_redirect_path=None, stderr_redirect_path=None):
    """
    Executes a shell command and checks for errors.
    Args:
        command (list): A list of strings representing the command and its arguments.
        description (str): A description of the command being run for logging.
        stdout_redirect_path (str, optional): Path to redirect stdout to.
        stderr_redirect_path (str, optional): Path to redirect stderr to.
    Returns:
        bool: True if command succeeded, False otherwise.
    """
    print(f"\n--- Running: {description} ---")
    print(f"Command: {' '.join(command)}")
    
    process_successful = False
    stdout_file = None
    stderr_file = None

    try:
        # Check for redirection operator in the command list
        # If found, remove it and the target file from the command list
        # and open the file for stdout redirection.
        actual_command_parts = list(command) # Create a mutable copy
        if ">" in actual_command_parts:
            output_file_index = actual_command_parts.index(">")
            output_filepath = actual_command_parts[output_file_index + 1]
            actual_command_parts = actual_command_parts[:output_file_index] # Command without redirection parts
            
            # Determine mode based on file extension for proper handling of gzipped files
            # For cat, 'wb' mode is safer as it's binary data
            mode = 'wb'
            
            stdout_file = open(output_filepath, mode)

        # Handle explicit stdout/stderr redirection via arguments
        if stdout_redirect_path and not stdout_file: # Only open if not already opened by '>'
            stdout_file = open(stdout_redirect_path, 'wb') # Use 'wb' for binary write (good for raw output or gzipped)
        
        if stderr_redirect_path:
            # If stderr is redirected to the same file as stdout, use the same file object
            if stdout_redirect_path and stdout_redirect_path == stderr_redirect_path and stdout_file:
                stderr_file = stdout_file
            else:
                stderr_file = open(stderr_redirect_path, 'wb') # Use 'wb' for binary write

        result = subprocess.run(
            actual_command_parts, # Use the modified command list
            check=True,
            stdout=stdout_file if stdout_file else subprocess.PIPE,
            stderr=stderr_file if stderr_file else subprocess.PIPE,
            # Use text=False if any redirection is happening (as we opened in 'wb' mode),
            # otherwise text=True for capturing to PIPE for console output.
            text=False if (stdout_file or stderr_file) else True
        )
        
        # Print captured stdout/stderr to console if not redirected to a file
        if not stdout_file and result.stdout:
            print("STDOUT:\n", result.stdout.decode(errors='ignore') if isinstance(result.stdout, bytes) else result.stdout)
        if not stderr_file and result.stderr:
            print("STDERR:\n", result.stderr.decode(errors='ignore') if isinstance(result.stderr, bytes) else result.stderr)

        print(f"--- {description} completed successfully. ---")
        process_successful = True
    except FileNotFoundError:
        print(f"Error: Command '{command[0]}' not found. Please ensure it's installed and in your PATH.")
    except subprocess.CalledProcessError as e:
        print(f"Error: {description} failed with exit code {e.returncode}")
        # Print captured stdout/stderr from the error object if not redirected to file
        if not stdout_file and e.stdout:
            print("STDOUT:\n", e.stdout.decode(errors='ignore') if isinstance(e.stdout, bytes) else e.stdout)
        if not stderr_file and e.stderr:
            print("STDERR:\n", e.stderr.decode(errors='ignore') if isinstance(e.stderr, bytes) else e.stderr)
    except Exception as e:
        print(f"An unexpected error occurred while running {description}: {e}")
    finally:
        # Close file handles if they were opened and are not the same object
        if stdout_file and stdout_file is not stderr_file:
            stdout_file.close()
        if stderr_file and stderr_file is not stdout_file: # Ensure not to double close if they are the same
            stderr_file.close()
    
    return process_successful

def ensure_dir(path):
    """Ensures that a directory exists, creating it if necessary."""
    os.makedirs(path, exist_ok=True)
    print(f"Ensured directory exists: {path}")

# --- Core Functions ---

def check_and_install_software(software_name, install_command):
    """
    Checks if a software is installed and attempts to install it if not.
    Assumes apt-get for installation. User might need sudo privileges.
    Args:
        software_name (str): The name of the software (e.g., "fastp").
        install_command (list): The apt-get command to install the software.
    Returns:
        bool: True if software is available/installed, False otherwise.
    """
    print(f"\nChecking for {software_name}...")
    try:
        # Check if the command exists in PATH
        subprocess.run(["which", software_name], check=True, capture_output=True)
        print(f"{software_name} is already installed.")
        return True
    except subprocess.CalledProcessError:
        print(f"{software_name} not found. Attempting to install...")
        print(f"Please enter your sudo password if prompted for '{' '.join(install_command)}'.")
        if run_command(install_command, f"Installing {software_name}"):
            print(f"{software_name} installed successfully.")
            return True
        else:
            print(f"Failed to install {software_name}. Please install it manually.")
            return False
    except Exception as e:
        print(f"An error occurred while checking/installing {software_name}: {e}")
        return False

def find_and_combine_fastq_files():
    """
    Finds FastQ files, identifies lanes, combines them per strain,
    and places them in Original_dir (archive) and WORKING_DIR (for processing).
    Returns:
        dict: A dictionary mapping strain_ID to a tuple (R1_path_in_WORKING_DIR, R2_path_in_WORKING_DIR).
    """
    print("\n--- Checking for FastQ files and combining lanes ---")
    # Only glob for files directly in BASE_DIR, not in subdirectories
    fastq_files = glob.glob(os.path.join(BASE_DIR, "*.fastq*"))
    
    if not fastq_files:
        print("No FastQ files found in the current directory. Please ensure your files are here.")
        return {}

    file_groups = defaultdict(lambda: defaultdict(list)) # {strain_id: {R1: [file_paths], R2: [file_paths]}}

    for fq_file in fastq_files:
        filename = os.path.basename(fq_file)
        
        # Regex to extract strain_ID and read_type from 'STRAINID_LXXX_R[12].fastq.gz' or 'STRAINID_R[12].fastq.gz'
        # This regex will capture everything before '_LXXX_R[12]' or '_R[12]' as the strain_id
        match = re.match(r"^(.*?)_L\d{3}_(R[12])(\.fastq\.gz|\.fastq)$", filename)
        if match:
            strain_id = match.group(1)
            read_type = match.group(2)
        else:
            # Fallback for files that might not have lane info but still have R1/R2
            match = re.match(r"^(.*?)_(R[12])(\.fastq\.gz|\.fastq)$", filename)
            if match:
                strain_id = match.group(1)
                read_type = match.group(2)
            else:
                print(f"Warning: Could not parse filename for strain ID and read type: {filename}. Skipping.")
                continue
        
        print(f"  Parsed '{filename}' -> Strain ID: '{strain_id}', Read Type: '{read_type}'")
        file_groups[strain_id][read_type].append(fq_file)

    combined_files_in_working_dir = {} # This will store paths to files in WORKING_DIR for downstream use
    ensure_dir(ORIGINAL_DIR)
    ensure_dir(WORKING_DIR)

    for strain_id, reads in file_groups.items():
        r1_files = sorted(reads.get("R1", []))
        r2_files = sorted(reads.get("R2", []))

        if not r1_files or not r2_files:
            print(f"Warning: Missing R1 or R2 files for strain {strain_id}. Skipping combination.")
            continue

        print(f"\nProcessing strain: {strain_id}")
        print(f"  R1 files ({len(r1_files)} lanes): {', '.join([os.path.basename(f) for f in r1_files])}")
        print(f"  R2 files ({len(r2_files)} lanes): {', '.join([os.path.basename(f) for f in r2_files])}")

        # Define paths for files in WORKING_DIR (for downstream processing)
        # These names will NOT contain L00X
        output_r1_working_path = os.path.join(WORKING_DIR, f"{strain_id}_R1.fastq.gz")
        output_r2_working_path = os.path.join(WORKING_DIR, f"{strain_id}_R2.fastq.gz")

        # Define paths for files in Original_dir (archive)
        # These names will also NOT contain L00X, as they are copies of the combined files
        output_r1_original_path = os.path.join(ORIGINAL_DIR, f"{strain_id}_R1.fastq.gz")
        output_r2_original_path = os.path.join(ORIGINAL_DIR, f"{strain_id}_R2.fastq.gz")

        # Combine R1 files into WORKING_DIR using cat
        # Ensure the output file is gzipped if input files are gzipped
        command_r1 = ["cat"] + r1_files
        if not run_command(command_r1, f"Combine R1 for {strain_id} into {WORKING_DIR}", stdout_redirect_path=output_r1_working_path):
            print(f"Error: Failed to combine R1 files for {strain_id}. Skipping this strain.")
            continue # Skip to next strain if combination fails

        # Combine R2 files into WORKING_DIR using cat
        command_r2 = ["cat"] + r2_files
        if not run_command(command_r2, f"Combine R2 for {strain_id} into {WORKING_DIR}", stdout_redirect_path=output_r2_working_path):
            print(f"Error: Failed to combine R2 files for {strain_id}. Skipping this strain.")
            continue # Skip to next strain if combination fails

        print(f"  Combined files for {strain_id} are in working directory: {os.path.basename(output_r1_working_path)}, {os.path.basename(output_r2_working_path)}")

        # Copy combined files from WORKING_DIR to ORIGINAL_DIR for archive
        try:
            shutil.copy(output_r1_working_path, output_r1_original_path)
            shutil.copy(output_r2_working_path, output_r2_original_path)
            print(f"  Archived combined files to Original_dir: {os.path.basename(output_r1_original_path)}, {os.path.basename(output_r2_original_path)}")
        except Exception as e:
            print(f"Error archiving combined files to Original_dir for {strain_id}: {e}")
            # Do not 'continue' here, as the files are already in working_dir and can be used downstream.
            # Just log the error.

        combined_files_in_working_dir[strain_id] = (output_r1_working_path, output_r2_working_path)
    
    print("\n--- Finished combining FastQ files ---")
    return combined_files_in_working_dir

def run_fastp(strain_id, r1_input, r2_input):
    """
    Runs fastp for trimming and quality filtering.
    Args:
        strain_id (str): The ID of the strain.
        r1_input (str): Path to the input R1 FastQ file.
        r2_input (str): Path to the input R2 FastQ file.
    Returns:
        tuple: (path_to_trimmed_r1, path_to_trimmed_r2, path_to_fastp_json_report, path_to_trim_stat_report) or (None, None, None, None) on failure.
    """
    print(f"\n--- Running fastp for {strain_id} ---")
    ensure_dir(TRIMMED_DIR)

    # Pre-check: Ensure input files are not empty
    if not os.path.exists(r1_input) or os.path.getsize(r1_input) == 0:
        print(f"Error: Input R1 file for fastp is missing or empty: {r1_input}. Skipping fastp for {strain_id}.")
        return None, None, None, None
    if not os.path.exists(r2_input) or os.path.getsize(r2_input) == 0:
        print(f"Error: Input R2 file for fastp is missing or empty: {r2_input}. Skipping fastp for {strain_id}.")
        return None, None, None, None

    output_r1 = os.path.join(TRIMMED_DIR, f"{strain_id}_trimmed_R1.fastq.gz")
    output_r2 = os.path.join(TRIMMED_DIR, f"{strain_id}_trimmed_R2.fastq.gz")
    json_report = os.path.join(TRIMMED_DIR, f"{strain_id}_fastp_report.json")
    html_report = os.path.join(TRIMMED_DIR, f"{strain_id}_fastp_report.html")
    trim_stat_report = os.path.join(TRIMMED_DIR, f"{strain_id}.trim_stat") # New trim stat file

    command = [
        "fastp",
        "-i", r1_input,
        "-o", output_r1,
        "-I", r2_input,
        "-O", output_r2,
        "-j", json_report,
        "-h", html_report,
        "-g",                      # Use -g for polyG trimming (equivalent to --trim_poly_g)
        "--poly_g_min_len", "10",  # Minimum length for polyG trimming
        "--detect_adapter_for_pe"  # Detect and trim adapters for paired-end reads
    ]

    # Run fastp, redirecting both stdout and stderr to the trim_stat_report file
    if run_command(command, f"fastp for {strain_id}", stdout_redirect_path=trim_stat_report, stderr_redirect_path=trim_stat_report):
        print(f"fastp completed for {strain_id}. Trimmed files: {output_r1}, {output_r2}")
        print(f"fastp raw output (stdout/stderr) saved to: {trim_stat_report}")
        return output_r1, output_r2, json_report, trim_stat_report
    else:
        print(f"fastp failed for {strain_id}.")
        print("Debugging Tip: If fastp consistently fails, check the content of the .trim_stat file for detailed error messages directly from fastp.")
        return None, None, None, None

def run_fastqc(strain_id, trimmed_r1, trimmed_r2):
    """
    Runs FastQC on trimmed FastQ files.
    Args:
        strain_id (str): The ID of the strain.
        trimmed_r1 (str): Path to the trimmed R1 FastQ file.
        trimmed_r2 (str): Path to the trimmed R2 FastQ file.
    Returns:
        bool: True if FastQC succeeded, False otherwise.
    """
    print(f"\n--- Running FastQC for {strain_id} ---")
    ensure_dir(FASTQC_REPORT_DIR)

    command = [
        "fastqc",
        trimmed_r1,
        trimmed_r2,
        "-o", FASTQC_REPORT_DIR
    ]

    if run_command(command, f"FastQC for {strain_id}"):
        print(f"FastQC completed for {strain_id}. Reports in {FASTQC_REPORT_DIR}")
        return True
    else:
        print(f"FastQC failed for {strain_id}.")
        return False

def rename_contigs(input_fasta_path, strain_id, output_fasta_path):
    """
    Renames contigs in a FASTA file from >NODES... to >Strain_ID_Contig0001.
    Args:
        input_fasta_path (str): Path to the input FASTA file.
        strain_id (str): The ID of the strain.
        output_fasta_path (str): Path where the renamed FASTA will be saved.
    Returns:
        bool: True if renaming succeeded, False otherwise.
    """
    print(f"\n--- Renaming contigs for {strain_id} in {os.path.basename(input_fasta_path)} ---")
    try:
        with open(input_fasta_path, 'r') as infile, open(output_fasta_path, 'w') as outfile:
            contig_count = 0
            for line in infile:
                if line.startswith(">"):
                    contig_count += 1
                    # Format contig number with leading zeros (e.g., 0001, 0010, 0100)
                    new_header = f">{strain_id}_Contig{contig_count:04d}\n"
                    outfile.write(new_header)
                else:
                    outfile.write(line)
        print(f"Contigs renamed and saved to {output_fasta_path}")
        return True
    except FileNotFoundError:
        print(f"Error: Input FASTA file not found: {input_fasta_path}")
        return False
    except Exception as e:
        print(f"An error occurred during contig renaming: {e}")
        return False

def filter_contigs_by_length(input_fasta, output_fasta, min_length):
    """
    Filters a FASTA file to keep only contigs greater than or equal to min_length.
    Uses awk for efficient filtering, robustly handling multi-line FASTA sequences.
    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path for the filtered output FASTA file.
        min_length (int): Minimum length of contigs to keep.
    Returns:
        bool: True if filtering succeeded, False otherwise.
    """
    print(f"\n--- Filtering contigs in {os.path.basename(input_fasta)} for length >= {min_length} bp ---")

    # This awk script correctly handles multi-line FASTA entries.
    # It collects header and sequence until a new header is found or end of file,
    # then checks the accumulated sequence length.
    awk_script = textwrap.dedent(f'''
        BEGIN {{
            RS=">"; # Record separator is '>'
            ORS=""; # Output record separator is empty, we will add newlines manually
        }}
        {{
            if (NR == 1) {{ # Skip the first empty record before the first '>'
                next;
            }}
            
            # Split the record by newline to separate header and sequence
            split($0, lines, "\\n");
            
            header = ">" lines[1]; # Reconstruct header with '>'
            
            # Reconstruct sequence from remaining lines
            sequence = "";
            for (i = 2; i <= length(lines); i++) {{
                sequence = sequence lines[i];
            }}
            
            # Remove any carriage returns from sequence data for accurate length
            gsub(/\\r/, "", sequence);
            
            if (length(sequence) >= {min_length}) {{
                print header "\\n" sequence "\\n";
            }}
        }}
    ''').strip() # .strip() removes any leading/trailing blank lines from textwrap

    awk_command = [
        "awk",
        awk_script,
    ]

    try:
        with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
            process = subprocess.run(
                awk_command,
                stdin=infile,
                stdout=outfile,
                stderr=subprocess.PIPE,
                check=True,
                text=True # Use text mode for awk input/output
            )
            if process.stderr:
                print("AWK STDERR:\n", process.stderr)
        print(f"Contigs filtered and saved to {output_fasta}")
        return True
    except FileNotFoundError:
        print(f"Error: awk command not found. Please ensure awk is installed and in your PATH.")
        return False
    except subprocess.CalledProcessError as e:
        print(f"Error: Contig filtering failed with exit code {e.returncode}")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        return False
    except Exception as e:
        print(f"An unexpected error occurred during contig filtering: {e}")
        return False


def run_spades_assembly(strain_id, trimmed_r1, trimmed_r2):
    """
    Runs genome assembly using SPAdes.
    Args:
        strain_id (str): The ID of the strain.
        trimmed_r1 (str): Path to the trimmed R1 FastQ file.
        trimmed_r2 (str): Path to the trimmed R2 FastQ file.
    Returns:
        str: Path to the assembled and filtered contigs FASTA file, or None on failure.
    """
    print(f"\n--- Running SPAdes assembly for {strain_id} ---")
    output_dir = os.path.join(CONTIGS_SPADES_DIR, strain_id)
    ensure_dir(output_dir)
    ensure_dir(CONTIGS_DIR) # Ensure the central Contigs_dir exists

    command = [
        "spades.py",
        "--cov-cutoff", "auto",
        "--pe1-1", trimmed_r1,
        "--pe1-2", trimmed_r2,
        "-o", output_dir,
        # Removed "--force-replace" as requested
    ]

    if run_command(command, f"SPAdes assembly for {strain_id}"):
        raw_contigs_fasta = os.path.join(output_dir, "contigs.fasta")
        if os.path.exists(raw_contigs_fasta):
            # 1. Rename contigs within the SPAdes output directory
            renamed_contigs_fasta = os.path.join(output_dir, f"{strain_id}_raw.fasta") # Temporarily rename to _raw
            if not rename_contigs(raw_contigs_fasta, strain_id, renamed_contigs_fasta):
                print(f"Failed to rename contigs for SPAdes assembly of {strain_id}.")
                return None
            
            # 2. Filter contigs by length
            filtered_contigs_fasta = os.path.join(output_dir, f"{strain_id}.fasta") # This will be the final filtered file
            if not filter_contigs_by_length(renamed_contigs_fasta, filtered_contigs_fasta, MIN_CONTIG_LENGTH):
                print(f"Failed to filter contigs by length for SPAdes assembly of {strain_id}.")
                return None

            # 3. Copy the *filtered* and renamed contigs to the central Contigs_dir
            try:
                shutil.copy(filtered_contigs_fasta, os.path.join(CONTIGS_DIR, f"{strain_id}.fasta"))
                print(f"Copied {strain_id}.fasta to {CONTIGS_DIR}")
            except Exception as e:
                print(f"Error copying {strain_id}.fasta to {CONTIGS_DIR} for SPAdes: {e}")
            
            return filtered_contigs_fasta # Return path to the filtered file for QUAST
        else:
            print(f"SPAdes assembly completed, but 'contigs.fasta' not found in {output_dir}.")
            return None
    else:
        print(f"SPAdes assembly failed for {strain_id}.")
        return None

def run_shovill_assembly(strain_id, trimmed_r1, trimmed_r2):
    """
    Runs genome assembly using Shovill.
    Args:
        strain_id (str): The ID of the strain.
        trimmed_r1 (str): Path to the trimmed R1 FastQ file.
        trimmed_r2 (str): Path to the trimmed R2 FastQ file.
    Returns:
        str: Path to the assembled contigs FASTA file, or None on failure.
    """
    print(f"\n--- Running Shovill assembly for {strain_id} ---")
    output_dir = os.path.join(CONTIGS_SHOVILL_DIR, strain_id)
    ensure_dir(output_dir)
    ensure_dir(CONTIGS_DIR) # Ensure the central Contigs_dir exists

    command = [
        "shovill",
        "--R1", trimmed_r1,
        "--R2", trimmed_r2,
        "--outdir", output_dir,
        "--minlen", str(MIN_CONTIG_LENGTH),
        "--mincov", "5.0",  # Added --mincov 5.0 as requested
        "--force"  # Force replace existing output directory
    ]

    if run_command(command, f"Shovill assembly for {strain_id}"):
        # Shovill outputs 'contigs.fa', need to check for this and rename
        shovill_raw_contigs_fa = os.path.join(output_dir, "contigs.fa") 
        contigs_fasta_standard_name = os.path.join(output_dir, "contigs.fasta") # Desired QUAST input name

        if os.path.exists(shovill_raw_contigs_fa):
            # Rename contigs.fa to contigs.fasta to standardize for downstream
            try:
                shutil.move(shovill_raw_contigs_fa, contigs_fasta_standard_name)
                print(f"Renamed {shovill_raw_contigs_fa} to {contigs_fasta_standard_name}")
            except Exception as e:
                print(f"Error renaming {shovill_raw_contigs_fa} to {contigs_fasta_standard_name}: {e}")
                return None
            
            # Now proceed with renaming contigs within the FASTA and copying
            renamed_contigs_fasta = os.path.join(output_dir, f"{strain_id}.fasta")
            if rename_contigs(contigs_fasta_standard_name, strain_id, renamed_contigs_fasta):
                # Copy the renamed contigs to the central Contigs_dir
                try:
                    shutil.copy(renamed_contigs_fasta, os.path.join(CONTIGS_DIR, f"{strain_id}.fasta"))
                    print(f"Copied {strain_id}.fasta to {CONTIGS_DIR}")
                except Exception as e:
                    print(f"Error copying {strain_id}.fasta to {CONTIGS_DIR} for Shovill: {e}")
                return renamed_contigs_fasta
            else:
                print(f"Failed to rename contigs for Shovill assembly of {strain_id}.")
                return None
        else:
            print(f"Shovill assembly completed, but 'contigs.fa' not found in {output_dir}.")
            return None
    else:
        print(f"Shovill assembly failed for {strain_id}.")
        return None

def run_quast(assembly_fasta_path, output_dir_base, strain_id):
    """
    Runs QUAST on an assembled FASTA file.
    Args:
        assembly_fasta_path (str): Path to the assembled FASTA file.
        output_dir_base (str): Base directory for QUAST output (e.g., Contigs_spades).
        strain_id (str): The ID of the strain.
    Returns:
        str: Path to the QUAST report.tsv, or None on failure.
    """
    print(f"\n--- Running QUAST for {strain_id} on {os.path.basename(assembly_fasta_path)} ---")
    quast_output_dir = os.path.join(output_dir_base, strain_id, "quast_report")
    ensure_dir(quast_output_dir)

    command = [
        "quast.py",
        assembly_fasta_path,
        "-o", quast_output_dir
    ]

    if run_command(command, f"QUAST for {strain_id}"):
        report_tsv_path = os.path.join(quast_output_dir, "report.tsv")
        if os.path.exists(report_tsv_path):
            print(f"QUAST report saved to {report_tsv_path}")
            return report_tsv_path
        else:
            print(f"QUAST completed, but 'report.tsv' not found in {quast_output_dir}.")
            return None
    else:
        print(f"QUAST failed for {strain_id}.")
        return None

def parse_fastp_report(json_report_path):
    """
    Parses a fastp JSON report to extract read counts.
    Args:
        json_report_path (str): Path to the fastp JSON report.
    Returns:
        tuple: (total_reads_before_trim, total_reads_after_trim, total_bases_after_trim) or (0, 0, 0) on failure.
    """
    try:
        with open(json_report_path, 'r') as f:
            report_data = json.load(f)
        
        before_trim_reads = report_data.get("summary", {}).get("before_filtering", {}).get("total_reads", 0)
        after_trim_reads = report_data.get("summary", {}).get("after_filtering", {}).get("total_reads", 0)
        total_bases_after_trim = report_data.get("summary", {}).get("after_filtering", {}).get("total_bases", 0)
        
        return before_trim_reads, after_trim_reads, total_bases_after_trim
    except FileNotFoundError:
        print(f"Error: fastp report not found: {json_report_path}")
        return 0, 0, 0
    except json.JSONDecodeError:
        print(f"Error: Could not decode fastp JSON report: {json_report_path}")
        return 0, 0, 0
    except Exception as e:
        print(f"An error occurred while parsing fastp report: {e}")
        return 0, 0, 0

def parse_trim_stat_report(trim_stat_report_path):
    """
    Parses the fastp .trim_stat file to extract total bases before and after filtering.
    Args:
        trim_stat_report_path (str): Path to the fastp .trim_stat file.
    Returns:
        tuple: (total_bases_before_filter, total_bases_after_filter) or (0, 0) on failure.
    """
    total_bases_before_filter = 0
    total_bases_after_filter = 0
    try:
        with open(trim_stat_report_path, 'r') as f:
            content = f.read()
            
            # Find all "total bases" under "Read1 before filtering" and "Read2 before filtering"
            before_filter_bases_matches = re.findall(r"Read[12] before filtering:\s*total reads:\s*\d+\s*total bases:\s*(\d+)", content)
            total_bases_before_filter = sum(int(b) for b in before_filter_bases_matches)

            # Find all "total bases" under "Read1 after filtering" and "Read2 after filtering"
            after_filter_bases_matches = re.findall(r"Read[12] after filtering:\s*total reads:\s*\d+\s*total bases:\s*(\d+)", content)
            total_bases_after_filter = sum(int(b) for b in after_filter_bases_matches)

    except FileNotFoundError:
        print(f"Error: fastp .trim_stat report not found: {trim_stat_report_path}")
    except Exception as e:
        print(f"An error occurred while parsing fastp .trim_stat report: {e}")
    return total_bases_before_filter, total_bases_after_filter

def parse_quast_report(quast_report_tsv_path):
    """
    Parses a QUAST report.tsv to extract contig number and genome size.
    Args:
        quast_report_tsv_path (str): Path to the QUAST report.tsv.
    Returns:
        tuple: (contigs_num, genome_size) or (0, 0) on failure.
    """
    contigs_num = 0
    genome_size = 0
    try:
        with open(quast_report_tsv_path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) >= 2:
                    # Use "# contigs (>= 0 bp)" for total contig count
                    if row[0].strip() == "# contigs (>= 0 bp)": 
                        contigs_num = int(row[1])
                    # Use "Total length (>= 0 bp)" for total genome size
                    elif row[0].strip() == "Total length (>= 0 bp)": 
                        genome_size = int(row[1])
        return contigs_num, genome_size
    except FileNotFoundError:
        print(f"Error: QUAST report not found: {quast_report_tsv_path}")
        return 0, 0
    except Exception as e:
        print(f"An error occurred while parsing QUAST report: {e}")
        return 0, 0

def generate_assembly_report(report_data, filename):
    """
    Generates a CSV assembly report.
    Args:
        report_data (dict): Dictionary containing all data for the report.
                            Keys are strain IDs, values are dicts of metrics.
        filename (str): The name of the CSV file to generate (e.g., "assembly_report_spades.csv").
    """
    report_path = os.path.join(BASE_DIR, filename)
    print(f"\n--- Generating Assembly Report: {filename} ---")
    headers = [
        "Strain ID",
        "Total Reads Before Trim",
        "Total Reads After Trim",
        "Total Bases Before Trim",
        "Total Bases After Trim",
        "Contigs Number",
        "Genome Size (bp)",
        "Coverage (X)"
    ]

    with open(report_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        for strain_id, data in report_data.items():
            total_reads_before_trim = data.get("total_reads_before_trim", 0)
            total_reads_after_trim = data.get("total_reads_after_trim", 0)
            total_bases_before_trim = data.get("total_bases_before_trim", 0)
            total_bases_after_trim = data.get("total_bases_after_trim", 0)
            contigs_num = data.get("contigs_num", 0)
            genome_size = data.get("genome_size", 0)
            

            coverage = 0
            if genome_size > 0 and total_bases_after_trim > 0:
                coverage = total_bases_after_trim / genome_size
            
            writer.writerow([
                strain_id,
                total_reads_before_trim,
                total_reads_after_trim,
                total_bases_before_trim,
                total_bases_after_trim,
                contigs_num,
                genome_size,
                f"{coverage:.2f}" # Format to two decimal places
            ])
    print(f"Assembly report generated: {report_path}")

# --- Main Workflow ---

def display_options():
    """Displays the available processing options to the user."""
    print("\n--- Available FastQ Processing Options ---")
    print("1. Check software availability")
    print("2. Check and combine reads of different lanes")
    print("3. Trimming - fastp")
    print("4. FastQC check on trimmed files")
    print("5. Assembly using SPAdes (includes QUAST)")
    print("6. Assembly using Shovill (includes QUAST)")
    print("7. Run QUAST (on existing SPAdes or Shovill assembly)")
    # Option 8 (Generate Assembly Report) is now automatic and removed from menu
    print("0. Exit")
    print("------------------------------------------")

def get_user_choices():
    """Gets user's choice of options to run."""
    while True:
        choices_str = input("Enter the numbers of the options you want to run (e.g., '1,2,3' or 'all' for all steps, '0' to exit): ").strip()
        if choices_str.lower() == '0':
            return []
        if choices_str.lower() == 'all':
            return list(range(1, 8)) # All options from 1 to 7 (excluding the now automatic option 8)
        
        try:
            choices = [int(x.strip()) for x in choices_str.split(',')]
            # Validate choices are within the valid range (1 to 7)
            if all(1 <= c <= 7 for c in choices):
                return sorted(list(set(choices))) # Return sorted unique choices
            else:
                print("Invalid option number detected. Please enter numbers between 1 and 7, or 'all', or '0'.")
        except ValueError:
            print("Invalid input format. Please enter comma-separated numbers, 'all', or '0'.")

def main_workflow():
    """Main function to orchestrate the FastQ processing workflow with user options."""
    
    print("Starting FastQ file processing workflow...")
    
    softwares_to_check = {
        "fastp": ["sudo", "apt-get", "install", "-y", "fastp"],
        "quast.py": ["sudo", "apt-get", "install", "-y", "quast"],
        "shovill": ["sudo", "apt-get", "install", "-y", "shovill"],
        "spades.py": ["sudo", "apt-get", "install", "-y", "spades"]
    }

    # Data structures to hold information for the final reports
    spades_assembly_report_data = {}
    shovill_assembly_report_data = {}
    
    # combined_files_for_processing will now store paths to files in WORKING_DIR
    combined_files_for_processing = {}
    
    # Ensure all main output directories exist at the start
    ensure_dir(ORIGINAL_DIR)
    ensure_dir(WORKING_DIR)
    ensure_dir(TRIMMED_DIR)
    ensure_dir(FASTQC_REPORT_DIR)
    ensure_dir(CONTIGS_SPADES_DIR)
    ensure_dir(CONTIGS_SHOVILL_DIR)
    ensure_dir(CONTIGS_DIR) # Ensure the new central Contigs_dir exists

    while True:
        display_options()
        selected_options = get_user_choices()

        if not selected_options:
            print("Exiting workflow.")
            break

        # combined_files_for_processing needs to be updated/re-checked if option 2 wasn't explicitly run.
        # This handles cases where user skips option 2 but expects to use existing files.
        if 2 not in selected_options:
            print("\n--- Checking for existing FastQ files in working_dir (Option 2 was not selected) ---")
            existing_working_files = glob.glob(os.path.join(WORKING_DIR, "*_R[12].fastq.gz"))
            if existing_working_files:
                current_combined_files = defaultdict(lambda: [None, None])
                for fq_file in existing_working_files:
                    filename = os.path.basename(fq_file)
                    match = re.match(r"^(.*)_(R[12])(\.fastq\.gz|\.fastq)$", filename) # Corrected regex
                    if match:
                        strain_id = match.group(1)
                        read_type = match.group(2)
                        if read_type == "R1":
                            current_combined_files[strain_id][0] = fq_file
                        else: # R2
                            current_combined_files[strain_id][1] = fq_file
                
                # Filter for complete pairs
                combined_files_for_processing = {
                    s_id: tuple(paths) for s_id, paths in current_combined_files.items()
                    if paths[0] and paths[1]
                }
                if not combined_files_for_processing:
                    print(f"No complete R1/R2 FastQ pairs found in {WORKING_DIR}. Cannot proceed with steps requiring them.")
                    # If files are missing and not running combine step, exit further processing options for this loop.
                    selected_options = [opt for opt in selected_options if opt <= 2] # Only allow software check and combine reads if it's the first time
                    continue
                else:
                    print(f"Found {len(combined_files_for_processing)} strains with FastQ files in {WORKING_DIR}.")
            else:
                print(f"No FastQ files found in {WORKING_DIR}. Please run option 2 first if you intend to process new files.")
                selected_options = [opt for opt in selected_options if opt <= 2]
                continue
        else: # Option 2 was selected, so run it
            print("\n--- Option 2: Checking and combining reads of different lanes ---")
            combined_files_for_processing = find_and_combine_fastq_files()
            if not combined_files_for_processing:
                print("No FastQ files found or processed in the current directory. Cannot proceed with further steps requiring these files.")
                selected_options = [opt for opt in selected_options if opt <= 2]
                continue


        if 1 in selected_options:
            # 1. Check and install necessary software
            print("\n--- Option 1: Checking software availability ---")
            all_software_installed = True
            for software, install_cmd in softwares_to_check.items():
                if not check_and_install_software(software, install_cmd):
                    all_software_installed = False
            if not all_software_installed:
                print("\nWarning: Not all required software could be installed. Subsequent steps might fail.")


        # Process each strain if combined files are available and relevant options are selected
        if combined_files_for_processing and any(opt in selected_options for opt in [3, 4, 5, 6, 7]):
            for strain_id, (r1_path, r2_path) in combined_files_for_processing.items():
                print(f"\nProcessing pipeline for strain: {strain_id}")
                
                # Initialize strain data for both reports if not already present
                if strain_id not in spades_assembly_report_data:
                    spades_assembly_report_data[strain_id] = {"strain_id": strain_id}
                if strain_id not in shovill_assembly_report_data:
                    shovill_assembly_report_data[strain_id] = {"strain_id": strain_id}

                trimmed_r1, trimmed_r2, fastp_json, trim_stat_report = None, None, None, None

                # Determine if trimming should be run or if existing trimmed files should be used
                if 3 in selected_options:
                    # 3. Trimming - fastp
                    print("\n--- Option 3: Trimming reads with fastp ---")
                    trimmed_r1, trimmed_r2, fastp_json, trim_stat_report = run_fastp(strain_id, r1_path, r2_path)
                    
                    if not (trimmed_r1 and trimmed_r2 and fastp_json and trim_stat_report):
                        print(f"Trimming failed for {strain_id}. Skipping subsequent steps for this strain.")
                        continue # Skip to next strain if trimming fails
                else: # Option 3 not selected, try to find existing trimmed files
                    trimmed_r1 = os.path.join(TRIMMED_DIR, f"{strain_id}_trimmed_R1.fastq.gz")
                    trimmed_r2 = os.path.join(TRIMMED_DIR, f"{strain_id}_trimmed_R2.fastq.gz")
                    fastp_json = os.path.join(TRIMMED_DIR, f"{strain_id}_fastp_report.json")
                    trim_stat_report = os.path.join(TRIMMED_DIR, f"{strain_id}.trim_stat")
                    
                    if not (os.path.exists(trimmed_r1) and os.path.exists(trimmed_r2) and os.path.exists(fastp_json) and os.path.exists(trim_stat_report)):
                        print(f"Warning: Trimmed files, fastp report, or trim_stat report not found for {strain_id} in {TRIMMED_DIR}. "
                              "Please run option 3 (Trimming) first, or ensure files exist if you intend to proceed.")
                        continue # Skip to next strain if trimmed files are missing
                    
                    print(f"Using existing trimmed files for {strain_id}.")

                # Always parse fastp report for assembly report data if trimmed files are available
                before_trim_reads, after_trim_reads, _ = parse_fastp_report(fastp_json)
                total_bases_before_trim_stat, total_bases_after_trim_stat = parse_trim_stat_report(trim_stat_report)

                # Store fastp data in both report dictionaries
                spades_assembly_report_data[strain_id]["total_reads_before_trim"] = before_trim_reads
                spades_assembly_report_data[strain_id]["total_reads_after_trim"] = after_trim_reads
                spades_assembly_report_data[strain_id]["total_bases_before_trim"] = total_bases_before_trim_stat
                spades_assembly_report_data[strain_id]["total_bases_after_trim"] = total_bases_after_trim_stat

                shovill_assembly_report_data[strain_id]["total_reads_before_trim"] = before_trim_reads
                shovill_assembly_report_data[strain_id]["total_reads_after_trim"] = after_trim_reads
                shovill_assembly_report_data[strain_id]["total_bases_before_trim"] = total_bases_before_trim_stat
                shovill_assembly_report_data[strain_id]["total_bases_after_trim"] = total_bases_after_trim_stat


                if 4 in selected_options:
                    # 4. FastQC check on trimmed files
                    print("\n--- Option 4: Running FastQC on trimmed files ---")
                    run_fastqc(strain_id, trimmed_r1, trimmed_r2)

                spades_contigs_fasta = None
                shovill_contigs_fasta = None
                
                # Check for existing contigs before attempting assembly
                existing_spades_fasta = os.path.join(CONTIGS_SPADES_DIR, strain_id, f"{strain_id}.fasta")
                existing_shovill_fasta = os.path.join(CONTIGS_SHOVILL_DIR, strain_id, f"{strain_id}.fasta")

                if 5 in selected_options:
                    # 5. Assembly using SPAdes (includes QUAST)
                    print("\n--- Option 5: Assembling with SPAdes ---")
                    spades_contigs_fasta = run_spades_assembly(strain_id, trimmed_r1, trimmed_r2)
                    
                    if spades_contigs_fasta:
                        # Run QUAST immediately after successful SPAdes assembly
                        quast_report_tsv = run_quast(spades_contigs_fasta, CONTIGS_SPADES_DIR, strain_id)
                        if quast_report_tsv:
                            contigs_num, genome_size = parse_quast_report(quast_report_tsv)
                            spades_assembly_report_data[strain_id]["contigs_num"] = contigs_num
                            spades_assembly_report_data[strain_id]["genome_size"] = genome_size
                        else:
                            print(f"QUAST report not available for SPAdes assembly of {strain_id}.")
                    else:
                        print(f"SPAdes assembly failed for {strain_id}. Skipping QUAST for this assembly.")
                elif os.path.exists(existing_spades_fasta):
                    spades_contigs_fasta = existing_spades_fasta
                    # If existing SPAdes assembly is used, try to parse its QUAST report for the data
                    quast_report_tsv = os.path.join(CONTIGS_SPADES_DIR, strain_id, "quast_report", "report.tsv")
                    if os.path.exists(quast_report_tsv):
                        contigs_num, genome_size = parse_quast_report(quast_report_tsv)
                        spades_assembly_report_data[strain_id]["contigs_num"] = contigs_num
                        spades_assembly_report_data[strain_id]["genome_size"] = genome_size


                if 6 in selected_options:
                    # 6. Assembly using Shovill (includes QUAST)
                    print("\n--- Option 6: Assembling with Shovill ---")
                    shovill_contigs_fasta = run_shovill_assembly(strain_id, trimmed_r1, trimmed_r2)
                    if shovill_contigs_fasta:
                        # Run QUAST immediately after successful Shovill assembly
                        quast_report_tsv = run_quast(shovill_contigs_fasta, CONTIGS_SHOVILL_DIR, strain_id)
                        if quast_report_tsv:
                            contigs_num, genome_size = parse_quast_report(quast_report_tsv)
                            # Store Shovill's QUAST results in the assembly report data
                            shovill_assembly_report_data[strain_id]["contigs_num"] = contigs_num
                            shovill_assembly_report_data[strain_id]["genome_size"] = genome_size
                            print(f"Shovill QUAST results for {strain_id}: Contigs={contigs_num}, Genome Size={genome_size}")
                        else:
                            print(f"QUAST report not available for Shovill assembly of {strain_id}.")
                    else:
                        print(f"Shovill assembly failed for {strain_id}.")
                elif os.path.exists(existing_shovill_fasta):
                    shovill_contigs_fasta = existing_shovill_fasta
                    # If existing Shovill assembly is used, try to parse its QUAST report for the data
                    quast_report_tsv = os.path.join(CONTIGS_SHOVILL_DIR, strain_id, "quast_report", "report.tsv")
                    if os.path.exists(quast_report_tsv):
                        contigs_num, genome_size = parse_quast_report(quast_report_tsv)
                        shovill_assembly_report_data[strain_id]["contigs_num"] = contigs_num
                        shovill_assembly_report_data[strain_id]["genome_size"] = genome_size


                if 7 in selected_options:
                    # 7. Run QUAST (on existing SPAdes or Shovill assembly)
                    print("\n--- Option 7: Running QUAST on existing assembly ---")
                    
                    quast_run_performed = False
                    if spades_contigs_fasta and os.path.exists(spades_contigs_fasta):
                        print(f"Running QUAST on SPAdes assembly for {strain_id}...")
                        quast_report_tsv = run_quast(spades_contigs_fasta, CONTIGS_SPADES_DIR, strain_id)
                        if quast_report_tsv:
                            contigs_num, genome_size = parse_quast_report(quast_report_tsv)
                            spades_assembly_report_data[strain_id]["contigs_num"] = contigs_num
                            spades_assembly_report_data[strain_id]["genome_size"] = genome_size
                            quast_run_performed = True
                        else:
                            print(f"QUAST report not available for SPAdes assembly of {strain_id}.")
                    
                    if shovill_contigs_fasta and os.path.exists(shovill_contigs_fasta) and not quast_run_performed:
                        # Only run QUAST on Shovill if SPAdes was not available/successful or if user prefers Shovill
                        print(f"Running QUAST on Shovill assembly for {strain_id}...")
                        quast_report_tsv = run_quast(shovill_contigs_fasta, CONTIGS_SHOVILL_DIR, strain_id)
                        if quast_report_tsv:
                            contigs_num, genome_size = parse_quast_report(quast_report_tsv)
                            shovill_assembly_report_data[strain_id]["contigs_num"] = contigs_num
                            shovill_assembly_report_data[strain_id]["genome_size"] = genome_size
                            quast_run_performed = True
                        else:
                            print(f"QUAST report not available for Shovill assembly of {strain_id}.")
                    
                    if not quast_run_performed:
                        print(f"No suitable existing SPAdes or Shovill assembly contigs found for {strain_id} to run QUAST on.")

    # Automatically generate the assembly reports at the end of the workflow
    if spades_assembly_report_data:
        generate_assembly_report(spades_assembly_report_data, "assembly_report_spades.csv")
    else:
        print("\nNo SPAdes assembly data collected. Skipping SPAdes assembly report generation.")
    
    if shovill_assembly_report_data:
        generate_assembly_report(shovill_assembly_report_data, "assembly_report_shovill.csv")
    else:
        print("\nNo Shovill assembly data collected. Skipping Shovill assembly report generation.")
    
    print("\nFastQ processing workflow completed.")

# Entry point for the script
if __name__ == "__main__":
    main_workflow()
