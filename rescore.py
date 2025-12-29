import pandas as pd 
import os 
import numpy as np




def list_docked_pdbqt_files(directory):
    """
    Lists all files in the specified directory that start with 'docked_' and end with '.pdbqt'.

    :param directory: The directory to search for files.
    :return: A list of filenames that match the criteria.
    """
    # List to store matching filenames
    matching_files = []

    # Iterate over all files in the specified directory
    for filename in os.listdir(directory):
        # Check if the filename starts with 'docked_' and ends with '.pdbqt'
        if filename.startswith('docked_') and filename.endswith('.pdbqt') and "ligand" not in filename:
            # Add the matching filename to the list
            matching_files.append(os.path.join(directory , filename))
    
    return matching_files



def get_score_from_a_file( file_path):
        """
        Make a list of a ligands information including its docking score.

        Inputs:
        :param str file_path: the path to the file to be scored

        Returns:
        :returns: list lig_info: a list containing all info from
            self.smiles_dict for a given ligand and the ligands short_id_name and
            the docking score from the best pose.
        """
        # grab the index of the ligand for the score
        basefile = os.path.basename(file_path)
        basefile_strip = basefile.replace(".pdbqt", "").replace("docked_", "")
        basefile_split = basefile.split("__")
        ligand_short_name = basefile_split[0].replace("docked_", "")

        affinity = None
        affinity_complex = None
        with open(file_path, "r") as f:
            for line in f.readlines():
                 if "REMARK VINA" in line:
                    line_stripped = line.replace("REMARK VINA RESULT:", "").replace(
                        "\n", ""
                    )
                
                    line_split = line_stripped.split()

                    if affinity is None:
                        affinity = float(line_split[0])
                    else:
                        if affinity > float(line_split[0]):
                            affinity = float(line_split[0])
        if affinity is None:
            # This file lacks a pose to use
            return None
        with open(os.path.split(file_path)[0]+"/complex_"+os.path.split(file_path)[1], "r") as f:
            for line in f.readlines():
               if "Affinity" in line:
                    line_stripped = line.replace("Affinity:", "").replace(
                        "\n", ""
                    )
                    line_split = line_stripped.split()

                    if affinity_complex is None:
                        affinity_complex = float(line_split[0])
                    else:
                        if affinity_complex > float(line_split[0]):
                            affinity_complex = float(line_split[0])

        # Obtain additional file
        r=0.001987
        t = 298.0


        kd1 = np.exp(affinity/(r*t))
        kd2 = np.exp(affinity_complex/(r*t))
        alpha = kd2/kd1
        file = open(os.path.split(file_path)[0]+"/complex_"+os.path.split(file_path)[1], "a")
        sc = -((-np.log(kd2)) * (1 - 1/(alpha+1)))
        file.write(str([affinity,affinity_complex,kd1 , kd2 , alpha , sc]))
        file.close()
        

        lig_info = [ligand_short_name, basefile_strip, kd2 , kd1 , alpha, sc]
        return lig_info
        

def process_docking_files(directory):
    """
    Processes all 'docked_*.pdbqt' files in the specified directory and saves the results to a CSV file.

    :param directory: The directory containing the docking files.
    :param smiles_dict: A dictionary containing ligand information.
    :return: None
    """
    # List all matching files
    files = list_docked_pdbqt_files(directory)

    # List to hold all lig_info
    all_lig_info = []

    # Iterate over each file
    for file_path in files:
        lig_info = get_score_from_a_file(file_path)
        if lig_info is not None:
            all_lig_info.append(lig_info)

    # Create DataFrame
    df = pd.DataFrame(all_lig_info, columns=['Ligand Short Name', 'Base File Name', 'KD1', 'KD2', 'Alpha', 'SC'])

    # Save to CSV
    df.to_csv(os.path.join(".", 'docking_scores.csv'), index=False)

direct = "/home/abdou/autgorw/complex/out/Run_28/generation_8/PDBs"

process_docking_files(direct)