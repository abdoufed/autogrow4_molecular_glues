import os
import re
import math
import csv
import biotite.structure.io.pdbqt as pdbqt
import biotite

def extract_first_score(file_path):
    """
    Extract the first docking score from the Vina result file.
    """
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('REMARK VINA RESULT:'):
                # Extract the score which is the second element in the line
                score = float(line.split()[3])
                return score
    return None

def extract_smiles(file_path):
    """
    Extract the first docking score from the Vina result file.
    """
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('REMARK SMILES'):
                # Extract the score which is the second element in the line
                smile = str(line.split()[2])
                return smile
    return None

def calculate_kd(binding_energy, temperature=298):
    """
    Calculate the dissociation constant (Kd) from binding energy.
    """
    R = 0.001987  # kcal/(mol·K)
    kd = math.exp(binding_energy / (R * temperature))
    return kd
def calculate_rmsd(ref , lig ):
    pdbqt_file = pdbqt.PDBQTFile.read(ref)
    reference = pdbqt_file.get_structure()[0]
    pdbqt_file = pdbqt.PDBQTFile.read(lig)
    ligand = pdbqt_file.get_structure()[0]

    return biotite.structure.rmsd(reference,ligand)


    pass
def process_files(folder_complex, folder_crbn):
    """
    Process files in both folders and calculate required values.
    """
    results = []

    for filename in os.listdir(folder_complex):
        if filename.endswith('.pdbqt'):
            file_path_complex = os.path.join(folder_complex, filename)
            file_path_crbn = os.path.join(folder_crbn, filename)
            
            if os.path.exists(file_path_crbn):
                score_complex = extract_first_score(file_path_complex)
                score_crbn = extract_first_score(file_path_crbn)
                smiles = extract_smiles(file_path_crbn)
                
                if score_complex is not None and score_crbn is not None:
                    kd1 = calculate_kd(score_crbn)
                    kd2 = calculate_kd(score_complex)
                    alpha =   kd1 / kd2
                    rmsd = calculate_rmsd(file_path_crbn,file_path_complex)
                    
                    results.append({
                        'Ligand': filename,
                        'Smiles': smiles,
                        'Kd1 (CRBN)': kd1,
                        'Kd2 (Complex)': kd2,
                        'Binding Energy CRBN': score_crbn,
                        'Binding Energy Complex': score_complex,
                        'Alpha (Kd2/Kd1)': alpha,
                        'rmsd' : rmsd
                    })
    
    return results

def write_to_csv(results, output_file):
    """
    Write the results to a CSV file.
    """
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Ligand','Smiles', 'Kd1 (CRBN)', 'Kd2 (Complex)', 'Binding Energy CRBN', 'Binding Energy Complex', 'Alpha (Kd2/Kd1)','rmsd']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for result in results:
            writer.writerow(result)

def main():
    folder_complex = input("Enter the directory path for 'complex': ")
    folder_crbn = input("Enter the directory path for 'crbn': ")
    output_file = input("Enter the output CSV file path: ")
    
    if not os.path.isdir(folder_complex):
        print("Invalid directory path for 'complex'")
        return
    
    if not os.path.isdir(folder_crbn):
        print("Invalid directory path for 'crbn'")
        return
    
    results = process_files(folder_complex, folder_crbn)
    write_to_csv(results, output_file)
    
    print(f"Results written to {output_file}")

if __name__ == "__main__":
    main()
