from rdkit import Chem
from rdkit.Chem import AllChem

def smi_to_sdf(input_file, output_file):
    try:
        with open(input_file, 'r', encoding='utf-8') as infile:
            unique_smiles = []
            seen_smiles = set()
            total_lines = 0
            for line in infile:
                total_lines += 1
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                smiles = line.split()[0]
                if smiles not in seen_smiles:
                    unique_smiles.append(smiles)
                    seen_smiles.add(smiles)
            # Now, convert unique SMILES to molecules and write to .sdf
            writer = Chem.SDWriter(output_file)
            count_written = 0
            invalid_smiles = 0
            for smiles in unique_smiles:
                mol = Chem.MolFromSmiles(smiles)
                hmol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(hmol)
                hmol.GetConformer(0)
                mp = AllChem.MMFFGetMoleculeProperties(hmol)
                ff = AllChem.MMFFGetMoleculeForceField(hmol, mp)
                Chem.rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(hmol)
                if hmol is not None:
                    hmol.SetProp('_Name', "molecule_"+str(count_written))
                    writer.write(hmol)
                    count_written += 1
                else:
                    invalid_smiles += 1
            # Logging
            print(f"Total lines read: {total_lines}")
            print(f"Unique SMILES: {len(unique_smiles)}")
            print(f"Invalid SMILES: {invalid_smiles}")
            print(f"Molecules written to SDF: {count_written}")
            writer.close()
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Replace 'input.smi' and 'output.sdf' with your file paths
    file = input("where is the file ")
    smi_to_sdf(file, 'library.sdf')