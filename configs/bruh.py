import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED, Crippen, Descriptors, Lipinski
import sascorer  # Make sure sascorer.py and fpscores.pkl.gz are in the same folder

# Define input SMILES
smiles_list = [
    "CC1=NC2=NC(=NN2C(NC2=CC=C(C=C2)S(F)(F)(F)(F)F)=C1)C(C)(F)F",
    "N(C=1N2C(=NC(C(C)(F)F)=N2)N=C(C)C1)C=3C=CC(C(F)(F)F)=NC3"
]

# Prepare output list
output = []

for smi in smiles_list:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        qed = QED.qed(mol)
        logp = Crippen.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        sa = sascorer.calculateScore(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        h_donors = Lipinski.NumHDonors(mol)

        output.append({
            "SMILES": smi,
            "QED": round(qed, 4),
            "SlogP": round(logp, 4),
            "Molecular weight": round(mw, 2),
            "SA score": round(sa, 4),
            "H-bond acceptors": h_acceptors,
            "H-bond donors": h_donors
        })

# Save to CSV
df = pd.DataFrame(output)
df.to_csv("molecule_scores.csv", index=False)
print("Saved to molecule_scores.csv")
