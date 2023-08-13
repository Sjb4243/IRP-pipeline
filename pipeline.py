import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import subprocess

def clean_df(raw_df):
    not_na_df = raw_df[raw_df.standard_value.notna()]
    not_na_df = not_na_df[raw_df.canonical_smiles.notna()]
    not_na_df.drop_duplicates(["canonical_smiles"], inplace=True)
    selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
    not_na_df = not_na_df[selection]
    not_na_df.to_csv("processed_chembl.csv", index = False)
    return not_na_df

def remove_charges(cleaned_df):
    smiles_only = cleaned_df["canonical_smiles"]
    smiles_only.to_csv("babel.smi", index=False, header =False)
    cmd = "obabel babel.smi -O charges.smi -p 7"
    subprocess.call(cmd, shell=True)
    smiles_charges = pd.read_csv("charges.smi", names = ["canonical_smiles"])
    charges = []
    for smi in smiles_charges["canonical_smiles"]:
        try:
            mol = Chem.MolFromSmiles(smi)
            charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        except:
            print(f"Skipping {smi}")
            charges.append(-4)
            continue
        charges.append(charge)
    cleaned_df["charge"] = charges
    testing = pd.Series(charges)
    testing.to_csv("charges.csv", index = False)
    cleaned_df = cleaned_df[cleaned_df["charge"] > -3]
    return cleaned_df

def fingerprinting(pIC50_df, input_file_name):
    selection = ['canonical_smiles', 'molecule_chembl_id']
    smi_df = pIC50_df[selection]
    smi_dir = "/home/sjb176/IRP/ML/data/"
    smi_dir_full = smi_dir + "molecule.smi"
    smi_df.to_csv(smi_dir_full, sep='\t', index=False, header=False)
    cmd = "java -Xms1G -Xmx1G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar " \
          "-removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml " \
           "-dir " +  smi_dir +  " -file " + smi_dir + input_file_name + "_fingerprint" +  ".csv"
    subprocess.call(cmd, shell=True)
    #For some reason the fingerprinting stuff reorders the chemblIDs
    pIC50_df = pIC50_df.sort_values(by=["molecule_chembl_id"])
    fingerprint_df = pd.read_csv(smi_dir + input_file_name + "_fingerprint" + ".csv")
    fingerprint_df = fingerprint_df.sort_values(by=["Name"])
    pIC50_values = pIC50_df["pIC50"].tolist()
    fingerprint_df["pIC50"] = pIC50_values
    fingerprint_df.drop("Name", axis = 1, inplace=True)
    return fingerprint_df

def full_pipeline(raw_data):
    df = pd.read_csv(raw_data)
    cleaned_df = clean_df(df)
    cleaned_smiles_df = clean_smiles(cleaned_df)
    no_charge_df = remove_charges(cleaned_smiles_df)
    normalised_df = norm_value(no_charge_df)
    pIC50_df = pIC50(normalised_df)
    pIC50_df.to_csv("pIC50.csv", index = False)
    input_file_name = os.path.normpath(raw_data).split(os.path.sep)[-1].split(".")[0]
    fingerprint_df = fingerprinting(pIC50_df, input_file_name)
    fingerprint_df.to_csv("final.csv", index = False)



def clean_smiles(test_df):
    df_no_smiles = test_df.drop(columns='canonical_smiles', axis = 1)
    df = test_df
    smiles = []
    for i in df.canonical_smiles.tolist():
        cpd = str(i).split('.')
        cpd_longest = max(cpd, key=len)
        smiles.append(cpd_longest)
    smiles = pd.Series(smiles, name='canonical_smiles')
    df_no_smiles["canonical_smiles"] = smiles.values
    return df_no_smiles

def pIC50(input):
    pIC50 = []
    for i in input['standard_value_norm']:
        molar = i * (10 ** -9)  # Converts nM to M
        pIC50.append(-np.log10(molar))
    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', axis = 1)
    return x

def norm_value(input):
    norm = []
    for i in input['standard_value']:
        if i > 100000000:
            i = 100000000
        norm.append(i)
    input['standard_value_norm'] = norm
    x = input.drop('standard_value', axis = 1)
    return x

def process_fingerprint(df_pIC50, fingerprint):
    fingerprint_df = pd.read_csv(fingerprint)
    fingerprint_df = fingerprint_df.drop(["Name"], axis = 1)
    pIC50 = df_pIC50['pIC50']
    combined_df = pd.concat([fingerprint_df, pIC50], axis = 1)
    combined_df.to_csv("final.csv", index = False)

def process_ligands(ligand_dir):
    ligands = pd.read_csv(ligand_dir)
    cleaned = clean_smiles(ligands)
    cleaned.to_csv('ligands.smi', sep='\t', index=False, header=False)


def main():
    raw_data = "/home/sjb176/IRP/ML/data/pubchem_chembl_combined.csv"
    ligands = "/home/sjb176/ligand_canonical.csv"
    fingerprint = "/home/sjb176/testing_fingerprint.csv"
    full_pipeline(raw_data)
    #process_ligands(ligands)





main()