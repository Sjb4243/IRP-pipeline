import pubchempy as pcp
import pandas as pd
import time
pubchem_df = pd.read_csv("pubchem.csv")
pubchem_df = pubchem_df[pubchem_df["Activity Name"] == "IC50"]
pubchem_df.rename(columns={"Activity Value [uM]": "standard_value"}, inplace=True)
no_na_pubchem_df = pubchem_df[pubchem_df.standard_value.notna()]
no_na_pubchem_df['CID'] = no_na_pubchem_df['CID'].astype(str)
no_na_pubchem_df['CID'] = no_na_pubchem_df['CID'].str.replace('.0', '')
no_na_pubchem_df.to_csv("processed_pubchem.csv", index = False)
smiles_list = []
for index, cid_id in enumerate(no_na_pubchem_df["CID"]):
    if index % 4 == 0:
        time.sleep(1)
    try:
        c = pcp.Compound.from_cid(cid_id)
        smiles_list.append(c.canonical_smiles)
        print(c.canonical_smiles)
        print(index)
    except:
        print(f"{cid_id} failed to request")
        smiles_list.append(" ")

no_na_pubchem_df["canonical_smiles"] = smiles_list
selection = ["CID", "canonical_smiles", "standard_value"]
no_na_pubchem_df = no_na_pubchem_df[selection]
no_na_pubchem_df.to_csv("pubchem_smiles.csv", index = False)
