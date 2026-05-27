"""
explore_pubchem.py
先探索PubChem API能給我們什麼，再決定要抓哪些欄位
"""

import requests
import json

BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def explore_compound(smiles: str, name: str = "compound"):
    """
    給一個SMILES，看PubChem能回傳什麼資訊
    """
    print(f"\n{'='*60}")
    print(f"探索: {name} | SMILES: {smiles}")
    print(f"{'='*60}")

    # ── Step 1: 用SMILES換取CID ──────────────────────────
    url = f"{BASE_URL}/compound/smiles/{requests.utils.quote(smiles)}/cids/JSON"
    r = requests.get(url)
    
    if r.status_code != 200:
        print(f"✗ 找不到這個compound，狀態碼：{r.status_code}")
        return None
    
    cid = r.json()["IdentifierList"]["CID"][0]
    print(f"\n✓ PubChem CID: {cid}")

    # ── Step 2: 抓這個CID的properties ───────────────────
    props = "MolecularFormula,MolecularWeight,IUPACName,InChIKey"
    url2  = f"{BASE_URL}/compound/cid/{cid}/property/{props}/JSON"
    r2    = requests.get(url2)
    props_data = r2.json()["PropertyTable"]["Properties"][0]
    
    print(f"\n── Properties ──")
    for k, v in props_data.items():
        print(f"  {k}: {v}")

    # ── Step 3: 抓Classification（這裡會有官能基資訊）───
    url3 = f"{BASE_URL}/compound/cid/{cid}/classification/JSON"
    r3   = requests.get(url3)
    
    if r3.status_code == 200:
        data = r3.json()
        # 只印前3層，避免資料太多
        print(f"\n── Classification (前5筆) ──")
        hierarchies = data.get("Hierarchies", {}).get("Hierarchy", [])
        for h in hierarchies[:5]:
            print(f"  Source: {h.get('SourceName', 'N/A')}")
            node = h.get("Node", [{}])[0]
            print(f"  Name:   {node.get('Information', {}).get('Name', 'N/A')}")
            print()
    
    # ── Step 4: 抓Computed Properties（包含官能基數量）──
    url4 = f"{BASE_URL}/compound/cid/{cid}/JSON"
    r4   = requests.get(url4)
    
    if r4.status_code == 200:
        full = r4.json()["PC_Compounds"][0]
        print(f"── Computed Properties ──")
        for prop in full.get("props", []):
            label = prop.get("urn", {}).get("label", "")
            name_p = prop.get("urn", {}).get("name", "")
            value = prop.get("value", {})
            # 只印跟官能基相關的
            if any(kw in label.lower() for kw in 
                   ["functional", "group", "bond", "ring", "stereo", "count"]):
                val = list(value.values())[0] if value else "N/A"
                print(f"  {label} / {name_p}: {val}")

    # ── Step 5: 把完整raw data存下來讓你自己看 ──────────
    output_file = f"data/{name}_pubchem_raw.json"
    with open(output_file, "w") as f:
        json.dump(r4.json(), f, indent=2)
    print(f"\n✓ 完整raw data已存到: {output_file}")
    
    return cid


if __name__ == "__main__":
    # 先用三個有代表性的compound測試
    # 各自有不同的官能基組合
    test_compounds = {
        "Aspirin":       "CC(=O)Oc1ccccc1C(=O)O",
        "Caffeine":      "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "Amitriptyline": "CN(C)CCCN1c2ccccc2CCc2ccccc21",
    }

    for name, smiles in test_compounds.items():
        explore_compound(smiles, name)

    print(f"\n{'='*60}")
    print("探索完成。去 data/ 資料夾看 _pubchem_raw.json 的結構")
    print("告訴我哪些欄位對你來說有感覺，我們再決定要抓什麼")