import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, '.')
from utils.fg_detector import detect_smarts

compounds = {
    'Vincristine':  'CCC1(CC2CC(C3=C(CN(C2)C1)C4=CC=CC=C4N3)(C(=O)OC)O)C(=O)OC',
    'EGCG':         'O=C(O[C@@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1',
    'Berberine':    'COc1ccc2CC3=CC=C[N+]3=Cc2c1OC',
    'Piperine':     'O=C(/C=C/C=C/c1ccc2c(c1)OCO2)N1CCCCC1',
    'Artemisinin':  '[C@@H]12[C@H](CC[C@@]3(C)[C@H]1CC[C@H]([C@@H]3C)OO2)C',
    'Camptothecin': 'CC[C@]1(O)C(=O)OCC2=C1C1=NC3=CC4=CC=CC=C4C=C3C=C1C=C2',  # PubChem CID 24360
    'Taxol':        'O=C(OC1C(=C2C(CC1(OC(=O)c1ccccc1)C)=C)C(=O)CC(O)(C2(C)C)C(=O)Oc1ccccc1)C(O)c1ccccc1',
    'Caffeine':     'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
    'Capsaicin':    'COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O',
    'Resveratrol':  'Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1',
    'Colchicine':   'COc1ccc2cc([C@@H]3CC(=O)[C@@H](NC(C)=O)C3)ccc2c1OC',
    'Psilocybin':   'COP(=O)(O)Oc1c[nH]c2cccc(CCN(C)C)c12',
    'Serotonin':    'NCCc1c[nH]c2cccc(O)c12',
    'Tryptamine':   'NCCc1c[nH]c2ccccc12',
    'Melatonin':    'COc1ccc2[nH]cc(CCNC(C)=O)c2c1',
    'Physostigmine':'CNC(=O)Oc1ccc2[nH]c(C)c(CN(C)C)c2c1',
    # ── New scaffold validation compounds ─────────────────────────────────────
    # Coumarin → should hit: Coumarin, Lactone, Ketone, Phenol, Phenyl ring
    'Warfarin':     'CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O',  # 4-hydroxycoumarin core; PubChem CID 54678486
    # Chromone/Flavonoid → should hit: Chromone, Phenol, Catechol, Phenyl ring
    'Quercetin':    'O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12',
    # Steroid → should hit: Steroid, Ketone, Hydroxyl (no phenol — A ring is cyclohexanone, not phenol)
    'Testosterone': 'CC12CCC(=O)C=C1CCC1C2CC[C@@H]1O',  # PubChem CID 6013
    # Macrolide → should hit: Macrolide, Lactone, Ether, Hydroxyl, Ketone
    'Erythromycin': 'CC[C@@H]1OC(=O)[C@H](C)[C@@H](O[C@@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H](N(C)C)[C@H]2O)[C@](C)(O)C[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@@]1(C)O',
}

print("%-18s  FGs detected" % "Compound")
print("-" * 80)
for name, smi in compounds.items():
    fgs = detect_smarts(smi)
    if not fgs:
        print("%-18s  [NONE]" % name)
    else:
        print("%-18s  %s" % (name, fgs))
