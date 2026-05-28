"""End-to-end validation: known compounds → target prediction."""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, '.')

from utils.target_predictor import predict, format_report

test_compounds = {
    # 已知 target 列在括號中
    'Testosterone':  ('CC12CCC(=O)C=C1CCC1C2CC[C@@H]1O',       'nuclear receptor (AR/GR)'),
    'Quercetin':     ('O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12', 'COX / kinase / ER'),
    'Artemisinin':   ('[C@@H]12[C@H](CC[C@@]3(C)[C@H]1CC[C@H]([C@@H]3C)OO2)C', 'antimalarial / heme'),
    'Warfarin':      ('CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O',  'VKORC1 / CYP450'),
    'Aspirin':       ('CC(=O)Oc1ccccc1C(=O)O',                  'COX-1/2'),
    'Caffeine':      ('Cn1cnc2c1c(=O)n(C)c(=O)n2C',             'adenosine receptor / PDE'),
    'Erythromycin':  ('CC[C@@H]1OC(=O)[C@H](C)[C@@H](O[C@@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H](N(C)C)[C@H]2O)[C@](C)(O)C[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@@]1(C)O', 'ribosome / mTOR'),
    'Vincristine':   ('CCC1(CC2CC(C3=C(CN(C2)C1)C4=CC=CC=C4N3)(C(=O)OC)O)C(=O)OC', 'tubulin'),
}

for name, (smi, known) in test_compounds.items():
    pred = predict(smi, top_residues=5)
    print(format_report(pred, compound_name=name))
    print(f"  [Known target: {known}]")
    print()
