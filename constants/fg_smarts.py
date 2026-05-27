"""SMARTS patterns for functional group visualization (highlighting in SVG output)."""

FG_SMARTS: dict[str, str] = {
    "Carboxylic acid (COOH)": "C(=O)[OH]",
    "Ester":                  "C(=O)OC",
    "Ether (C-O-C)":          "COC",
    "Aromatic ring":          "c1ccccc1",
    "Amide (CONH)":           "C(=O)N",
    "Amine (NH2)":            "N",
    "Carbonyl (C=O)":         "C=O",
    "Hydroxyl (OH)":          "[OH]",
    "Halogen":                "[F,Cl,Br,I]",
    "Nitrile (CN)":           "C#N",
    "Nitro (NO2)":            "[N+](=O)[O-]",
    "Sulfonamide":            "S(=O)(=O)N",
    "Imidazole":              "c1cnc[nH]1",
    "Benzene":                "c1ccccc1",
    "Tertiary amine":         "N(C)(C)C",
}
