"""
Mouse liver scRNA-seq / snRNA-seq annotation markers (mouse gene symbols).

Structure:  PANEL -> POPULATION -> {"genes": [...], "ref": "...", "note": "..."}
Panels are designed to be plotted separately (one dotplot each), coarse -> fine:

  L1_overview      : all major compartments (run on the full object)
  L2_mesenchyme    : HSC / fibroblast / activated / mural / mesothelial (subset on Pdgfrb+ etc.)
  L2_macrophage    : KC states, moKC, LAM, capsule, monocytes (subset on Adgre1/Csf1r/Fcgr1+)
  L2_lymphoid_DC   : optional finer immune panel
  QC               : stress / ambient / doublet flags

Annotation logic (boolean rules, read alongside the dotplots):

  Mesenchyme
    HSC          = Pdgfrb+ & (Lrat | Reln | Ecm1 | Colec11)+ & Myh11-
      PaHSC      = HSC & Ngfr hi & Itgb3 hi                    (portal)
      CaHSC      = HSC & Adamtsl2 hi & Rspo3 hi                (central)
      cyHSC      = HSC & Col1a1 lo & (Hgf, Rspo3, Cxcl12) hi   (low-activation)
      myHSC/aHSC = HSC identity retained (Lrat/Des) & Col1a1/Col3a1/Acta2/Timp1 hi
    Fibroblast   = Pdgfra+ & (Cd34 | Ly6a | Pi16 | Dpt | Clec3b)+ & Lrat- & Reln-
      PF         = Fibroblast & (Thy1, Eln, Fbln2) & low collagen program
      aPF        = Fibroblast & (Msln, Muc16, Gpc3, Aspn, Bnc1) & Col1a1 hi   (cholestatic)
    Mesothelium  = Msln+ & Upk3b+ & Lrrn4+ & Col1a1 lo          (do NOT call aPF if Upk3b+)
    VSMC         = Myh11+ & Cnn1+

  Macrophages
    pan-Mac      = Adgre1 / Csf1r / Fcgr1 / C1qa +
    emKC         = Clec4f & Vsig4 & Cd5l & Timd4 +
    moKC         = Clec4f & Vsig4 & Clec1b + ,  Timd4 lo/-
    LAM          = Clec4f- & (Trem2, Gpnmb, Spp1, Cd9, Cd63) +
    LAM-like KC  = Clec4f+ & Trem2/Gpnmb/Cd9 +   (KC adopting LAM program)
    Capsule Mac  = Clec4f- & Timd4- & Cx3cr1+    (rarely captured by perfusion)
    Mo-Mac trans = Ccr2 + Cx3cr1 + MHCII hi, KC & LAM genes low

All gene symbols should be verified against adata.var_names (filter_to_var does this).
"""

MARKERS = {
    # ------------------------------------------------------------------ L1
    "L1_overview": {
        "Hepatocyte":            {"genes": ["Alb", "Apoa2", "Ttr", "Fabp1"],
                                  "ref": "general", "note": "Alb/Apoa2 are major ambient RNA; judge by fraction, not presence"},
        "Hep_periportal":        {"genes": ["Cyp2f2", "Hal", "Sds", "Pck1"],
                                  "ref": "Halpern 2017 Nature"},
        "Hep_pericentral":       {"genes": ["Glul", "Cyp2e1", "Oat", "Cyp1a2"],
                                  "ref": "Halpern 2017 Nature"},
        "Cholangiocyte":         {"genes": ["Krt19", "Epcam", "Sox9", "Krt7"],
                                  "ref": "general", "note": "Spp1 also high here; not used to avoid LAM overlap"},
        "LSEC":                  {"genes": ["Kdr", "Stab2", "Clec4g", "Fcgr2b", "Dnase1l3"],
                                  "ref": "Halpern 2018 Nat Biotech; Guilliams 2022 Cell"},
        "LSEC_periportal":       {"genes": ["Efnb2", "Ltbp4", "Msr1", "Dll4"],
                                  "ref": "Halpern 2018 Nat Biotech"},
        "LSEC_pericentral_CVEC": {"genes": ["Kit", "Rspo3", "Wnt2", "Wnt9b", "Thbd"],
                                  "ref": "Halpern 2018 Nat Biotech", "note": "Rspo3 also CaHSC/cyHSC"},
        "Arterial_EC":           {"genes": ["Gja5", "Sox17", "Hey1", "Sema3g"],
                                  "ref": "general"},
        "Lymphatic_EC":          {"genes": ["Prox1", "Mmrn1", "Ccl21a", "Pdpn", "Flt4"],
                                  "ref": "general"},
        "Mesenchyme_pan":        {"genes": ["Pdgfrb", "Col1a2", "Dcn", "Pdgfra"],
                                  "ref": "Dobie 2019 Cell Rep"},
        "Immune_pan":            {"genes": ["Ptprc"], "ref": "general"},
        "Macrophage_pan":        {"genes": ["Adgre1", "Csf1r", "Fcgr1", "C1qa", "Cd68"],
                                  "ref": "general"},
        "Monocyte":              {"genes": ["Ly6c2", "Ccr2", "Plac8", "Chil3"],
                                  "ref": "general"},
        "cDC1":                  {"genes": ["Xcr1", "Clec9a", "Cadm1"], "ref": "general"},
        "cDC2":                  {"genes": ["Cd209a", "Mgl2", "Clec10a", "Sirpa"], "ref": "general"},
        "pDC":                   {"genes": ["Siglech", "Bst2", "Ccr9", "Ly6d"], "ref": "general"},
        "Neutrophil":            {"genes": ["S100a8", "S100a9", "Retnlg", "Csf3r"], "ref": "general"},
        "Mast_Basophil":         {"genes": ["Cpa3", "Mcpt8", "Prss34", "Kit"],
                                  "ref": "general", "note": "Cpa3 mast; Mcpt8/Prss34 basophil"},
        "B_cell":                {"genes": ["Cd79a", "Ms4a1", "Cd19"], "ref": "general"},
        "Plasma_cell":           {"genes": ["Jchain", "Mzb1", "Xbp1"], "ref": "general"},
        "T_cell":                {"genes": ["Cd3e", "Cd4", "Cd8a", "Trdc"], "ref": "general"},
        "NK_ILC1":               {"genes": ["Ncr1", "Klrb1c", "Gzma"], "ref": "general"},
        "Proliferating":         {"genes": ["Mki67", "Top2a"], "ref": "general"},
    },

    # ------------------------------------------------------------------ L2 mesenchyme
    "L2_mesenchyme": {
        "Mesenchyme_pan":        {"genes": ["Pdgfrb", "Col1a2", "Dcn"],
                                  "ref": "Dobie 2019 Cell Rep"},
        # --- HSC identity
        "HSC_core":              {"genes": ["Lrat", "Reln", "Ecm1", "Colec11", "Gdf2", "Vipr1", "Des"],
                                  "ref": "Dobie 2019; Guilliams 2022 (Gdf2/BMP9 from HSC)"},
        "HSC_portal_PaHSC":      {"genes": ["Ngfr", "Itgb3"],
                                  "ref": "Dobie 2019 Cell Rep"},
        "HSC_central_CaHSC":     {"genes": ["Adamtsl2", "Rspo3", "Lpar1"],
                                  "ref": "Dobie 2019 Cell Rep",
                                  "note": "Lpar1 = therapeutic target on collagen-producing CaHSC; enrichment, not exclusive"},
        "HSC_cytokine_cyHSC":    {"genes": ["Hgf", "Rspo3", "Cxcl12"],
                                  "ref": "Filliol 2022 Nature",
                                  "note": "Relative state: requires Col1a1 lo. Overlaps quiescent HSC"},
        "HSC_myofibro_myHSC":    {"genes": ["Col1a1", "Col3a1", "Acta2", "Timp1", "Lox", "Loxl2", "Tgfbi"],
                                  "ref": "Filliol 2022 Nature; Dobie 2019",
                                  "note": "Call aHSC only if HSC_core genes retained (esp. Lrat/Des)"},
        "Senescent_mesench":     {"genes": ["Cdkn1a", "Cdkn2a", "Plaur"],
                                  "ref": "review (Semin Liver Dis 2025)", "note": "state; optional"},
        # --- Fibroblasts
        "Fibro_universal":       {"genes": ["Pi16", "Dpt", "Ly6a", "Cd34", "Clec3b"],
                                  "ref": "Buechler 2021 Nature; Dobie 2019 (CD34+ portal fibroblasts)"},
        "Portal_fibroblast":     {"genes": ["Thy1", "Eln", "Fbln2", "Gli1", "Mfap4"],
                                  "ref": "Katsumata 2017 Hepatology (Thy1); Gli1 periportal mesenchyme; Dobie 2019 (Mfap4)"},
        "Activated_PF_aPF":      {"genes": ["Msln", "Muc16", "Gpc3", "Aspn", "Bnc1", "Fbln2"],
                                  "ref": "Iwaisako 2014 PNAS; Koyama 2017 JCI; Nishio 2021 PNAS",
                                  "note": "Cholestatic models (BDL, DDC, Mdr2-/-) only; not induced by CCl4. Exclude Upk3b+ (mesothelium)"},
        "Activated_ECM_shared":  {"genes": ["Col1a1", "Col3a1", "Postn", "Cthrc1", "Fn1", "Spp1"],
                                  "ref": "general; Tsukui 2020 Nat Commun (Cthrc1, cross-organ)",
                                  "note": "Use as activation SCORE across HSC and PF, not as identity"},
        # --- other mesenchyme
        "Mesothelial_capsule":   {"genes": ["Upk3b", "Upk1b", "Lrrn4", "Msln", "Wt1"],
                                  "ref": "Bhattacharjee 2021 JCI (PF/mesothelial markers)"},
        "VSMC":                  {"genes": ["Myh11", "Cnn1", "Tagln", "Acta2"],
                                  "ref": "Dobie 2019 (Calponin1 in vessel walls)"},
        "Proliferating":         {"genes": ["Mki67", "Top2a"], "ref": "general"},
    },

    # ------------------------------------------------------------------ L2 macrophages
    "L2_macrophage": {
        "Macrophage_pan":        {"genes": ["Adgre1", "Csf1r", "Fcgr1", "C1qa"],
                                  "ref": "general"},
        "KC_core":               {"genes": ["Clec4f", "Vsig4", "Cd5l", "Folr2", "Marco", "Clec1b"],
                                  "ref": "Guilliams 2022 Cell (VSIG4, FOLR2); Scott 2016"},
        "KC_identity_TF":        {"genes": ["Id3", "Spic", "Nr1h3"],
                                  "ref": "Bonnardel 2019 Immunity"},
        "KC_embryonic_Timd4":    {"genes": ["Timd4"],
                                  "ref": "Scott 2016; Bonnardel 2019; Tran 2020 Immunity",
                                  "note": "emKC = KC_core+ Timd4+ ; moKC = KC_core+ Timd4 lo/-"},
        "KC2_candidate":         {"genes": ["Cd36", "Esam", "Mrc1"],
                                  "ref": "Bleriot 2021 Immunity",
                                  "note": "Guilliams 2022 suggests KC2 may be KC-LSEC doublets: check Kdr/Stab2 co-expression"},
        "LAM_Trem2":             {"genes": ["Trem2", "Gpnmb", "Spp1", "Cd9", "Cd63", "Fabp5", "Lgals3"],
                                  "ref": "Remmerie 2020; Seidman 2020; Xiong 2019; Guilliams 2022; Ganguly 2024 PNAS",
                                  "note": "LAM = KC_core- ; LAM-like KC = KC_core+ with this program"},
        "MoMac_transitional":    {"genes": ["Ccr2", "Cx3cr1", "H2-Aa", "Cd74", "Itgax"],
                                  "ref": "Daemen 2021 Cell Rep (C-LAM / Cx3cr1-Ccr2 MdM)"},
        "Capsule_macrophage":    {"genes": ["Cx3cr1", "Cd207", "Itgax"],
                                  "ref": "Sierro 2017 Immunity",
                                  "note": "Clec4f-/Timd4-; Cd207 also on KCs, Cx3cr1 discriminates"},
        "Monocyte_Ly6Chi":       {"genes": ["Ly6c2", "Ccr2", "Plac8", "Chil3", "F13a1"],
                                  "ref": "general"},
        "Monocyte_Ly6Clo":       {"genes": ["Nr4a1", "Ace", "Treml4", "Itgal"],
                                  "ref": "general"},
        "IFN_state":             {"genes": ["Isg15", "Ifit3", "Irf7", "Rsad2"],
                                  "ref": "general", "note": "state, overlay on any macrophage type"},
        "Inflammatory_state":    {"genes": ["Il1b", "Tnf", "Cxcl2", "Nfkbia"],
                                  "ref": "general", "note": "also induced by enzymatic digestion; compare with QC panel"},
        "Proliferating":         {"genes": ["Mki67", "Top2a"], "ref": "general"},
        "LSEC_doublet_check":    {"genes": ["Kdr", "Stab2"], "ref": "Guilliams 2022 Cell"},
    },

    # ------------------------------------------------------------------ L2 lymphoid/DC (optional)
    "L2_lymphoid_DC": {
        "CD4_T":                 {"genes": ["Cd3e", "Cd4", "Il7r"], "ref": "general"},
        "Treg":                  {"genes": ["Foxp3", "Il2ra", "Ctla4"], "ref": "general"},
        "CD8_T":                 {"genes": ["Cd8a", "Cd8b1", "Gzmk"], "ref": "general"},
        "gdT":                   {"genes": ["Trdc", "Tcrg-C1", "Il17a"], "ref": "general"},
        "NKT_MAIT":              {"genes": ["Zbtb16", "Klrb1c", "Cd3e"], "ref": "general"},
        "NK":                    {"genes": ["Ncr1", "Eomes", "Itga2", "Gzma"], "ref": "general"},
        "ILC1":                  {"genes": ["Ncr1", "Itga1", "Cd200r1", "Tnfsf10"],
                                  "ref": "general", "note": "liver-resident; Eomes- Itga2-"},
        "ILC2":                  {"genes": ["Gata3", "Il1rl1", "Areg"], "ref": "general"},
        "B_cell":                {"genes": ["Cd79a", "Ms4a1", "Cd19"], "ref": "general"},
        "Plasma_cell":           {"genes": ["Jchain", "Mzb1", "Xbp1"], "ref": "general"},
        "cDC1":                  {"genes": ["Xcr1", "Clec9a", "Cadm1"], "ref": "general"},
        "cDC2":                  {"genes": ["Cd209a", "Mgl2", "Clec10a"], "ref": "general"},
        "mregDC":                {"genes": ["Ccr7", "Fscn1", "Cacnb3"], "ref": "general"},
        "pDC":                   {"genes": ["Siglech", "Bst2", "Ccr9"], "ref": "general"},
    },

    # ------------------------------------------------------------------ QC
    "QC": {
        "Digestion_stress":      {"genes": ["Fos", "Jun", "Egr1", "Hspa1a"], "ref": "Guilliams 2022 (digestion signature)"},
        "Ambient_hepatocyte":    {"genes": ["Alb", "Apoa1", "Apoc3"], "ref": "general"},
        "Erythroid":             {"genes": ["Hba-a1", "Hbb-bs"], "ref": "general"},
        "Mito":                  {"genes": ["mt-Co1", "mt-Nd1"], "ref": "general"},
    },
}


# ------------------------------------------------------------------ helpers
def get_dotplot_dict(panel):
    """{population: [genes]} for one panel, e.g. sc.pl.dotplot(var_names=...)."""
    return {pop: d["genes"] for pop, d in MARKERS[panel].items()}


def filter_to_var(marker_dict, var_names, verbose=True):
    """Drop genes absent from the object (avoids dotplot KeyErrors); report what was dropped."""
    var = set(var_names)
    out, missing = {}, {}
    for pop, genes in marker_dict.items():
        keep = [g for g in genes if g in var]
        miss = [g for g in genes if g not in var]
        if keep:
            out[pop] = keep
        if miss:
            missing[pop] = miss
    if verbose and missing:
        print("Missing genes:", missing)
    return out


def to_long_table():
    """Long format: panel, population, gene, ref, note (for R/Seurat or a CSV)."""
    rows = []
    for panel, pops in MARKERS.items():
        for pop, d in pops.items():
            for g in d["genes"]:
                rows.append({"panel": panel, "population": pop, "gene": g,
                             "ref": d.get("ref", ""), "note": d.get("note", "")})
    return rows


# ------------------------------------------------------------------ usage
if __name__ == "__main__":
    import csv
    rows = to_long_table()
    with open("mouse_liver_markers.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["panel", "population", "gene", "ref", "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rows written")

    # Scanpy example:
    # import scanpy as sc
    # m = filter_to_var(get_dotplot_dict("L2_mesenchyme"), adata.var_names)
    # sc.pl.dotplot(adata_mes, var_names=m, groupby="leiden", standard_scale="var",
    #               dendrogram=False, swap_axes=False)
    #
    # Seurat example (R), reading the CSV:
    # mk <- read.csv("mouse_liver_markers.csv")
    # g  <- unique(mk$gene[mk$panel == "L2_macrophage"])
    # DotPlot(obj_mac, features = intersect(g, rownames(obj_mac))) + RotatedAxis()
