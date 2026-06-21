#!/usr/bin/env python3
"""Build per-BioProject DESeq2 sample sheets from ENA metadata.

The output is deliberately conservative: every project gets a full metadata sheet,
and projects with a plausible class field also get a DESeq2-ready Run/condition CSV.
The candidate field and review status are written to a summary table.
"""

from __future__ import annotations

import argparse
import csv
import re
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd


ANALYSIS = Path("/pc2/users/o/omiks001/hpc-prf-omiks/ja/analysis")
OUT_DIR = ANALYSIS / "deseq2_metadata"
FULL_DIR = OUT_DIR / "full_metadata"
SHEET_DIR = OUT_DIR / "sample_sheets"
REVIEW_DIR = OUT_DIR / "review"
MANUAL_NOT_FOR_DESEQ2 = {"PRJNA480287"}

ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
RATE_LIMIT_SECONDS = 0.25

FIELDS = [
    "run_accession",
    "sample_accession",
    "secondary_sample_accession",
    "sample_alias",
    "sample_title",
    "sample_description",
    "experiment_title",
    "experiment_alias",
    "study_title",
    "study_alias",
    "library_strategy",
    "library_source",
    "instrument_model",
    "scientific_name",
    "sex",
    "host_sex",
    "submitted_host_sex",
    "age",
    "dev_stage",
    "disease",
    "host_status",
    "host_phenotype",
    "host_genotype",
    "cell_line",
    "cell_type",
    "tissue_type",
    "tissue_lib",
    "isolation_source",
    "strain",
    "sub_strain",
    "experimental_factor",
    "control_experiment",
    "description",
]

CLASS_FIELD_PRIORITY = [
    "control_experiment",
    "experimental_factor",
    "disease",
    "host_status",
    "host_phenotype",
    "host_genotype",
    "cell_type",
    "cell_line",
    "tissue_type",
    "isolation_source",
    "strain",
    "sub_strain",
    "sample_title",
    "sample_description",
]

BAD_VALUES = {
    "",
    "na",
    "n/a",
    "not applicable",
    "not available",
    "not collected",
    "missing",
    "unknown",
    "none",
    "null",
    "-",
}

CONTROL_PATTERNS = [
    r"\bcontrol\b",
    r"\bctrl\b",
    r"\bvehicle\b",
    r"\buntreated\b",
    r"\bmock\b",
    r"\bhealthy\b",
    r"\bnormal\b",
    r"\bwt\b",
    r"\bwild[ -]?type\b",
    r"\bbaseline\b",
    r"\bday[ _-]?0\b",
    r"\b0h\b",
]
TREATMENT_PATTERNS = [
    r"\btreated\b",
    r"\btreatment\b",
    r"\binfected\b",
    r"\bdisease\b",
    r"\btumou?r\b",
    r"\bcancer\b",
    r"\bcase\b",
    r"\bko\b",
    r"\bknock[ -]?out\b",
    r"\bmutant\b",
    r"\bstimulated\b",
]


def manual_project_condition(project_id: str, row: pd.Series) -> tuple[str, str] | None:
    """Return (condition_raw, condition) for manually curated project rules."""
    if project_id == "PRJNA1014106":
        return "ataxia PBMC", "ataxia_pbmc"

    if project_id == "PRJNA1014965":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        if alias.startswith("MEC_N_"):
            return "Normal", "normal"
        if alias.startswith("MEC_T_"):
            return "Tumor", "tumor"

    if project_id == "PRJNA1105191":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        group = re.sub(r"_R\d+$", "", alias)
        group_map = {
            "SAEC": ("SAEC", "saec"),
            "IPSC_C6": ("IPSC_C6", "ipsc_c6"),
            "IPSC_C7": ("IPSC_C7", "ipsc_c7"),
            "H1417": ("H1417", "h1417"),
            "H69": ("H69", "h69"),
            "H841": ("H841", "h841"),
            "DMS53": ("DMS53", "dms53"),
            "DMS53_siCtrl": ("DMS53 siCtrl", "dms53_sictrl"),
            "DMS53_siNFIC": ("DMS53 siNFIC", "dms53_sinfic"),
            "H841_siCtrl": ("H841 siCtrl", "h841_sictrl"),
            "H841_siNFIC": ("H841 siNFIC", "h841_sinfic"),
        }
        if group in group_map:
            return group_map[group]

    if project_id == "PRJNA1120369":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_title", ""))
        if alias.startswith("Mega_01"):
            return "Mega_01", "mega_01"
        if alias.startswith("Mega_02"):
            return "Mega_02", "mega_02"
        if alias.startswith("WT"):
            return "WT", "wt"

    if project_id == "PRJNA1127555":
        title = clean_value(row.get("sample_title", ""))
        if title == "Naive":
            return "Naive", "naive"
        if title == "Engager":
            return "Engager", "engager"
        if title == "Tumor":
            return "Tumor", "tumor"

    if project_id == "PRJNA1133701":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        group = re.sub(r"_\d+$", "", alias)
        group = group.replace("Etopoisde", "Etoposide")
        tp53_match = re.fullmatch(r"(TP53WT|TP53KO|TP53rescue) (control|Decitabine|Etoposide|D_E)", group)
        if tp53_match:
            genotype, treatment = tp53_match.groups()
            raw = f"{genotype} {treatment}"
            return raw, normalize_condition(raw)

    if project_id == "PRJNA1165739":
        cell_line = clean_value(row.get("cell_line", ""))
        cell_map = {
            "Kasumi-1": ("Kasumi-1 / CRL-2724", "kasumi1_crl2724"),
            "HL-60": ("HL-60 / CCL-240", "hl60_ccl240"),
            "THP-1": ("THP-1 / TIB-202", "thp1_tib202"),
            "K562": ("K562 / CCL-243", "k562_ccl243"),
        }
        if cell_line in cell_map:
            return cell_map[cell_line]

    if project_id == "PRJNA1175639":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        if alias.startswith("hBA-d0"):
            return "Brown preadipocytes Day 0", "brown_preadipocyte_d0"
        if alias.startswith("hBA-diff"):
            return "Brown differentiated", "brown_differentiated"
        if alias.startswith("hWA-d0"):
            return "White preadipocytes Day 0", "white_preadipocyte_d0"
        if alias.startswith("hWA-diff"):
            return "White differentiated", "white_differentiated"

    if project_id == "PRJNA1176539":
        title = clean_value(row.get("sample_title", "")).lower()
        ratio_match = re.search(r"\b0\.\d+\b", title)
        if ratio_match:
            ratio = ratio_match.group(0).replace(".", "p")
            return f"Raji/Ramos spike-in ratio {ratio_match.group(0)}", f"raji_ramos_ratio_{ratio}"
        benchmark_map = [
            ("raji cells", "Raji cells", "raji_cells"),
            ("ramos cells", "Ramos cells", "ramos_cells"),
            ("mrc-5 cells", "MRC-5 cells", "mrc5_cells"),
            ("medium", "cell culture medium", "cell_culture_medium"),
            ("pbs", "NTC / PBS", "ntc_pbs"),
            ("water", "water", "water"),
        ]
        for needle, raw, condition in benchmark_map:
            if needle in title:
                return raw, condition

    if project_id == "PRJNA1185243":
        alias = clean_value(row.get("experiment_alias", ""))
        subject_match = re.search(r"WU(\d+-\d+)_d57", alias)
        if subject_match:
            subject = subject_match.group(1).replace("-", "_")
            return f"WU{subject_match.group(1)} day 57 lymph node", f"wu{subject}_d57_lymph_node"

    if project_id == "PRJNA316201":
        return "T47D breast cancer cell line", "t47d_breast_cancer_cell_line"

    text = " ".join(
        clean_value(row.get(col, ""))
        for col in ["experiment_title", "experiment_alias", "sample_title", "sample_description", "tissue_type"]
    ).lower()

    if project_id == "PRJNA321028":
        stage_map = [
            ("stage iii", "Stage III Endometrial Cancer", "stage_iii_endometrial_cancer"),
            ("stage ii", "Stage II Endometrial Cancer", "stage_ii_endometrial_cancer"),
            ("stage i", "Stage I Endometrial Cancer", "stage_i_endometrial_cancer"),
        ]
        for needle, raw, condition in stage_map:
            if needle in text:
                return raw, condition

    if project_id == "PRJNA321087":
        # Cancer must be checked before generic epithelium: some cancer samples
        # have tissue_type=Epithelium but experiment titles identify "cancer".
        histology_map = [
            ("cancer", "cancer", "cancer"),
            ("gland", "gland", "gland"),
            ("epithelium", "epithelium", "epithelium"),
            ("muscle", "muscle", "muscle"),
        ]
        for needle, raw, condition in histology_map:
            if re.search(rf"\b{needle}\b", text):
                return raw, condition

    if project_id == "PRJNA321967":
        germ_cell_map = [
            ("sertoli", "Sertoli cells", "sertoli_cells"),
            ("primary_spermatocytes", "primary spermatocytes", "primary_spermatocytes"),
            ("primary spermatocytes", "primary spermatocytes", "primary_spermatocytes"),
            ("spermatids", "spermatids", "spermatids"),
            ("undifferentiated_spermatogonia", "undifferentiated spermatogonia", "undifferentiated_spermatogonia"),
            ("undifferentiated spermatogonia", "undifferentiated spermatogonia", "undifferentiated_spermatogonia"),
        ]
        for needle, raw, condition in germ_cell_map:
            if needle in text:
                return raw, condition

    if project_id == "PRJNA352875":
        if "rho 0" in text and "tnf" in text:
            return "rho0 + TNF", "rho0_tnf"
        if "without mitochondrial dna" in text and "tnf" in text:
            return "rho0 + TNF", "rho0_tnf"
        if "rho 0" in text or "without mitochondrial dna" in text:
            return "rho0 / without mitochondrial DNA", "rho0"
        if "tnf" in text:
            return "TNF-treated control", "tnf_treated_control"
        if re.search(r"\bcontrol\b|\bcon replicate\b", text):
            return "control", "control"

    if project_id == "PRJNA373978":
        if "leptotene-zygotene" in text or "leptotene zygotene" in text:
            return "leptotene-zygotene spermatocytes", "leptotene_zygotene_spermatocytes"
        if "pachytene" in text:
            return "pachytene spermatocytes", "pachytene_spermatocytes"

    if project_id == "PRJNA378952":
        if "wdlps" in text or "well-differentiated liposarcoma" in text:
            return "WDLPS", "wdlps"
        if "lsc" in text or "lung sarcomatoid carcinoma" in text:
            return "LSC", "lsc"

    if project_id == "PRJNA381115":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        alias = re.sub(r"^mRNA[_-]", "", alias, flags=re.IGNORECASE)
        alias = alias.replace("mRNA_", "")
        if alias:
            return alias, normalize_condition(alias)

    if project_id == "PRJNA381757":
        if "primed pluripotent" in text:
            return "primed pluripotent stem cells", "primed_pluripotent_stem_cells"
        if "naive induced pluripotent" in text:
            return "naive induced pluripotent stem cells", "naive_induced_pluripotent_stem_cells"
        if "fibroblast" in text:
            return "fibroblast", "fibroblast"

    if project_id == "PRJNA384289":
        if "primary tumour" in text or "primary colorectal tumour" in text:
            return "primary colorectal tumour", "primary_colorectal_tumour"
        if "circulating tumour cells" in text or re.search(r"sample\d+c", text):
            return "blood / circulating tumour cells", "circulating_tumour_cells"

    if project_id == "PRJNA386992":
        if "remission" in text:
            return "remission", "remission"
        if "aml case" in text or "aml-dmin" in text:
            return "AML case", "aml_case"

    if project_id == "PRJNA397941":
        # Human naive pluripotency culture/reprogramming states. Keep compound
        # conditions distinct because they represent different protocols.
        if "no klf4" in text:
            if "d6" in text:
                return "E8 primed D6 no Klf4", "e8_primed_d6_no_klf4"
            if "d2" in text:
                return "E8 primed D2 no Klf4", "e8_primed_d2_no_klf4"
            return "E8 primed no Klf4", "e8_primed_no_klf4"
        if "with klf4" in text:
            if "d6" in text:
                return "E8 primed D6 with Klf4", "e8_primed_d6_with_klf4"
            if "d2" in text:
                return "E8 primed D2 with Klf4", "e8_primed_d2_with_klf4"
            return "E8 primed with Klf4", "e8_primed_with_klf4"
        if "mrna primed" in text:
            return "mRNA primed", "mrna_primed"
        if "mrna rt2ilgoy" in text:
            return "mRNA rt2iLGoY", "mrna_rt2ilgoy"
        if "ct2ilgoy oskm" in text:
            return "ct2iLGoY OSKM", "ct2ilgoy_oskm"
        if "ct2ilgoy k" in text:
            return "ct2iLGoY K", "ct2ilgoy_k"
        if "rt2ilgoy" in text or "rt2ilgöy" in text:
            return "rt2iLGoY", "rt2ilgoy"
        if "rrset" in text or "rrset" in text.replace("rset", "rset"):
            return "rRSeT", "rrset"
        if "rnhsm" in text:
            return "rNHSM", "rnhsm"
        if "r5ilaf" in text:
            return "r5iLAF", "r5ilaf"
        if "crset" in text:
            return "cRSeT", "crset"
        if "cnhsm" in text:
            return "cNHSM", "cnhsm"
        if "e8 primed" in text or re.search(r"\bprimed\b", text):
            return "E8 primed", "e8_primed"

    if project_id == "PRJNA416439":
        if "adjacent normal liver" in text:
            return "adjacent normal liver tissue", "adjacent_normal_liver_tissue"
        if "hepatoblastoma" in text:
            return "hepatoblastoma tissue", "hepatoblastoma"

    if project_id == "PRJNA419934":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        alias = re.sub(r"_(?:sc)?rnaseq$", "", alias, flags=re.IGNORECASE)
        alias_lc = alias.lower()
        if alias_lc in {"parental", "caov3"}:
            return "parental / CAOV3", "parental_caov3"
        if re.fullmatch(r"s0[1-4]", alias_lc):
            return f"sensitive {alias.upper()}", f"sensitive_{alias_lc}"
        resistant = re.fullmatch(r"(r\d+)_step(\d+)(a?)", alias_lc)
        if resistant:
            clone, step, suffix = resistant.groups()
            return f"{clone.upper()} step{step}{suffix}", f"{clone}_step{step}{suffix}"

    if project_id == "PRJNA451395":
        if "a2027g" in text:
            return "SMC1A c.A2027G", "smc1a_c_a2027g"
        if "wild-type" in text or "wild type" in text or "smc1awt" in text:
            return "SMC1A wild-type", "smc1a_wild_type"
        if "control" in text:
            return "Control", "control"

    if project_id == "PRJNA479479":
        is_p1 = "p1" in text or "early passage" in text
        is_p4 = "p4" in text or "late passage" in text
        is_tcp = "tcp" in text or "tissue culture plastic" in text
        is_scm = "scm" in text or "synoviocyte matrix" in text
        if is_p1 and is_tcp:
            return "P1 TCP", "p1_tcp"
        if is_p4 and is_tcp:
            return "P4 TCP", "p4_tcp"
        if is_p1 and is_scm:
            return "P1 SCM", "p1_scm"
        if is_p4 and is_scm:
            return "P4 SCM", "p4_scm"

    if project_id == "PRJNA480287":
        return "not for DESeq2", "not_for_deseq2"

    if project_id == "PRJNA509121":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        age_raw = "young" if alias.upper().startswith("Y") else "older"
        if "postCON" in alias:
            return f"{age_raw} 5h post concentric exercise", f"{age_raw}_5h_post_concentric"
        if "postECC" in alias:
            return f"{age_raw} 5h post eccentric exercise", f"{age_raw}_5h_post_eccentric"
        if "_BL" in alias or alias.endswith("BL"):
            return f"{age_raw} baseline", f"{age_raw}_baseline"

    if project_id == "PRJNA528522":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        prefix = alias.split("_", 1)[0].upper()
        hepatic_map = {
            "PSC": ("human pluripotent stem cell", "pluripotent_stem_cell"),
            "HE": ("hPSC hepatic endoderm", "hepatic_endoderm"),
            "MH": ("hPSC mature hepatocytes", "mature_hepatocytes"),
            "EM": ("organoid expansion medium", "organoid_expansion_medium"),
            "DM": ("organoid differentiation medium", "organoid_differentiation_medium"),
            "HM": ("organoid hepatic medium", "organoid_hepatic_medium"),
        }
        if prefix in hepatic_map:
            return hepatic_map[prefix]

    if project_id == "PRJNA544334":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        if alias.startswith("XHS"):
            leukapheresis_map = {
                "EMP": ("leukapheresis EMP", "leukapheresis_emp"),
                "LMPP": ("leukapheresis LMPP", "leukapheresis_lmpp"),
                "90neg": ("leukapheresis CD90-negative MPP", "leukapheresis_cd90neg_mpp"),
                "90pos": ("leukapheresis CD90-positive HSC", "leukapheresis_cd90pos_hsc"),
                "CD34P": ("bone marrow CD34+ HSPC", "bone_marrow_cd34p_hspc"),
                "CD133P": ("bone marrow CD133+ HSC", "bone_marrow_cd133p_hsc"),
                "CD38L": ("bone marrow CD38low HSC", "bone_marrow_cd38low_hsc"),
                "CD90P": ("bone marrow CD90+ HSC", "bone_marrow_cd90p_hsc"),
            }
            suffix = alias.rsplit("_", 1)[-1]
            if suffix in leukapheresis_map:
                return leukapheresis_map[suffix]
        parts = alias.split("_")
        if len(parts) >= 3:
            subtype = parts[-1].lower()
            return f"bone marrow {subtype.upper()}", f"bone_marrow_{subtype}"

    if project_id == "PRJNA610985":
        if "siartd1" in text and "dox" in text:
            return "siARTD1 knockdown + doxycycline", "siartd1_knockdown_doxycycline"
        if "simock" in text and "dox" in text:
            return "siMock + doxycycline", "simock_doxycycline"
        if "simock" in text and ("untr" in text or "untreated" in text):
            return "siMock untreated", "simock_untreated"

    if project_id == "PRJNA623750":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        alias_lc = alias.lower()
        if alias_lc.startswith("a549_treated"):
            return "A549 treated", "a549_treated"
        if alias_lc.startswith("a549_untreated"):
            return "A549 untreated", "a549_untreated"
        if alias_lc.startswith("hct116_treated"):
            return "HCT116 treated", "hct116_treated"
        if alias_lc.startswith("hct116_untreated"):
            return "HCT116 untreated", "hct116_untreated"

    if project_id == "PRJNA625628":
        exp_title = clean_value(row.get("experiment_title", "")).lower()
        if " from lcl" in exp_title:
            return "LCL / blood-derived lymphoblastoid cells", "lcl"
        if " from fb" in exp_title:
            return "FB / skin fibroblasts", "fb"

    if project_id == "PRJNA634200":
        sample_title = clean_value(row.get("sample_title", ""))
        if sample_title.endswith("_RN"):
            return "replacement normal", "replacement_normal"
        if sample_title.endswith("_RT"):
            return "replacement tumor", "replacement_tumor"
        if sample_title.endswith("_DN"):
            return "desmoplastic normal", "desmoplastic_normal"
        if sample_title.endswith("_DT"):
            return "desmoplastic tumor", "desmoplastic_tumor"

    if project_id == "PRJNA664293":
        if "egfp-cag27" in text and "db213" in text:
            return "EGFP-CAG27 + DB213", "egfp_cag27_db213"
        if "egfp-cag78" in text and "db213" in text:
            return "EGFP-CAG78 + DB213", "egfp_cag78_db213"
        if "untransfected" in text and "db213" in text:
            return "untransfected SK-N-MC + DB213", "untransfected_db213"
        if "egfp-cag27" in text:
            return "EGFP-CAG27", "egfp_cag27"
        if "egfp-cag78" in text:
            return "EGFP-CAG78", "egfp_cag78"
        if "untransfected" in text:
            return "untransfected SK-N-MC", "untransfected"

    if project_id == "PRJNA673295":
        alias = clean_value(row.get("experiment_alias", "")).lower()
        if alias.startswith("tio2_ame"):
            return "AME + TiO2", "ame_tio2"
        if alias.startswith("tio2"):
            return "TiO2", "tio2"
        if alias.startswith("ame"):
            return "AME", "ame"
        if alias.startswith("com"):
            return "control / com", "control"

    if project_id == "PRJNA673745":
        alias = clean_value(row.get("experiment_alias", "")).lower()
        if alias.startswith("z_e"):
            return "ZnO2 + EGCG", "zno2_egcg"
        if alias.startswith("zno2"):
            return "ZnO2", "zno2"
        if alias.startswith("egcg"):
            return "EGCG", "egcg"
        if alias.startswith("con"):
            return "control / con", "control"

    if project_id == "PRJNA734133":
        if "kaiso ko" in text:
            return "Caki1 Kaiso KO", "caki1_kaiso_ko"
        if "caki1 wt" in text:
            return "Caki1 WT", "caki1_wt"

    if project_id == "PRJNA752868":
        if "naive treg" in text:
            return "naive Treg", "naive_treg"
        if "cd4 mature naive" in text:
            return "CD4 mature naive", "cd4_mature_naive"

    if project_id == "PRJNA767228":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        cell_line = alias.split("_", 1)[0].lower()
        breast_lines = {
            "mcf102a": ("MCF102A", "mcf102a"),
            "mcf7": ("MCF7", "mcf7"),
            "skbr3": ("SKBR3", "skbr3"),
            "mdamb231": ("MDAMB231", "mdamb231"),
            "mdamb361": ("MDAMB361", "mdamb361"),
        }
        if cell_line in breast_lines:
            return breast_lines[cell_line]

    if project_id == "PRJNA786266":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        alias_lc = alias.lower()
        if alias_lc.startswith("normal"):
            return "Normal endometrium", "normal_endometrium"
        if alias_lc.startswith("aeh"):
            return "AEH", "aeh"
        if alias_lc.startswith("eec"):
            return "EEC", "eec"

    if project_id == "PRJNA788948":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        alias_lc = alias.lower()
        primary_map = {
            "cd34_p6_gt39m": ("CD34+ cells P6 GT39m", "cd34_patient_gt39m"),
            "cd34_normal_donor": ("CD34+ cells normal donor", "cd34_normal_donor"),
            "cd14_p6_gt39m": ("CD14+ monocytes P6 GT39m", "cd14_monocytes"),
            "cd3_p6_gt39m": ("CD3+ T cells P6 GT39m", "cd3_t_cells"),
            "nk_p6_gt39m": ("NK cells P6 GT39m", "nk_cells"),
            "cd19_p6_gt39m": ("CD19+ B cells P6 GT39m", "cd19_b_cells"),
            "pmn_p6_gt39m": ("PMN cells P6 GT39m", "pmn_cells"),
        }
        if alias_lc in primary_map:
            return primary_map[alias_lc]
        is_macrophage = "macrophage" in alias_lc
        copy_text = " ".join(
            clean_value(row.get(col, ""))
            for col in ["sample_description", "sample_title", "experiment_title"]
        ).lower()
        copy_match = re.search(r"\b(no|\d+)\s+copy\b", copy_text)
        if copy_match:
            copy = "0" if copy_match.group(1) == "no" else copy_match.group(1)
            prefix = "iPSC-derived macrophage" if is_macrophage else "iPSC clone"
            condition_prefix = "ipsc_derived_macrophage" if is_macrophage else "ipsc_clone"
            return f"{prefix} {copy} IL2RG copies", f"{condition_prefix}_{copy}copy"

    if project_id == "PRJNA838478":
        sample_alias = clean_value(row.get("sample_alias", "")).lower()
        sample_title = clean_value(row.get("sample_title", "")).lower()
        if sample_title.startswith("rna-patient"):
            fraction = "cell" if "cell" in sample_alias else "supernatant" if "sup" in sample_alias else "unknown_fraction"
            return f"hydronephrosis patient {fraction}", f"hydronephrosis_patient_{fraction}"
        if sample_title.startswith("rna-control"):
            fraction = "cell" if "cell" in sample_alias else "supernatant" if "sup" in sample_alias else "unknown_fraction"
            return f"control {fraction}", f"control_{fraction}"
        return "unlabeled amniotic fluid", "unlabeled_amniotic_fluid"

    if project_id == "PRJNA876028":
        if "fto" in text:
            return "Kasumi-1 FTO overexpression", "kasumi1_fto_overexpression"
        if "ctrl" in text or "wildtype" in text:
            return "Kasumi-1 control / wildtype", "kasumi1_control"

    if project_id == "PRJNA903521":
        alias = clean_value(row.get("experiment_alias", "")).lower()
        if alias == "scrna_seq_control":
            return "scRNA control lung stem-cell pool", "scrna_control_lung_stem_pool"
        if alias.startswith("scrna_seq_ipf"):
            return "scRNA IPF lung stem-cell pool", "scrna_ipf_lung_stem_pool"
        if "ipf_nm_sc" in alias:
            return "IPF normal clone stem-cell state", "ipf_normal_clone_sc"
        if "ipf_nm_ali" in alias:
            return "IPF normal clone adult-cell ALI state", "ipf_normal_clone_ali"
        if "ipf_fm_sc" in alias:
            return "IPF fibrosis clone stem-cell state", "ipf_fibrosis_clone_sc"
        if "ipf_fm_ali" in alias:
            return "IPF fibrosis clone adult-cell ALI state", "ipf_fibrosis_clone_ali"

    if project_id == "PRJNA976462":
        alias = clean_value(row.get("experiment_alias", ""))
        sample_alias = clean_value(row.get("sample_alias", ""))
        patient = ""
        if alias.lower().startswith("pt"):
            patient = alias
        else:
            match = re.search(r"pt(\d+)", sample_alias, flags=re.IGNORECASE)
            if match:
                patient = f"Pt{match.group(1)}"
        if patient:
            return patient, normalize_condition(patient)

    if project_id == "PRJNA987832":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        condition = re.sub(r"_\d+$", "", alias)
        ferroptosis_map = {
            "293FT": ("293FT", "293ft"),
            "293FT-Erastin-R": ("293FT-Erastin-R", "293ft_erastin_r"),
            "293FT-RSL3-R": ("293FT-RSL3-R", "293ft_rsl3_r"),
            "U2OS": ("U2OS", "u2os"),
            "U2OS-Erastin-R": ("U2OS-Erastin-R", "u2os_erastin_r"),
            "U2OS-RSL3-R": ("U2OS-RSL3-R", "u2os_rsl3_r"),
        }
        if condition in ferroptosis_map:
            return ferroptosis_map[condition]

    if project_id == "PRJNA996357":
        alias = clean_value(row.get("experiment_alias", "")) or clean_value(row.get("sample_alias", ""))
        if alias.startswith("CSF_"):
            return alias, normalize_condition(alias)

    return None


def clean_value(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    text = re.sub(r"\s+", " ", text)
    return text


def normalize_condition(value: str) -> str:
    text = clean_value(value)
    text = re.sub(r"[^A-Za-z0-9]+", "_", text.lower()).strip("_")
    return text or "unknown"


def meaningful(value: object) -> bool:
    text = clean_value(value).lower()
    return text not in BAD_VALUES


def fetch_ena_project(project_id: str, retries: int = 3) -> pd.DataFrame:
    params = urllib.parse.urlencode(
        {
            "accession": project_id,
            "result": "read_run",
            "fields": ",".join(FIELDS),
            "format": "tsv",
            "limit": "0",
        }
    )
    url = f"{ENA_URL}?{params}"
    last_error = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "omiks-deseq2-metadata/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return pd.read_csv(resp, sep="\t", dtype=str).fillna("")
        except Exception as exc:
            last_error = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"failed to fetch {project_id}: {last_error}")


def simple_value_counts(df: pd.DataFrame, field: str) -> Counter[str]:
    values = []
    for raw in df[field].tolist():
        if meaningful(raw):
            values.append(clean_value(raw))
    return Counter(values)


def specificity_score(values: Counter[str], n: int) -> float:
    if not values:
        return -1.0
    classes = len(values)
    largest = max(values.values())
    smallest = min(values.values())
    singletons = sum(1 for v in values.values() if v == 1)
    # Prefer a small number of replicated classes; penalize near-unique labels.
    return (
        20.0
        - classes * 1.8
        - singletons * 2.5
        + min(smallest, 5) * 1.2
        - (largest / max(n, 1)) * 1.5
    )


def choose_class_field(df: pd.DataFrame) -> tuple[str, Counter[str], str]:
    best_field = ""
    best_counts: Counter[str] = Counter()
    best_score = -999.0
    best_status = "no_candidate"
    n = len(df)

    for field in CLASS_FIELD_PRIORITY:
        if field not in df.columns:
            continue
        counts = simple_value_counts(df, field)
        if len(counts) < 2:
            continue
        if sum(counts.values()) < max(4, int(n * 0.5)):
            continue
        if len(counts) > min(12, max(4, n // 2)):
            continue
        score = specificity_score(counts, n)
        priority_bonus = (len(CLASS_FIELD_PRIORITY) - CLASS_FIELD_PRIORITY.index(field)) * 0.05
        score += priority_bonus
        if score > best_score:
            best_field = field
            best_counts = counts
            best_score = score

    if best_field:
        if len(best_counts) == 2:
            best_status = "two_class_candidate"
        else:
            best_status = "multi_class_candidate"
    return best_field, best_counts, best_status


def label_role(condition: str, all_conditions: list[str]) -> str:
    text = condition.lower().replace("_", " ")
    control = any(re.search(pattern, text) for pattern in CONTROL_PATTERNS)
    treatment = any(re.search(pattern, text) for pattern in TREATMENT_PATTERNS)
    if control and not treatment:
        return "control"
    if treatment and not control:
        return "treatment"
    if len(all_conditions) == 2:
        other = [c for c in all_conditions if c != condition]
        if other and any(re.search(pattern, other[0].lower().replace("_", " ")) for pattern in CONTROL_PATTERNS):
            return "treatment"
    return "class"


def build_project_sheet(project_id: str, current_runs: pd.DataFrame, skip_fetch: bool = False) -> dict[str, object]:
    if skip_fetch:
        ena = pd.DataFrame()
    else:
        ena = fetch_ena_project(project_id)
        time.sleep(RATE_LIMIT_SECONDS)

    if ena.empty:
        ena = current_runs.rename(columns={"SRR_ID": "run_accession"}).copy()
    else:
        ena = ena.rename(columns={"run_accession": "SRR_ID"})
        ena = current_runs[["SRR_ID", "project_id"]].merge(ena, on="SRR_ID", how="left")

    ena["project_id"] = project_id
    for field in FIELDS:
        col = "SRR_ID" if field == "run_accession" else field
        if col not in ena.columns:
            ena[col] = ""
    ena = ena.fillna("")

    manual = ena.apply(lambda row: manual_project_condition(project_id, row), axis=1)
    if manual.notna().all():
        ena["condition_raw"] = manual.map(lambda item: item[0])
        ena["condition"] = manual.map(lambda item: item[1])
        class_field = "manual_override"
        counts = Counter(ena["condition"])
        if len(counts) == 1:
            status = "manual_not_for_deseq2" if project_id in MANUAL_NOT_FOR_DESEQ2 else "manual_one_class"
        elif len(counts) == 2:
            status = "manual_two_class"
        else:
            status = "manual_multi_class"
    else:
        class_field, counts, status = choose_class_field(ena)
    if class_field and "condition" not in ena.columns:
        raw = ena[class_field].map(clean_value)
        ena["condition_raw"] = raw
        ena["condition"] = raw.map(normalize_condition)
    elif "condition" not in ena.columns:
        ena["condition_raw"] = "unknown"
        ena["condition"] = "unknown"

    conditions = sorted(c for c in ena["condition"].unique() if c != "unknown")
    ena["condition_role"] = ena["condition"].map(lambda c: label_role(c, conditions) if c != "unknown" else "unknown")
    ena["class_field"] = class_field
    ena["metadata_status"] = status

    full_file = FULL_DIR / f"{project_id}_metadata.tsv"
    ena.to_csv(full_file, sep="\t", index=False)

    sheet_cols = ["SRR_ID", "condition", "condition_raw", "condition_role", "class_field", "metadata_status"]
    sheet = ena[sheet_cols].rename(columns={"SRR_ID": "Run"}).copy()
    sheet_file = SHEET_DIR / f"{project_id}.csv"
    sheet.to_csv(sheet_file, index=False)

    review_cols = [
        "SRR_ID",
        "condition",
        "condition_raw",
        "condition_role",
        "class_field",
        "sample_title",
        "sample_description",
        "experiment_title",
        "disease",
        "host_status",
        "host_phenotype",
        "host_genotype",
        "cell_type",
        "cell_line",
        "tissue_type",
        "isolation_source",
        "experimental_factor",
        "control_experiment",
    ]
    review_file = REVIEW_DIR / f"{project_id}_review.tsv"
    ena[[c for c in review_cols if c in ena.columns]].to_csv(review_file, sep="\t", index=False)

    condition_counts = Counter(ena["condition"])
    role_counts = Counter(ena["condition_role"])
    usable_two_class = (
        status in {"two_class_candidate", "manual_two_class"}
        and len([c for c in condition_counts if c != "unknown"]) == 2
        and min(v for c, v in condition_counts.items() if c != "unknown") >= 3
    )

    return {
        "project_id": project_id,
        "n_runs": len(ena),
        "class_field": class_field or "NA",
        "metadata_status": status,
        "n_classes": len([c for c in condition_counts if c != "unknown"]),
        "condition_counts": "; ".join(f"{k}:{v}" for k, v in sorted(condition_counts.items())),
        "role_counts": "; ".join(f"{k}:{v}" for k, v in sorted(role_counts.items())),
        "usable_two_class_min3": "yes" if usable_two_class else "no",
        "sample_sheet": str(sheet_file),
        "review_file": str(review_file),
        "full_metadata": str(full_file),
    }


def load_current_projects(project_filter: str = "") -> pd.DataFrame:
    current = pd.read_csv(ANALYSIS / "current_srr_platform_check.tsv", sep="\t", dtype=str).fillna("")
    eval_srrs = pd.read_csv(ANALYSIS / "per_srr_eval.tsv", sep="\t", usecols=["SRR_ID", "project_id"], dtype=str)
    current = eval_srrs.merge(current, on=["SRR_ID", "project_id"], how="left")
    current = current.loc[(current["library_strategy"] == "RNA-Seq") & (current["is_illumina"].str.lower().isin(["yes", "true", "1"]))]
    if project_filter:
        current = current.loc[current["project_id"] == project_filter]
    return current[["SRR_ID", "project_id"]].drop_duplicates()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default="", help="Only process one BioProject")
    parser.add_argument("--skip-fetch", action="store_true", help="Do not call ENA; mainly for testing")
    args = parser.parse_args()

    for directory in [OUT_DIR, FULL_DIR, SHEET_DIR, REVIEW_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    current = load_current_projects(args.project)
    rows = []
    for project_id, group in current.groupby("project_id", sort=True):
        print(f"Processing {project_id} ({len(group)} RNA-seq SRRs)...", flush=True)
        try:
            rows.append(build_project_sheet(project_id, group, skip_fetch=args.skip_fetch))
        except Exception as exc:
            rows.append(
                {
                    "project_id": project_id,
                    "n_runs": len(group),
                    "class_field": "NA",
                    "metadata_status": f"fetch_or_parse_failed: {exc}",
                    "n_classes": 0,
                    "condition_counts": "",
                    "role_counts": "",
                    "usable_two_class_min3": "no",
                    "sample_sheet": "",
                    "review_file": "",
                    "full_metadata": "",
                }
            )

    summary = pd.DataFrame(rows)
    summary_file = OUT_DIR / "deseq2_metadata_summary.tsv"
    summary.to_csv(summary_file, sep="\t", index=False)
    print(f"Written: {summary_file}")
    print(summary[["project_id", "n_runs", "class_field", "metadata_status", "n_classes", "usable_two_class_min3"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
