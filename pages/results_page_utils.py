"""Python script containing static variables used in CRISPRme's result page. 
"""

# number of entries in report table (for each table page)
PAGE_SIZE = 10  
# number of barplots in each row of Populations Distributions
BARPLOT_LEN = 4  
# column names for reference report
COL_REF = [
    "Bulge Type", 
    "crRNA", 
    "Off target_motif", 
    "Reference sequence", 
    "Chromosome",
    "Position", 
    "Direction", 
    "Mismatches",
    "Bulge Size", 
    "PAM gen", 
    "Samples", 
    "Variant",
    "CFD", 
    "CFD ref", 
    "Highest CFD Risk Score",
    "AF", 
    "Annotation Type"
]
# reference column types
COL_REF_TYPE = [
    "text", 
    "text", 
    "text", 
    "text", 
    "text", 
    "numeric",
    "numeric", 
    "text", 
    "numeric", 
    "numeric", 
    "text", 
    "text", 
    "text",
    "numeric", 
    "numeric", 
    "numeric", 
    "numeric", 
    "text"
]
# reference columns renaming
COL_REF_RENAME = {
    0: "Bulge Type", 
    1: "crRNA", 
    2: "Off target motif", 
    3: "Reference sequence", 
    4: "Chromosome", 
    5: "Position", 
    6: "Cluster Position", 
    7: "Direction",
    8: "Mismatches", 
    9: "Bulge Size", 
    10: "Total", 
    11: "PAM gen", 
    12: "Variant Unique", 
    13: "Samples", 
    14: "Annotation Type", 
    15: "Real Guide",
    16: "rsID", 
    17: "AF", 
    18: "Variant", 
    19: "#Seq in cluster", 
    20: "CFD", 
    21: "CFD ref", 
    22: "Highest CFD Risk Score"
}
# reference and non reference columns
COL_BOTH = [
    "Highest_CFD_Strand", 
    "Chromosome", 
    "Highest_CFD_start_coordinate",
    "Highest_CFD_aligned_spacer+PAM",
    "Highest_CFD_aligned_protospacer+PAM_REF", 
    "Highest_CFD_aligned_protospacer+PAM_ALT",
    "Highest_CFD_mismatches", 
    "Highest_CFD_bulges", 
    "Highest_CFD_mismatches+bulges",
    "Highest_CFD_bulge_type", 
    "Highest_CFD_PAM_gen", 
    "Highest_CFD_score", 
    "Highest_CFD_score_REF",
    "Highest_CFD_risk_score", 
    "Not_found_in_REF", 
    "Highest_CFD_variant_info_genome",
    "Highest_CFD_variant_MAF", 
    "Highest_CFD_variant_rsID",
    "Highest_CFD_variant_samples", 
    "Other_motifs", 
    "Annotation_ENCODE"
]
# reference and non reference column types
COL_BOTH_TYPE = [
    "text", 
    "text", 
    "numeric", 
    "text",
    "text", 
    "text",
    "numeric", 
    "numeric", 
    "numeric",
    "text", 
    "text", 
    "numeric", 
    "numeric",
    "numeric", 
    "text", 
    "text", 
    "numeric", 
    "text",
    "text", 
    "numeric", 
    "text"
]
# reference and non reference column renaming 
COL_BOTH_RENAME = {
    0: "Highest_CFD_Strand", 
    1: "Chromosome", 
    2: "Highest_CFD_start_coordinate",
    3: "Highest_CFD_aligned_spacer+PAM", 
    4: "Highest_CFD_aligned_protospacer+PAM_REF",
    5: "Highest_CFD_aligned_protospacer+PAM_ALT", 
    6: "Highest_CFD_mismatches",
    7: "Highest_CFD_bulges", 
    8: "Highest_CFD_mismatches+bulges",
    9: "Highest_CFD_bulge_type", 
    10: "Highest_CFD_PAM_gen", 
    11: "Highest_CFD_score",
    12: "Highest_CFD_score_REF", 
    13: "Highest_CFD_risk_score",
    14: "Not_found_in_REF", 
    15: "Highest_CFD_variant_info_genome",
    16: "Highest_CFD_variant_MAF", 
    17: "Highest_CFD_variant_rsID",
    18: "Highest_CFD_variant_samples", 
    19: "Other_motifs", 
    37: "Annotation_ENCODE"
}
# genome database fields
GENOME_DATABASE = [
    "Reference", "Enriched", "Samples", "Dictionary", "Annotation"
]
# guide column name
GUIDE_COLUMN = "Spacer+PAM"
# chromosome column name
CHR_COLUMN = "Chromosome"
# position column name
POS_COLUMN = "Start_coordinate_(highest_CFD)"
# mismatches column name
MM_COLUMN = "Mismatches_(highest_CFD)"
# bulges column name
BLG_COLUMN = "Bulges_(highest_CFD)"
# total column name
TOTAL_COLUMN = "Mismatches+bulges_(highest_CFD)"
# bulge type column name
BLG_T_COLUMN = "Bulge_type_(highest_CFD)"
# CFD score column name
CFD_COLUMN = "CFD_score_(highest_CFD)"
# CFD risk score column name
RISK_COLUMN = "CFD_risk_score_(highest_CFD)"
# variant samples column name
SAMPLES_COLUMN = "Variant_samples_(highest_CFD)"

