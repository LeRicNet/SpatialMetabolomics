#' @title Metabolic pathway gene sets
#' @description Functions to retrieve curated metabolic pathway gene sets

#' Get oxidative phosphorylation genes
#'
#' @param species Species name ("mouse" or "human")
#' @return Character vector of gene symbols
#' @export
#'
#' @examples
#' # Get mouse OXPHOS genes
#' oxphos_genes <- getOXPHOSGenes("mouse")
#' length(oxphos_genes)  # Should be ~140 genes
getOXPHOSGenes <- function(species = "mouse") {
  species <- tolower(species)

  if (species == "mouse") {
    # Complex I genes (NADH dehydrogenase)
    complex_i <- c("Ndufa1", "Ndufa2", "Ndufa3", "Ndufa4", "Ndufa5", "Ndufa6",
                   "Ndufa7", "Ndufa8", "Ndufa9", "Ndufa10", "Ndufa11", "Ndufa12",
                   "Ndufa13", "Ndufab1", "Ndufb1", "Ndufb2", "Ndufb3", "Ndufb4",
                   "Ndufb5", "Ndufb6", "Ndufb7", "Ndufb8", "Ndufb9", "Ndufb10",
                   "Ndufb11", "Ndufc1", "Ndufc2", "Ndufs1", "Ndufs2", "Ndufs3",
                   "Ndufs4", "Ndufs5", "Ndufs6", "Ndufs7", "Ndufs8", "Ndufv1",
                   "Ndufv2", "Ndufv3")

    # Complex II genes (Succinate dehydrogenase)
    complex_ii <- c("Sdha", "Sdhb", "Sdhc", "Sdhd", "Sdhaf1", "Sdhaf2", "Sdhaf3", "Sdhaf4")

    # Complex III genes (Cytochrome bc1 complex)
    complex_iii <- c("Uqcrc1", "Uqcrc2", "Uqcrfs1", "Uqcrb", "Uqcrq", "Uqcrh",
                     "Uqcr10", "Uqcr11", "Cyc1", "Cycs", "Uqcc1", "Uqcc2", "Uqcc3")

    # Complex IV genes (Cytochrome c oxidase)
    complex_iv <- c("Cox4i1", "Cox4i2", "Cox5a", "Cox5b", "Cox6a1", "Cox6a2",
                    "Cox6b1", "Cox6b2", "Cox6c", "Cox7a1", "Cox7a2", "Cox7a2l",
                    "Cox7b", "Cox7b2", "Cox7c", "Cox8a", "Cox8b", "Cox8c",
                    "Cox10", "Cox11", "Cox14", "Cox15", "Cox16", "Cox17",
                    "Cox18", "Cox19", "Cox20")

    # Complex V genes (ATP synthase)
    complex_v <- c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atp5f1",
                   "Atp5g1", "Atp5g2", "Atp5g3", "Atp5h", "Atp5i", "Atp5j",
                   "Atp5j2", "Atp5k", "Atp5l", "Atp5o", "Atp5s", "Atp5pb",
                   "Atp5pf", "Atp5pg", "Atp5pd", "Atp5me", "Atp5mf", "Atp5mg",
                   "Atp5if1", "Atpaf1", "Atpaf2")

    return(c(complex_i, complex_ii, complex_iii, complex_iv, complex_v))

  } else if (species == "human") {
    # Human orthologs - uppercase nomenclature
    complex_i <- c("NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6",
                   "NDUFA7", "NDUFA8", "NDUFA9", "NDUFA10", "NDUFA11", "NDUFA12",
                   "NDUFA13", "NDUFAB1", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4",
                   "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10",
                   "NDUFB11", "NDUFC1", "NDUFC2", "NDUFS1", "NDUFS2", "NDUFS3",
                   "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1",
                   "NDUFV2", "NDUFV3")

    complex_ii <- c("SDHA", "SDHB", "SDHC", "SDHD", "SDHAF1", "SDHAF2", "SDHAF3", "SDHAF4")

    complex_iii <- c("UQCRC1", "UQCRC2", "UQCRFS1", "UQCRB", "UQCRQ", "UQCRH",
                     "UQCR10", "UQCR11", "CYC1", "CYCS", "UQCC1", "UQCC2", "UQCC3")

    complex_iv <- c("COX4I1", "COX4I2", "COX5A", "COX5B", "COX6A1", "COX6A2",
                    "COX6B1", "COX6B2", "COX6C", "COX7A1", "COX7A2", "COX7A2L",
                    "COX7B", "COX7B2", "COX7C", "COX8A", "COX8C",
                    "COX10", "COX11", "COX14", "COX15", "COX16", "COX17",
                    "COX18", "COX19", "COX20")

    complex_v <- c("ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
                   "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MF",
                   "ATP5MG", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO",
                   "ATP5IF1", "ATPAF1", "ATPAF2")

    return(c(complex_i, complex_ii, complex_iii, complex_iv, complex_v))
  } else {
    stop("Species must be 'mouse' or 'human'")
  }
}

#' Get glycolysis pathway genes
#'
#' @param species Species name ("mouse" or "human")
#' @return Character vector of gene symbols
#' @export
#'
#' @examples
#' # Get mouse glycolysis genes
#' glycolysis_genes <- getGlycolysisGenes("mouse")
getGlycolysisGenes <- function(species = "mouse") {
  species <- tolower(species)

  if (species == "mouse") {
    return(c(
      # Glucose transporters
      "Slc2a1", "Slc2a3", "Slc2a4",
      # Glycolytic enzymes
      "Hk1", "Hk2", "Hk3", "Gck",             # Hexokinase
      "Gpi1",                                   # Glucose-6-phosphate isomerase
      "Pfkl", "Pfkm", "Pfkp", "Pfkfb1",       # Phosphofructokinase
      "Pfkfb2", "Pfkfb3", "Pfkfb4",
      "Aldoa", "Aldob", "Aldoc",               # Aldolase
      "Tpi1",                                   # Triose phosphate isomerase
      "Gapdh", "Gapdhs",                       # GAPDH
      "Pgk1", "Pgk2",                          # Phosphoglycerate kinase
      "Pgam1", "Pgam2", "Pgam5",               # Phosphoglycerate mutase
      "Eno1", "Eno2", "Eno3",                  # Enolase
      "Pkm", "Pklr",                           # Pyruvate kinase
      "Ldha", "Ldhb", "Ldhc"                   # Lactate dehydrogenase
    ))
  } else if (species == "human") {
    return(c(
      # Glucose transporters
      "SLC2A1", "SLC2A3", "SLC2A4",
      # Glycolytic enzymes
      "HK1", "HK2", "HK3", "GCK",
      "GPI",
      "PFKL", "PFKM", "PFKP", "PFKFB1",
      "PFKFB2", "PFKFB3", "PFKFB4",
      "ALDOA", "ALDOB", "ALDOC",
      "TPI1",
      "GAPDH", "GAPDHS",
      "PGK1", "PGK2",
      "PGAM1", "PGAM2", "PGAM5",
      "ENO1", "ENO2", "ENO3",
      "PKM", "PKLR",
      "LDHA", "LDHB", "LDHC"
    ))
  } else {
    stop("Species must be 'mouse' or 'human'")
  }
}

#' Get mitochondrial biogenesis genes
#'
#' @param species Species name ("mouse" or "human")
#' @return Character vector of gene symbols
#' @export
getMitoBiogenesisGenes <- function(species = "mouse") {
  species <- tolower(species)

  if (species == "mouse") {
    return(c(
      # Master regulators
      "Ppargc1a", "Ppargc1b", "Pprc1",        # PGC-1 family
      "Nrf1", "Nfe2l2",                        # Nuclear respiratory factors
      "Tfam", "Tfb1m", "Tfb2m",                # Mitochondrial transcription factors
      "Polrmt",                                 # Mitochondrial RNA polymerase
      "Mterf1", "Mterf2", "Mterf3", "Mterf4", # Mitochondrial transcription terminators
      # Co-activators and regulators
      "Sirt1", "Sirt3", "Sirt5",               # Sirtuins
      "Esrra", "Esrrb", "Esrrg",               # Estrogen-related receptors
      "Gabpa", "Gabpb1", "Gabpb2",             # GA-binding proteins
      "Yy1",                                    # Yin Yang 1
      "Camk4",                                  # Ca2+/calmodulin-dependent kinase
      "Prkaa1", "Prkaa2",                      # AMPK catalytic subunits
      "Creb1", "Creb3", "Creb5",               # CREB family
      "Myc", "Mycn"                            # MYC family
    ))
  } else if (species == "human") {
    return(c(
      "PPARGC1A", "PPARGC1B", "PPRC1",
      "NRF1", "NFE2L2",
      "TFAM", "TFB1M", "TFB2M",
      "POLRMT",
      "MTERF1", "MTERF2", "MTERF3", "MTERF4",
      "SIRT1", "SIRT3", "SIRT5",
      "ESRRA", "ESRRB", "ESRRG",
      "GABPA", "GABPB1", "GABPB2",
      "YY1",
      "CAMK4",
      "PRKAA1", "PRKAA2",
      "CREB1", "CREB3", "CREB5",
      "MYC", "MYCN"
    ))
  } else {
    stop("Species must be 'mouse' or 'human'")
  }
}

#' Get mitochondrial dynamics genes
#'
#' @param species Species name ("mouse" or "human")
#' @return Character vector of gene symbols
#' @export
getMitoDynamicsGenes <- function(species = "mouse") {
  species <- tolower(species)

  if (species == "mouse") {
    # Fission genes
    fission <- c("Dnm1l", "Fis1", "Mff", "Mief1", "Mief2", "Inf2", "Gdap1")

    # Fusion genes
    fusion <- c("Mfn1", "Mfn2", "Opa1", "Opa3", "Slc25a46", "Yme1l1")

    # Mitophagy genes
    mitophagy <- c("Pink1", "Park2", "Park7", "Bnip3", "Bnip3l", "Fundc1",
                   "Atg5", "Atg7", "Atg12", "Map1lc3a", "Map1lc3b", "Sqstm1",
                   "Optn", "Nbr1", "Ambra1", "Ulk1", "Becn1")

    # Mitochondrial motility
    motility <- c("Miro1", "Miro2", "Trak1", "Trak2", "Kif5b", "Dync1h1")

    return(c(fission, fusion, mitophagy, motility))
  } else if (species == "human") {
    fission <- c("DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2", "INF2", "GDAP1")
    fusion <- c("MFN1", "MFN2", "OPA1", "OPA3", "SLC25A46", "YME1L1")
    mitophagy <- c("PINK1", "PRKN", "PARK7", "BNIP3", "BNIP3L", "FUNDC1",
                   "ATG5", "ATG7", "ATG12", "MAP1LC3A", "MAP1LC3B", "SQSTM1",
                   "OPTN", "NBR1", "AMBRA1", "ULK1", "BECN1")
    motility <- c("RHOT1", "RHOT2", "TRAK1", "TRAK2", "KIF5B", "DYNC1H1")

    return(c(fission, fusion, mitophagy, motility))
  } else {
    stop("Species must be 'mouse' or 'human'")
  }
}

#' Get ATP synthesis and consumption genes
#'
#' @param species Species name ("mouse" or "human")
#' @return Character vector of gene symbols
#' @export
getATPGenes <- function(species = "mouse") {
  species <- tolower(species)

  if (species == "mouse") {
    # ATP synthesis (see Complex V in OXPHOS)
    synthesis <- c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atp5f1",
                   "Atp5g1", "Atp5g2", "Atp5g3", "Atp5h", "Atp5i", "Atp5j",
                   "Atp5j2", "Atp5k", "Atp5l", "Atp5o", "Atp5s",
                   # Adenine nucleotide translocases
                   "Slc25a4", "Slc25a5", "Slc25a6", "Slc25a31",
                   # Creatine kinases
                   "Ckb", "Ckm", "Ckmt1", "Ckmt2")

    # Major ATP consumers
    consumers <- c(
      # Na+/K+-ATPase
      "Atp1a1", "Atp1a2", "Atp1a3", "Atp1a4", "Atp1b1", "Atp1b2", "Atp1b3",
      # Ca2+-ATPase
      "Atp2a1", "Atp2a2", "Atp2a3", "Atp2b1", "Atp2b2", "Atp2b3", "Atp2b4",
      # V-ATPase subunits
      "Atp6v1a", "Atp6v1b1", "Atp6v1b2", "Atp6v1c1", "Atp6v1d", "Atp6v1e1",
      # Protein synthesis
      "Eef1a1", "Eef1a2", "Eef2",
      # Chaperones
      "Hspa1a", "Hspa1b", "Hspa8", "Hspa9", "Hsp90aa1", "Hsp90ab1"
    )

    return(c(synthesis, consumers))
  } else if (species == "human") {
    synthesis <- c("ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
                   "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MF",
                   "ATP5MG", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO",
                   "SLC25A4", "SLC25A5", "SLC25A6", "SLC25A31",
                   "CKB", "CKM", "CKMT1A", "CKMT1B", "CKMT2")

    consumers <- c("ATP1A1", "ATP1A2", "ATP1A3", "ATP1A4", "ATP1B1", "ATP1B2", "ATP1B3",
                   "ATP2A1", "ATP2A2", "ATP2A3", "ATP2B1", "ATP2B2", "ATP2B3", "ATP2B4",
                   "ATP6V1A", "ATP6V1B1", "ATP6V1B2", "ATP6V1C1", "ATP6V1D", "ATP6V1E1",
                   "EEF1A1", "EEF1A2", "EEF2",
                   "HSPA1A", "HSPA1B", "HSPA8", "HSPA9", "HSP90AA1", "HSP90AB1")

    return(c(synthesis, consumers))
  } else {
    stop("Species must be 'mouse' or 'human'")
  }
}

#' Get small vessel disease signature genes
#'
#' @param species Species name ("mouse" or "human")
#' @return Character vector of gene symbols
#' @export
getSVDGenes <- function(species = "mouse") {
  species <- tolower(species)

  if (species == "mouse") {
    # Hypoxia response
    hypoxia <- c("Hif1a", "Epas1", "Hif3a", "Vhl", "Egln1", "Egln2", "Egln3",
                 "Vegfa", "Vegfb", "Vegfc", "Pgk1", "Ldha", "Slc2a1", "Slc2a3",
                 "Eno1", "Pdk1", "Bnip3", "Bnip3l", "P4ha1", "P4ha2", "Plod1", "Plod2")

    # Vascular integrity and blood-brain barrier
    vascular <- c("Col4a1", "Col4a2", "Col3a1", "Col1a1", "Col1a2",
                  "Notch3", "Jag1", "Jag2", "Hey1", "Hey2", "Heyl",
                  "Pecam1", "Cdh5", "Cldn5", "Ocln", "Tjp1", "Tjp2",
                  "Acta2", "Tagln", "Myh11", "Cnn1", "Smtn",
                  "Pdgfrb", "Pdgfb", "Angpt1", "Angpt2", "Tek")

    # Inflammation and immune response
    inflammation <- c("Il1b", "Il6", "Il18", "Tnf", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
                      "Cxcl1", "Cxcl2", "Cxcl10", "Cxcl12",
                      "Icam1", "Vcam1", "Sele", "Selp",
                      "Nfkb1", "Rela", "Tlr2", "Tlr4", "Myd88")

    # Oxidative stress and antioxidant response
    oxidative <- c("Sod1", "Sod2", "Sod3", "Cat", "Gpx1", "Gpx3", "Gpx4",
                   "Prdx1", "Prdx2", "Prdx3", "Prdx5", "Prdx6",
                   "Gsr", "Gss", "Txn1", "Txn2", "Txnrd1", "Txnrd2",
                   "Nqo1", "Hmox1", "Gclc", "Gclm")

    # White matter and oligodendrocyte markers
    white_matter <- c("Mbp", "Plp1", "Mag", "Mog", "Cnp", "Olig1", "Olig2",
                      "Sox10", "Nkx2-2", "Pdgfra", "Cspg4")

    return(c(hypoxia, vascular, inflammation, oxidative, white_matter))

  } else if (species == "human") {
    hypoxia <- c("HIF1A", "EPAS1", "HIF3A", "VHL", "EGLN1", "EGLN2", "EGLN3",
                 "VEGFA", "VEGFB", "VEGFC", "PGK1", "LDHA", "SLC2A1", "SLC2A3",
                 "ENO1", "PDK1", "BNIP3", "BNIP3L", "P4HA1", "P4HA2", "PLOD1", "PLOD2")

    vascular <- c("COL4A1", "COL4A2", "COL3A1", "COL1A1", "COL1A2",
                  "NOTCH3", "JAG1", "JAG2", "HEY1", "HEY2", "HEYL",
                  "PECAM1", "CDH5", "CLDN5", "OCLN", "TJP1", "TJP2",
                  "ACTA2", "TAGLN", "MYH11", "CNN1", "SMTN",
                  "PDGFRB", "PDGFB", "ANGPT1", "ANGPT2", "TEK")

    inflammation <- c("IL1B", "IL6", "IL18", "TNF", "CCL2", "CCL3", "CCL4", "CCL5",
                      "CXCL1", "CXCL2", "CXCL10", "CXCL12",
                      "ICAM1", "VCAM1", "SELE", "SELP",
                      "NFKB1", "RELA", "TLR2", "TLR4", "MYD88")

    oxidative <- c("SOD1", "SOD2", "SOD3", "CAT", "GPX1", "GPX3", "GPX4",
                   "PRDX1", "PRDX2", "PRDX3", "PRDX5", "PRDX6",
                   "GSR", "GSS", "TXN", "TXN2", "TXNRD1", "TXNRD2",
                   "NQO1", "HMOX1", "GCLC", "GCLM")

    white_matter <- c("MBP", "PLP1", "MAG", "MOG", "CNP", "OLIG1", "OLIG2",
                      "SOX10", "NKX2-2", "PDGFRA", "CSPG4")

    return(c(hypoxia, vascular, inflammation, oxidative, white_matter))
  } else {
    stop("Species must be 'mouse' or 'human'")
  }
}

#' List all available metabolic pathways
#'
#' @param species Species name ("mouse" or "human")
#' @return Data frame with pathway names and gene counts
#' @export
#'
#' @examples
#' # List all pathways
#' pathways <- listMetabolicPathways("mouse")
#' print(pathways)
listMetabolicPathways <- function(species = "mouse") {
  pathways <- list(
    OXPHOS = getOXPHOSGenes(species),
    Glycolysis = getGlycolysisGenes(species),
    Mito_Biogenesis = getMitoBiogenesisGenes(species),
    Mito_Dynamics = getMitoDynamicsGenes(species),
    ATP_Metabolism = getATPGenes(species),
    SVD_Signature = getSVDGenes(species)
  )

  df <- data.frame(
    pathway = names(pathways),
    n_genes = sapply(pathways, length),
    stringsAsFactors = FALSE
  )

  return(df)
}

#' Get all metabolic pathways as a list
#'
#' @param species Species name ("mouse" or "human")
#' @param pathways Character vector of pathway names to include
#' @return Named list of gene sets
#' @export
#'
#' @examples
#' # Get all pathways
#' all_pathways <- getAllMetabolicPathways("mouse")
#'
#' # Get specific pathways
#' energy_pathways <- getAllMetabolicPathways("mouse",
#'                                           c("OXPHOS", "Glycolysis", "ATP_Metabolism"))
getAllMetabolicPathways <- function(species = "mouse",
                                    pathways = c("OXPHOS", "Glycolysis", "Mito_Biogenesis",
                                                 "Mito_Dynamics", "ATP_Metabolism", "SVD_Signature")) {

  pathway_list <- list()

  for (pathway in pathways) {
    if (pathway == "OXPHOS") {
      pathway_list[[pathway]] <- getOXPHOSGenes(species)
    } else if (pathway == "Glycolysis") {
      pathway_list[[pathway]] <- getGlycolysisGenes(species)
    } else if (pathway == "Mito_Biogenesis") {
      pathway_list[[pathway]] <- getMitoBiogenesisGenes(species)
    } else if (pathway == "Mito_Dynamics") {
      pathway_list[[pathway]] <- getMitoDynamicsGenes(species)
    } else if (pathway == "ATP_Metabolism") {
      pathway_list[[pathway]] <- getATPGenes(species)
    } else if (pathway == "SVD_Signature") {
      pathway_list[[pathway]] <- getSVDGenes(species)
    } else {
      warning("Unknown pathway: ", pathway)
    }
  }

  return(pathway_list)
}
