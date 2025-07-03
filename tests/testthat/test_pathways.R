# Test metabolic pathway functions

test_that("Pathway gene retrieval works correctly", {
  # Test mouse pathways
  oxphos_mouse <- getOXPHOSGenes("mouse")
  expect_true(is.character(oxphos_mouse))
  expect_true(length(oxphos_mouse) > 100)
  expect_true(all(grepl("^[A-Z][a-z]", oxphos_mouse)))  # Mouse gene naming

  # Check for known genes
  expect_true("Ndufa1" %in% oxphos_mouse)
  expect_true("Atp5a1" %in% oxphos_mouse)
  expect_true("Cox4i1" %in% oxphos_mouse)

  # Test human pathways
  oxphos_human <- getOXPHOSGenes("human")
  expect_true(is.character(oxphos_human))
  expect_true(all(grepl("^[A-Z0-9]+$", oxphos_human)))  # Human gene naming
  expect_true("NDUFA1" %in% oxphos_human)

  # Test glycolysis
  glyc_mouse <- getGlycolysisGenes("mouse")
  expect_true(length(glyc_mouse) >= 20)
  expect_true("Hk1" %in% glyc_mouse)
  expect_true("Gapdh" %in% glyc_mouse)
  expect_true("Pkm" %in% glyc_mouse)

  # Test error handling
  expect_error(getOXPHOSGenes("invalid"), "Species must be")
})

test_that("Mitochondrial pathways are comprehensive", {
  # Biogenesis
  biogen <- getMitoBiogenesisGenes("mouse")
  expect_true("Ppargc1a" %in% biogen)  # PGC-1α
  expect_true("Tfam" %in% biogen)
  expect_true(length(biogen) >= 15)

  # Dynamics
  dynamics <- getMitoDynamicsGenes("mouse")

  # Check fission genes
  expect_true("Dnm1l" %in% dynamics)  # DRP1
  expect_true("Fis1" %in% dynamics)

  # Check fusion genes
  expect_true("Mfn1" %in% dynamics)
  expect_true("Mfn2" %in% dynamics)
  expect_true("Opa1" %in% dynamics)

  # Check mitophagy
  expect_true("Pink1" %in% dynamics)
  expect_true("Park2" %in% dynamics)  # Parkin
})

test_that("ATP metabolism genes are complete", {
  atp_genes <- getATPGenes("mouse")

  # Should include ATP synthase subunits
  expect_true(any(grepl("^Atp5", atp_genes)))

  # Should include major ATP consumers
  expect_true(any(grepl("^Atp1a", atp_genes)))  # Na/K-ATPase
  expect_true(any(grepl("^Atp2a", atp_genes)))  # SERCA

  expect_true(length(atp_genes) >= 40)
})

test_that("SVD signature includes all components", {
  svd_genes <- getSVDGenes("mouse")

  # Hypoxia genes
  expect_true("Hif1a" %in% svd_genes)
  expect_true("Vegfa" %in% svd_genes)

  # Vascular genes
  expect_true("Col4a1" %in% svd_genes)
  expect_true("Notch3" %in% svd_genes)

  # Inflammation
  expect_true("Il1b" %in% svd_genes)
  expect_true("Tnf" %in% svd_genes)

  # Oxidative stress
  expect_true("Sod1" %in% svd_genes)
  expect_true("Sod2" %in% svd_genes)

  expect_true(length(svd_genes) >= 80)
})

test_that("Pathway listing functions work", {
  # List pathways
  pathway_df <- listMetabolicPathways("mouse")
  expect_s3_class(pathway_df, "data.frame")
  expect_equal(ncol(pathway_df), 2)
  expect_equal(nrow(pathway_df), 6)
  expect_true(all(c("pathway", "n_genes") %in% colnames(pathway_df)))

  # Get all pathways
  all_paths <- getAllMetabolicPathways("mouse")
  expect_type(all_paths, "list")
  expect_equal(length(all_paths), 6)
  expect_true(all(c("OXPHOS", "Glycolysis", "SVD_Signature") %in% names(all_paths)))

  # Test subsetting
  subset_paths <- getAllMetabolicPathways("mouse", pathways = c("OXPHOS", "Glycolysis"))
  expect_equal(length(subset_paths), 2)
})

test_that("Pathway genes don't overlap inappropriately", {
  all_paths <- getAllMetabolicPathways("mouse")

  # OXPHOS and glycolysis should be mostly distinct
  oxphos <- all_paths$OXPHOS
  glycolysis <- all_paths$Glycolysis

  overlap <- intersect(oxphos, glycolysis)
  expect_true(length(overlap) / length(oxphos) < 0.1)  # Less than 10% overlap

  # Biogenesis should not overlap much with dynamics
  biogen <- all_paths$Mito_Biogenesis
  dynamics <- all_paths$Mito_Dynamics

  overlap_mito <- intersect(biogen, dynamics)
  expect_true(length(overlap_mito) < 5)  # Very few shared genes
})

test_that("Human-mouse gene conversion is consistent", {
  # Get both species
  oxphos_mouse <- getOXPHOSGenes("mouse")
  oxphos_human <- getOXPHOSGenes("human")

  # Should have similar numbers
  expect_true(abs(length(oxphos_mouse) - length(oxphos_human)) < 10)

  # Check case conversion pattern
  mouse_upper <- toupper(oxphos_mouse)

  # Many should match when converted
  matches <- sum(mouse_upper %in% oxphos_human)
  expect_true(matches / length(oxphos_mouse) > 0.8)  # >80% match
})

test_that("Pathway coverage calculation works", {
  # Create test object with known genes
  genes <- c(getOXPHOSGenes("mouse")[1:50],
             getGlycolysisGenes("mouse")[1:10],
             paste0("Unknown", 1:40))

  # Calculate coverage
  pathways <- getAllMetabolicPathways("mouse", c("OXPHOS", "Glycolysis"))

  coverage_oxphos <- sum(pathways$OXPHOS %in% genes) / length(pathways$OXPHOS)
  coverage_glyc <- sum(pathways$Glycolysis %in% genes) / length(pathways$Glycolysis)

  expect_true(coverage_oxphos > 0.3)  # At least 30% coverage
  expect_true(coverage_glyc > 0.4)    # At least 40% coverage
})
