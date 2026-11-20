test_that("Seurat and MAST complete all three models without mutating input", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("MAST")
  set.seed(42)
  counts <- matrix(rpois(300*80,lambda=3),nrow=300,
                   dimnames=list(paste0("gene",1:300),paste0("cell",1:80)))
  counts[1:10,41:80] <- counts[1:10,41:80]+10
  obj <- Seurat::CreateSeuratObject(counts)
  obj <- Seurat::NormalizeData(obj,verbose=FALSE)
  obj$condition <- rep(c("A","B"),each=40)
  obj$donor <- rep(paste0("d",1:4),20)
  obj$S.Score <- rnorm(80)
  obj$G2M.Score <- rnorm(80)
  obj$theta <- runif(80,0,2*pi)
  before <- obj
  ans <- run_cc_sensitivity(obj,"condition","A","B",donor.var="donor")
  expect_s3_class(ans,"cc_sensitivity")
  expect_identical(obj,before)
  expect_equal(nrow(ans$susceptibility),2)
  expect_true(all(vapply(ans$models,function(v) "donor" %in% v,logical(1))))
  manual <- Seurat::FindMarkers(obj,group.by="condition",ident.1="A",ident.2="B",
    assay="RNA",test.use="MAST",latent.vars="donor",min.pct=.1,
    logfc.threshold=0,base=2,verbose=FALSE)
  expect_equal(ans$deg_tables$direct,manual)
})

test_that("tricycle representation uses cell names and preserves identities", {
  for (p in c("Seurat","tricycle","SingleCellExperiment","SummarizedExperiment"))
    skip_if_not_installed(p)
  e <- new.env()
  utils::data("neurosphere_example",package="tricycle",envir=e)
  sce <- e$neurosphere_example
  counts <- SummarizedExperiment::assay(sce,"counts")
  symbols <- as.character(SummarizedExperiment::rowData(sce)$Gene)
  keep <- !is.na(symbols) & nzchar(symbols) & !duplicated(symbols)
  counts <- counts[keep,,drop=FALSE]
  rownames(counts) <- symbols[keep]
  obj <- Seurat::NormalizeData(Seurat::CreateSeuratObject(counts),verbose=FALSE)
  genes <- new.env()
  utils::data("cc.genes.updated.2019",package="Seurat",envir=genes)
  sf <- rownames(obj)[toupper(rownames(obj)) %in% genes$cc.genes.updated.2019$s.genes]
  gf <- rownames(obj)[toupper(rownames(obj)) %in% genes$cc.genes.updated.2019$g2m.genes]
  before <- obj
  ans <- add_cc_representation(obj,species="mouse",s.features=sf,g2m.features=gf)
  expect_identical(obj,before)
  expect_identical(SeuratObject::Idents(ans),SeuratObject::Idents(obj))
  expect_equal(ans$cc_sin,sin(ans$theta))
  expect_true(all(is.finite(ans$theta)))
  expect_error(add_cc_representation(ans,species="mouse"),"Existing")
})
