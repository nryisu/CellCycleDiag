library(CellCycleDiag)
internal <- function(n) getFromNamespace(n, "CellCycleDiag")
equal <- function(a, b) stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-12)))
fails <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
d <- data.frame(group = rep(c("A", "B"), each = 4),
                S.Score = c(0,0,1,1,0,1,1,2),
                G2M.Score = c(0,0,0,0,1,1,1,1),
                theta = c(0,0,pi,pi,0,pi,pi,pi))
original <- d
sc <- calc_D_score(d, "group")
equal(sc$D_score_S, .25)
equal(sc$D_score_G2M, 1)
equal(sc$D_score, .625)
equal(calc_D_theta(d, "group")$D_theta, .25)
equal(calc_D_score(d, "group", group_levels = c("B", "A"))$D_score, sc$D_score)
equal(calc_D_theta(d, "group", group_levels = c("B", "A"))$D_theta, .25)
same <- rbind(d[1:4,], transform(d[1:4,], group = "B"))
equal(calc_D_score(same, "group")$D_score, 0)
equal(calc_D_theta(same, "group")$D_theta, 0)
wrap <- data.frame(group=c("A","A","B","B"), theta=c(0, pi, 2*pi, -pi))
equal(calc_D_theta(wrap, "group")$D_theta, 0)
boundary <- data.frame(group=c("A","B"), theta=c(0, pi/6))
equal(calc_D_theta(boundary,"group")$D_theta, 0) # right-closed bins
boundary$theta[2] <- pi/6 + 1e-10
equal(calc_D_theta(boundary,"group")$D_theta, 1)
rot <- d
rot$theta <- rot$theta + .3
equal(calc_D_theta(rot,"group",origin=.3)$D_theta, .25)
set.seed(7)
for (n in c(3, 10, 100)) {
  z <- data.frame(group=rep(c("A","B"),each=n), S.Score=sample(1:8,2*n,TRUE),
                  G2M.Score=rnorm(2*n), theta=runif(2*n,-8,8))
  ks <- suppressWarnings(stats::ks.test(z$S.Score[1:n], z$S.Score[n+1:n])$statistic)
  equal(calc_D_score(z,"group")$D_score_S, unname(ks))
  theta_bins <- cut(z$theta %% (2*pi),seq(0,2*pi,length.out=13),include.lowest=TRUE)
  p <- prop.table(table(z$group, theta_bins), 1)
  equal(calc_D_theta(z,"group")$D_theta, sum(abs(p[1,]-p[2,]))/2)
}
fails(calc_D_score(d, "absent"), "Missing")
fails(calc_D_score(d,"group",group_levels=c("A","A")), "distinct")
fails(calc_D_score(d,"group",group_levels=c("A","C")), "Each")
fails(calc_D_score(transform(d,S.Score=NA_real_),"group"), "finite")
fails(calc_D_theta(transform(d,theta=Inf),"group"), "finite")
fails(calc_D_theta(d,"group",n_bins=2.5), "range")
fails(calc_D_theta(d,"group",n_bins=1), "range")
fails(calc_D_theta(d,"group",origin=NA_real_), "range")
fails(calc_D_score(d,"group",s_col="S.Score",g2m_col="S.Score"), "distinct")
extra <- rbind(d, transform(d[1,,drop=FALSE], group=NA, theta=NA))
fails(calc_D_score(extra,"group"), "Missing group")
res <- diagnose_cc_imbalance(extra,"group","A","B")
stopifnot(inherits(res,"cc_diagnostic"), res$settings$excluded_cells==1)
equal(res$metrics$D_score, .625)
equal(d, original)
stopifnot(any(grepl("biologically required",capture.output(print(res)))))
f <- tempfile(fileext=".pdf")
grDevices::pdf(f)
out <- plot_cc_diagnostic(res)
grDevices::dev.off()
equal(out$value,c(.625,.25))
fails(plot_cc_diagnostic(res,"susceptibility"), "require")
cmp <- internal("compare_deg_sets")
v <- cmp(c("a","a","b"), c("b","c"))
equal(v$jaccard_distance,2/3)
equal(v$removed_fraction, .5)
stopifnot(is.na(cmp(character(),character())$jaccard_distance))
equal(cmp(character(),"a")$jaccard_distance,1)
stopifnot(is.na(cmp(character(),"a")$removed_fraction))
tab <- data.frame(avg_log2FC=c(.25,.26,-.3,.5),p_val_adj=c(.01,.049,.05,.01),
                  pct.1=c(.1,.1,.2,.099),pct.2=0,row.names=letters[1:4])
equal(internal(".sig_genes")(tab,.05,.25,.1),"b")
fails(internal(".sig_genes")(data.frame(),.05,.25,.1),"no tested")
fails(internal(".sig_genes")(transform(tab,p_val_adj=NA_real_),.05,.25,.1),"finite")
design <- data.frame(group=rep(c("A","B"),each=6), donor=rep(c("d1","d2","d3"),4),
                     score=seq_len(12)^2)
internal(".check_design")(design,"group",c("donor","score"),"test")
fails(internal(".check_design")(transform(design,donor=group),"group","donor","test"),"rank deficient")
fails(internal(".check_design")(transform(design,score=0),"group","score","test"),"Constant")
set.seed(12)
before <- .Random.seed
internal(".with_seed")(1,runif(3))
equal(.Random.seed,before)
fails(internal(".with_seed")(1,stop("deliberate")),"deliberate")
equal(.Random.seed,before)
calls <- list()
mock <- function(...) { calls[[length(calls)+1L]] <<- list(...); tab }
models <- list(direct="donor",score_adjusted=c("donor","S.Score","G2M.Score"),
               cyclic_adjusted=c("donor","sin","cos"))
got <- internal(".fit_three")(d,"group","A","B","RNA",models,.1,NULL,1,FALSE,mock)
equal(names(got),names(models))
for (i in 1:3) {
  equal(calls[[i]]$latent.vars,models[[i]])
  equal(calls[[i]]$object,d)
  stopifnot(calls[[i]]$base==2,calls[[i]]$logfc.threshold==0,calls[[i]]$test.use=="MAST")
}
fails(internal(".fit_three")(d,"group","A","B","RNA",models,.1,NULL,1,FALSE,
  function(...) stop("backend unavailable")),"direct model failed")
cat("Core numeric, edge-case, model-dispatch, plotting and RNG tests passed.\n")
