context("fusionImage")

test_that("Ergas spat gives correct results with PCA", {
    pca <- pca_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
    res <- ergas_spat(original=mis, modified=pca, pan=pan, method="bilinear")
    expect_equal(round(res, 5), 0.74118)
})


test_that("Ergas spec gives correct results with PCA", {
    pca <- pca_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
    res <- ergas_spec(original=mis, modified=pca, pan=pan, method="bilinear")
    expect_equal(round(res, 5), 2.1473)
})


test_that("Ergas spat gives correct results with HPF", {
    hpf <- hpf_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
    res <- ergas_spat(original=mis, modified=hpf, pan=pan, method="bilinear")
    expect_equal(round(res, 5), 1.44371)
})


test_that("Ergas spec gives correct results with HPF", {
    hpf <- hpf_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
    res <- ergas_spec(original=mis, modified=hpf, pan=pan, method="bilinear")
    expect_equal(round(res, 5), 0.98474)
})


test_that("Ergas spat gives correct results with GS", {
    gs <- gs_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
    res <- ergas_spat(original=mis, modified=gs, pan=pan, method="bilinear")
    expect_equal(round(res, 5), 0.82001)
})


test_that("Ergas spec gives correct results with GS", {
    gs <-  gs_fusion(mis=mis, pan=pan, method="bilinear", bits=16)
    res <- ergas_spec(original=mis, modified=gs, pan=pan, method="bilinear")
    expect_equal(round(res, 5), 2.24995)
})
