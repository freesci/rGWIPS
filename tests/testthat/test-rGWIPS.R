

context("Loading GWIPS data")

test_that("loadGWIPSdata function fails on wrong file paths",{

  expect_error(loadGWIPSdata("wrongpath", system.file(package="rGWIPS", "extdata","Liu13_All.RiboProElongRev.bw")),
               "file.exists(forward) is not TRUE", fixed=TRUE)

  expect_error(loadGWIPSdata(system.file(package="rGWIPS", "extdata","Liu13_All.RiboProElongRev.bw"),"wrongpath"),
               "file.exists(reverse) is not TRUE", fixed=TRUE)

})
