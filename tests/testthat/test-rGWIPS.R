

context("Loading GWIPS data")

test_that("loadGWIPSdata function correctly loads data",{

  expect_error(loadGWIPSdata("wrongpath", system.file(package="rGWIPS", "extdata","Liu13_All.RiboProElongRev.bw")),
               "file.exists(forward) is not TRUE", fixed=TRUE)

  expect_error(loadGWIPSdata(system.file(package="rGWIPS", "extdata","Liu13_All.RiboProElongRev.bw"),"wrongpath"),
               "file.exists(reverse) is not TRUE", fixed=TRUE)

  for_bw = system.file(package = "rGWIPS", "extdata", "Liu13_All.RiboProElongFor.bw")
  rev_bw = system.file(package = "rGWIPS", "extdata", "Liu13_All.RiboProElongRev.bw")

  loadGWIPSdata(forward = for_bw, reverse = rev_bw)

  expect_is(gwips_forw, "GRanges")
  expect_is(gwips_rev, "GRanges")

})
