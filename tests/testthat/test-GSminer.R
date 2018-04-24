library(GSminer)
library(testthat)
context("test removeDup")

 

test_that("remove Duplicate rows", {
	testdata <- data.frame("a" = c("ab", "cd", "ef"), "b" = c("cd", "ab", "ef"), stringsAsFactors = FALSE) 
	result <- data.frame("a" = c("ab", "ef"), "b" = c("cd", "ef"), stringsAsFactors = FALSE) 
    expect_true(all.equal(result, removeDup(testdata),  check.attributes = F))
})

