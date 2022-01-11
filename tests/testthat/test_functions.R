
library(edlibR)
library(testthat)


## from https://pypi.org/project/edlib/
## """
## import edlib
## edlib.align("ACTG", "CACTRT", mode="HW", task="path")
## # {'locations': [(1, 3), (1, 4)], 'cigar': '3=1I', 'alphabetLength': 5, 'editDistance': 1}
## """

test_that("basic align() functionality, query=ACTG, target=CACTRT", {
	algn1 = align("ACTG", "CACTRT", mode="HW", task="path")
	expect_true(is.list(algn1))
	expect_equal(length(algn1), 5)
	expect_equal(algn1$editDistance, 1)
	expect_equal(algn1$alphabetLength, 5)
	expect_equal(length(algn1$locations), 2)
	expect_equal(algn1$cigar, "3=1I")
	expect_equal(algn1$cigarFormat, "extended")
})

## from https://pypi.org/project/edlib/
## """
## import edlib
## edlib.align("elephant", "telephone")
## # {'locations': [(None, 8)], 'cigar': None, 'alphabetLength': 8, 'editDistance': 3}
## """

test_that("basic align() functionality, query=elephant, target=telephone", {
	expect_equal(length(align("elephant", "telephone")), 5)
	expect_equal(align("elephant", "telephone")$editDistance, 3)
	expect_equal(align("elephant", "telephone")$alphabetLength, 8)
	expect_equal(length(align("elephant", "telephone")$locations), 1)
	expect_true(is.null(align("elephant", "telephone")$cigar))
	expect_equal(align("elephant", "telephone")$cigarFormat, "extended")
})

## Check additionalEqualities functionality

test_that("Check additionalEqualities functionality", {
	algn2 = align("ACTG", "CACTRT", mode="HW", task="path", additionalEqualities=list(c("R", "A"), c("R", "G")))
	expect_true(is.list(algn2))
	expect_equal(length(algn2), 5)
	expect_equal(algn2$editDistance, 0)
	expect_equal(algn2$alphabetLength, 5)
	expect_equal(algn2$locations[[1]][1], 1)
	expect_equal(algn2$locations[[1]][2], 4)
	expect_equal(algn2$cigar, "4=")
	expect_equal(algn2$cigarFormat, "extended")
})

## Check locations working correctly
test_that("Check locations working correctly", {
	algn3 = align("AACG", "TCAACCTG", mode = "HW", task = "path")
	expect_true(is.list(algn3))
	expect_equal(length(algn3), 5)
	expect_equal(algn3$editDistance, 1)
	expect_equal(algn3$alphabetLength, 4)
	expect_equal(algn3$locations[[1]][1], 2)
	expect_equal(algn3$locations[[1]][2], 4)
	expect_equal(algn3$locations[[2]][1], 2)
	expect_equal(algn3$locations[[2]][2], 5)
	expect_equal(algn3$cigar, "3=1I")
	expect_equal(algn3$cigarFormat, "extended")
})


## Check getNiceAlignment()
test_that("Check getNiceAlignment()", {
	query = "elephant"
	target = "telephone"
	result = align(query, target, task = "path")
	nice = getNiceAlignment(result, query, target)
	expect_equal(nice$query_aligned, "-elephant")
	expect_equal(nice$matched_aligned, "-|||||.|.")
	expect_equal(nice$target_aligned, "telephone")
})




