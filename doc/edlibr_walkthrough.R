## -----------------------------------------------------------------------------
library(edlibR)

algn1 = align("ACTG", "CACTRT", mode="HW", task="path")
print(algn1)

## -----------------------------------------------------------------------------
algn2 = align("elephant", "telephone")
print(algn2)

## -----------------------------------------------------------------------------
algn3 = align("ACTG", "CACTRT", mode="HW", task="path")
print(algn3)

## -----------------------------------------------------------------------------
## the previous example with additionalEqualities 
algn4 = align("ACTG", "CACTRT", mode="HW", task="path", additionalEqualities=list(c("R", "A"), c("R", "G")))
print(algn4)

## -----------------------------------------------------------------------------
library(edlibR)

query = "elephant"
target = "telephone"
result = align(query, target, task = "path")
nice_algn = getNiceAlignment(result, query, target)
print(nice_algn)

## -----------------------------------------------------------------------------
library(edlibR)
## example above from getNiceAlignment()

query = "elephant"
target = "telephone"
result = align(query, target, task = "path")
nice_algn = getNiceAlignment(result, query, target)
nice_print(nice_algn)

