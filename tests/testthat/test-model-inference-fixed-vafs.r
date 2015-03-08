# Test model inference using fixed VAFs of the cluster, no variant list with
# their VAFs provided

library(testthat)
library(clonevol)

#test_check("clonevol")

context('Test model inference using fixed VAFs of the clusters')

test_that('*** One sample with 100% purity, two clones of same 0.5 VAF,
expect to have 2 models in mono and 2 in polyclonal\n', {
  clones = data.frame(cluster=c(1,2,3,4,5,6), primary=c(0.5, 0, 0.5, 0, 0.3, 0))
  models = c('monoclonal', 'polyclonal')
  num.true.models = as.integer(c(2, 2))
  names(num.true.models) = models
  expect_identical(run.test(clones, models), num.true.models)
})


test_that('*** Test 0. 1 sample with 80% purity, 1 clone, expect 1 model
in both mono and polyclonal\n',{
  clones = data.frame(cluster=c(1), primary=c(0.4))
  num.true.models = as.integer(c(1, 1))
  models = c('monoclonal', 'polyclonal')
  names(num.true.models) = models
  expect_identical(run.test(clones, models), num.true.models)
})

test_that('*** Test 1. 1 sample with 100% purity, expect to have 2 models
in both mono and polyclonal\n',{
  clones = data.frame(cluster=c(2,3,4), primary=c(0.5,0.1,0.05))
  num.true.models = c('monoclonal'=2, 'polyclonal'=2)
  num.true.models = as.integer(num.true.models)
  models = c('monoclonal', 'polyclonal')
  names(num.true.models) = models
  expect_identical(run.test(clones, models), num.true.models)
})


test_that('*** Test 2. 1 sample with 70% purity, expect to have 2 models
in mono and 6 in polyclonal\n',{
  clones = data.frame(cluster=c(1,2,3), primary=c(0.35,0.1,0.05))
  num.true.models = c('monoclonal'=2, 'polyclonal'=6)
  num.true.models = as.integer(num.true.models)
  models = c('monoclonal', 'polyclonal')
  names(num.true.models) = models
  expect_identical(run.test(clones, models), num.true.models)
})

test_that('*** Test 3. 1 sample with 100% purity, two clones of same VAF,
expect to have 3 models in mono and 6 in polyclonal\n',{
  clones = data.frame(cluster=c(1,2,3), primary=c(0.5,0.2,0.2))
  num.true.models = c('monoclonal'=3, 'polyclonal'=3)
  num.true.models = as.integer(num.true.models)
  models = c('monoclonal', 'polyclonal')
  names(num.true.models) = models
  expect_identical(run.test(clones, models), num.true.models)
})


test_that('*** Test 4. 3 samples with 100% purity
expect to have 1 model in mono and 1 in polyclonal\n',{
  clones = threeSamples
  num.true.models = c('monoclonal'=1, 'polyclonal'=1)
  num.true.models = as.integer(num.true.models)
  models = c('monoclonal', 'polyclonal')
  names(num.true.models) = models
  expect_identical(run.test(clones, models), num.true.models)
})


test_that('*** Test 5. AML31 primary/relapse samples with VAF scaled to 0.5 (100% purity),
expect 5 models in mono and 5 in polyclonal\n',{
  clones = aml31
  num.true.models = c('monoclonal'=5, 'polyclonal'=5)
  num.true.models = as.integer(num.true.models)
  models = c('monoclonal', 'polyclonal')
  names(num.true.models) = models
  expect_identical(run.test(clones, models), num.true.models)
})


test_that('*** Test 6. AML31 primary/relapse samples with 90% purity in primary,
and 40% purity in relapse. Expect 5 models in mono and 37 in polyclonal\n',{
  clones = aml31
  clones$primary = clones$primary*0.9
  clones$relapse = clones$relapse*0.4
  num.true.models = c('monoclonal'=5, 'polyclonal'=37)
  num.true.models = as.integer(num.true.models)
  models = c('monoclonal', 'polyclonal')
  names(num.true.models) = models
  expect_identical(run.test(clones, models), num.true.models)
})

