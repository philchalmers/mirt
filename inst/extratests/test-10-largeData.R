#
# "largeData" context is a set of tests checking consistency of item parameters
# and factor scores estimetes obtained in mirt and the other software 
# (currently Mplus only) using big real-world datasets
#
context('largeData')

test_that('large1dim', {
  require(plyr)
  mirtCluster(3)
  
  data <- dataComplete[, grep('^i[0-9]', names(dataComplete))]
  itemTypes = ifelse(
    sapply(data, max) > 1,
    'graded', 
    '2PL'
  )
  model = mirt(data, 1, itemtype=itemTypes, se=T)
  
  mplusParams = paramComplete # are part of data/dataComplete.RData
  mirtParams = mod2values(model)[, c('item', 'name', 'value')]
  mirtParams$item = as.character(mirtParams$item)
  mirtParams$name = as.character(mirtParams$name)
  mirtParams$name = sub('^d$', 'd1', mirtParams$name)
  params = join(mplusParams, mirtParams, by=c('item', 'name'), type='left', match='first')
  for(n in c('^a1$', '^d')){
    filter = grepl(n, params$name)
    comp = abs(cor(params[filter, 3], params[filter, 5]))
    plot(params[filter, 3], params[filter, 5], main=paste0(n, ': cor(', round(comp, 3), ')'))
    expect_equal(comp, 1, tolerance = 0.05)
  }

  mirtScores = fscores(model, full.scores=T, method='EAP')
  comp = cor(mirtScores[, 1], dataComplete$theta)
  plot(mirtScores[, 1], dataComplete$theta)
  expect_equal(comp, 1, tolerance = 0.01)  
})

test_that('sparse1dim', {
  require(plyr)
  mirtCluster(3)
  
  data <- dataSparse[, grep('^i[0-9]', names(dataSparse))]
  itemTypes = ifelse(
    sapply(data, max, na.rm=T) > 1,
    'graded', 
    '2PL'
  )
  model = mirt(data, 1, itemtype=itemTypes, se=T)
  
  mplusParams = paramSparse # are part of data/dataSparse.RData
  mirtParams = mod2values(model)[, c('item', 'name', 'value')]
  mirtParams$item = as.character(mirtParams$item)
  mirtParams$name = as.character(mirtParams$name)
  mirtParams$name = sub('^d$', 'd1', mirtParams$name)
  params = join(mplusParams, mirtParams, by=c('item', 'name'), type='left', match='first')
  for(n in c('^a1$', '^d')){
    filter = grepl(n, params$name)
    comp = abs(cor(params[filter, 3], params[filter, 5]))
    plot(params[filter, 3], params[filter, 5], main=paste0(n, ': cor(', round(comp, 3), ')'))
    expect_equal(comp, 1, tolerance = 0.05)
  }
  
  mirtScores = fscores(model, full.scores=T, method='EAP')
  comp = cor(mirtScores[, 1], dataSparse$theta)
  plot(mirtScores[, 1], dataSparse$theta)
  expect_equal(comp, 1, tolerance = 0.01)  
})

test_that('large2dim', {
  require(plyr)
  mirtCluster(3)
  
  data <- data2Dim[, grep('^[ij][0-9]', names(data2Dim))]
  itemTypes = ifelse(
    sapply(data, max, na.rm=T) > 1,
    'graded', 
    '2PL'
  )
  modelDef = mirt.model('
    theta1 = 1-12
    theta2 = 13-52
    COV = theta1 * theta2
  ')
  model = mirt(data, modelDef, itemtype=itemTypes, se=T)
  
  mplusParams = param2Dim # are part of data/dataSparse.RData
  mirtParams = mod2values(model)[, c('item', 'name', 'value')]
  mirtParams$item = as.character(mirtParams$item)
  mirtParams$name = as.character(mirtParams$name)
  mirtParams$name = sub('^d$', 'd1', mirtParams$name)
  params = join(mplusParams, mirtParams, by=c('item', 'name'), type='left', match='first')
  for(n in c('^a1$', '^d')){
    filter = grepl(n, params$name)
    comp = abs(cor(params[filter, 3], params[filter, 5]))
    plot(params[filter, 3], params[filter, 5], main=paste0(n, ': cor(', round(comp, 3), ')'))
    expect_equal(comp, 1, tolerance = 0.05)
  }
  filter = grepl('COV_21', params$name)
  expect_equal(params[filter, 3], params[filter, 5], tolerance = 0.01)
  
  mirtScores = fscores(model, full.scores=T, method='EAP')
  
  comp = cor(mirtScores[, 1], data2Dim$theta1)
  plot(mirtScores[, 1], data2Dim$theta1)
  expect_equal(comp, 1, tolerance = 0.01)  

  comp = cor(mirtScores[, 2], data2Dim$theta2)
  plot(mirtScores[, 2], data2Dim$theta2)
  expect_equal(comp, 1, tolerance = 0.01)  
})