test_that("simu_db() returns a data.frame", {
  simu_db() %>% expect_s3_class('data.frame')
})

test_that("Column names are as exepected", {
  simu_db(multi_imp = T) %>%
    expect_named(c('Peptide', 'Group', 'Sample', 'Output', 'Draw'),
                 ignore.order = T)
})


