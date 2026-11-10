# single factor
set.seed(12345)
a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d <- matrix(rnorm(15,0,.7),ncol=1)
itemtype <- rep('2PL', nrow(a))
N <- 1000
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
data <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))

# test plot
empirical_plot(data)
empirical_plot(data, group=group)

empirical_plot(data, type = 'cor')
empirical_plot(data, type = 'cor', group=group)

empirical_plot(data, type = 'hist')
empirical_plot(data, type = 'hist', group=group)

empirical_plot(data, type = 'hist', ts_fun=rowMeans)
empirical_plot(data, type = 'hist', ts_fun=rowMeans, group=group)

empirical_plot(data, type = 'hist', breaks=20)
empirical_plot(data, type = 'hist', breaks=20, group=group)

empirical_plot(data, type = 'bar')
empirical_plot(data, type = 'bar', group=group)

empirical_plot(data, type = 'discrim')
empirical_plot(data, type = 'discrim', group=group)

empirical_plot(data, type = 'discrim', sort=TRUE)
empirical_plot(data, type = 'discrim', sort=TRUE, group=group)

empirical_plot(data, type = 'difficulty')
empirical_plot(data, type = 'difficulty', group=group)

empirical_plot(data, type = 'difficulty', sort=TRUE)
empirical_plot(data, type = 'difficulty', sort=TRUE, group=group)

empirical_plot(data, type = 'discrim_diff', group=group)

empirical_plot(data, type = 'alpha_rm')
empirical_plot(data, type = 'alpha_rm', group=group)

empirical_plot(data, type = 'freq')
empirical_plot(data, type = 'freq', group=group)
