library(mvtnorm)
set.seed(123)

# generating a covariance matrix
covariance <- matrix(c(5,3,0,-3,0, 3,5,0,-3,0, 0,0,5,0,0,
    -3,-3,0,6,0, 0,0,0,0,3), nrow=5)
covariance
    # row1 col 2 is 3- variables 1 and 2 have a cov of 3
    # row 1 col 1- variance of variable 1 is 5 (cov with yourself)

data <- rmvnorm(1000, sigma = covariance)
    # note thediff variables need the same variances for pca to work
    # so you must z transform them at this point

colnames(data) <- c("a", "b", "c", "d", "e")
colnames(data)

head(data)

# run pca
pca <- prcomp(data)
biplot(pca)
biplot(pca, choices = 2:3)
pca # the loadings, i.e. the correlations between the pc axis and the variable
    # pc1 is correlated with a and b and negatively with d
        # can see this on arrow directions on main biplot ^

plot(pca$x[,2] ~ data[,3], xlab="'c' variable", ylab="PC2")

summary(pca) # sds of each axis, i.e. how important is each axis
    # want a high sd as means axis explains a lot of variance in data

# how many axes should we include? Hard to answer w PCA- do a scree plot
plot(pca)
    # find inflection point in sds- when it bends
    # take the straight line plus the next one- here, take pca1-3 (?)

pca <- prcomp(data, scale = TRUE)
pca