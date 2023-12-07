################################################################
############# Define functions needed for analysis #############
################################################################

# functions:
# squared exp kernel:
K <- function(xi, xj, sigma_f, l) {
  return((sigma_f^2)*exp((-1/(2*l^2))*((xi-xj)^2) ))
}
# partial kernel wrt l
K_prime_l <- function(xi, xj, sigma_f, l) {
  return((sigma_f^2)*exp((-1/(2*l^2))*( (xi-xj)^2) )*( ( (xi-xj)^2)/(l^3) ) 
  )
}
# partial kernel wrt sigma
K_prime_sigma <- function(xi, xj, sigma_f, l){
  return(2*sigma_f*exp( (-1/(2*l^2))*((xi-xj)^2) ) )
} 
# compute covariance matrix
get_sigmad <- function(K, x, h, l){
  Sigma <- sapply(x, function(s1) {
    sapply(x, function(s2) {
      K(s1, s2, h, l)
    })
  })
  return(Sigma)
}
# find gradient for l
gradient <- function(l, t, path, sigma_f){
  
  Sigma <- get_sigmad(K, t, sigma_f, l)
  
  Sigma_prime_l <- get_sigmad(K_prime_l, t, sigma_f, l)
  
  K_inverse <- FastGP::rcppeigen_invert_matrix(Sigma)
  # alpha is K^(-1)y as defined in (27 environmetrics)
  alpha <- K_inverse%*%y
  gradient <- (1/2)*psych::tr(  ( ( alpha %*% t(alpha) ) - 
                                    K_inverse) %*% Sigma_prime_l )
  return(gradient)
}

get_gradiented_l <- function(y, l, x){
  y <- y
  x <- x
  L_1 <- l
  S_1 <- 1
  N <- 12
  continue <- TRUE
  while(continue){
    # old <- dmvnorm(y, mean=rep(0, N), 
    #                get_sigmad(K, x, S_1, L_1),
    #                log=T)
    # compute log likelihood:
    Sigma <- get_sigmad(K, x, S_1, L_1)
    K_inverse <- FastGP::rcppeigen_invert_matrix(Sigma)
    old <- (-1/2)*(y %*% K_inverse %*% y + log(det(Sigma)) + N*log(2*pi))
    # 
    # temp variables:
    move_by_l <- gradient(L_1, x, y, S_1) 
    L_1 <- L_1 + .001*move_by_l
    
    
    
    
    
    # new <- dmvnorm(y, mean = rep(0, N), 
    #                get_sigmad(K, x, S_1, L_1),
    #                log = T)
    Sigma <- get_sigmad(K, x, S_1, L_1)
    K_inverse <- FastGP::rcppeigen_invert_matrix(Sigma)
    new <-  (-1/2)*(y %*% K_inverse %*% y + log(det(Sigma)) + N*log(2*pi))
    
    #print(new)
    continue <- abs(new-old) > .00000001
  }
  return(c(L_1, S_1, new) ) 
  
}


# log like
log_like_mixture <- function(data_1, N, K, t, l_1, l_2, pi_1, pi_2, s_1, s_2){
  # calculates log likelihood
  # need dataframe, 
  # N=number of realizations, 
  # K=kernel, and parameters
  pi_1 <- pi_1
  pi_1 <- pi_1
  l_1  <- l_1
  l_2  <- l_2
  s_1  <- s_1
  s_2  <- s_2
  
  like_1 <- apply(data_1, MARGIN = 2, dmvnorm, mean = rep(0, N), 
                  sigma = get_sigmad(K, t, s_1, l_1)) 
  
  like_2 <- apply(data_1, MARGIN = 2, dmvnorm, mean = rep(0, N), 
                  sigma = get_sigmad(K, t, s_2, l_2)) 
  # checked: log likelihood
  zz <- sum(log(pi_1*like_1 + pi_2*like_2))
  return(zz)
}


################################################################
############# Calculate Optimal Ls #############################
################################################################

years1 <- c("1991","1992", "1993", "1994", 
            "1995","1996","1997", "1998", "1999", "2000")
pacman::p_load(dplyr, tidyr, stringr, ggplot2, lubridate, mvtnorm, latex2exp, 
               FastGP, psych, NMOF, gtools, ggrepel)

tofino <- read.csv( "./data/Tofino_A.csv") %>% 
  dplyr::select(Date.Time, Year, Month, Total.Rain..mm.) %>% 
  mutate(Date.Time = paste0(Date.Time, "-1")) %>% 
  mutate(Date.Time = myd(Date.Time) ) %>% 
  mutate(el_nino_season = ifelse(Month %in% c(10, 11, 12, 1, 2, 3), 1, 0))

gp_data_rain <- tofino %>% 
  filter(Year %in% years1) %>% 
  dplyr::select(Total.Rain..mm., Year, Month) %>% 
  spread(Year, Total.Rain..mm.) %>% dplyr::select(-Month) 


tofino_seasoned <- (gp_data_rain - rowMeans(gp_data_rain)) %>% 
  scale() %>% as.data.frame()

# tofino l_opt
tofino_l_opt <-  c("Year","L_start", "L_final", "S", "Likelihood")
L = seq(from=1, to=2, by=0.01)
# for (year in colnames(tofino_seasoned)) {
#   for(i in L){
#     y <- tofino_seasoned[,year]
#     x <- x
#     L_1 <- i
#     S_1 <- 1
#     res = c()
#     tryCatch(
#       expr={
#         res = get_gradiented_l(y,L_1,x)
#       }, 
#       error = function(e){
#         print(e)
#         res <<- c(NA, 1, NA)
#       }
#     )
#     tofino_l_opt <-  rbind(tofino_l_opt,c(year, L_1, res))
#   }
# }

# loads in tofino_l_opt
load("./L_opts/tofino_l_opt.RData")

# vancouver l_opt

van <- read.csv( "./Vancouver_International_Airport.csv") %>% 
  dplyr::select(Date.Time, Year, Month, Total.Rain..mm.) %>% 
  mutate(Date.Time = paste0(Date.Time, "-1")) %>% 
  mutate(Date.Time = myd(Date.Time) ) %>% 
  filter(Year %in% years1) %>% 
  dplyr::select(Total.Rain..mm., Year, Month) %>% 
  spread(Year, Total.Rain..mm.) %>% dplyr::select(-Month) 

van_seasoned <-  (van - rowMeans(van)) %>% 
  scale() %>% as.data.frame()

x <- 1:12
van_l_opt <-  c("Year","L_start", "L_final", "S", "Likelihood")
L = seq(from=0, to=2, by=0.01)
# for (year in colnames(van_seasoned)) {
#   print(paste("Year:",year))
#   for(i in L){
#     print(paste0((round(i/length(L),4))*100,"% done"))
#     y <- van_seasoned[,year]
#     x <- x
#     L_1 <- i
#     S_1 <- 1
#     res = c()
#     tryCatch(
#       expr={
#         res = get_gradiented_l(y,L_1,x)
#       }, 
#       error = function(e){
#         print(e)
#         res <<- c(NA, 1, NA)
#       }
#     )
#     van_l_opt <-  rbind(van_l_opt,c(year, L_1, res))
#   }
# }

# loads in van_l_opt
load("./L_opts/van_l_opt.RData")


# port hardy l_opt

port_h <- read.csv( "./Port_Hardy_A.csv") %>% 
  dplyr::select(Date.Time, Year, Month, Total.Rain..mm.) %>% 
  mutate(Date.Time = paste0(Date.Time, "-1")) %>% 
  mutate(Date.Time = myd(Date.Time) ) %>% 
  filter(Year %in% years1) %>% 
  dplyr::select(Total.Rain..mm., Year, Month) %>% 
  spread(Year, Total.Rain..mm.) %>% dplyr::select(-Month)

port_seasoned <- (port_h - rowMeans(port_h)) %>% 
  scale() %>% as.data.frame()

x <- 1:12
ph_l_opt <-  c("Year","L_start", "L_final", "S", "Likelihood")
L = seq(from=0, to=2, by=0.01)
# for (year in colnames(port_seasoned)) {
#   print(paste("Year:",year))
#   for(i in L){
#     print(paste0((round(i/length(L),4))*100,"% done"))
#     y <- port_seasoned[,year]
#     x <- x
#     L_1 <- i
#     S_1 <- 1
#     res = c()
#     tryCatch(
#       expr={
#         res = get_gradiented_l(y,L_1,x)
#       }, 
#       error = function(e){
#         print(e)
#         res <<- c(NA, 1, NA)
#       }
#     )
#     ph_l_opt <-  rbind(ph_l_opt,c(year, L_1, res))
#   }
# }

# loads in ph_l_opt
load("./L_opts/ph_l_opt.RData")

# victoria l_opt

victoria <- read.csv( "data/Victoria_International_Airport.csv") %>% 
  dplyr::select(Date.Time, Year, Month, Total.Rain..mm.) %>% 
  mutate(Date.Time = paste0(Date.Time, "-1")) %>% 
  mutate(Date.Time = myd(Date.Time) ) %>% 
  filter(Year %in% years1) %>% 
  dplyr::select(Total.Rain..mm., Year, Month) %>% 
  spread(Year, Total.Rain..mm.) %>% dplyr::select(-Month) 

vic_seasoned <- (victoria - rowMeans(victoria)) %>% 
  scale() %>% as.data.frame()

# vic_l_opt <-  c("Year","L_start", "L_final", "S", "Likelihood")
# L = seq(from=0, to=2, by=0.01)
# for (year in colnames(vic_seasoned)) {
#    print(paste("Year:",year))
#   for(i in L){
#     print(paste0((round(i/length(L),4))*100,"% done"))
#     y <- vic_seasoned[,year]
#     x <- x
#     L_1 <- i
#     S_1 <- 1
#     res = c()
#     tryCatch(
#       expr={
#         res = get_gradiented_l(y,L_1,x)
#       }, 
#       error = function(e){
#         print(e)
#         res <<- c(NA, 1, NA)
#       }
#     )
#     vic_l_opt <-  rbind(vic_l_opt,c(year, L_1, res))
#   }
# }

# loads in vic_l_opt
load("./L_opts/vic_l_opt.RData")












































van_l_opt 


port_l_opt 


vic_l_opt

# load optimized HP l scales:
# load("optimizedhp.RData")

data_1 <- tofino_seasoned
# set dimensions:
N <- nrow(data_1)
t <- 1:N
# optimize GP hyperparameters:
ML_solutions <- tofino_l_opt %>% as.data.frame()

pi_1 = .9
pi_2 = .1

l_1 = .1
l_2 = .1
s_1 <- 1
s_2 <- 1
continue <- T
counter_ <- 1
log_like_mixture(data_1, N, K, t, l_1, l_2, pi_1, pi_2, s_1, s_2)
while (continue) {
  # contsruct sigma matrix
  Sigma_1 <- get_sigmad(K, t, s_1, l_1)
  Sigma_2 <- get_sigmad(K, t, s_2, l_2)
  # compute log likelihood:
  old <- log_like_mixture(data_1, N, K, t, l_1, l_2, pi_1, pi_2, s_1, s_2)
  
  # cluster 1 responsibilities
  r_j1 <- 
    (pi_1*apply(data_1, MARGIN = 2,  dmvnorm, mean = rep(0, N), sigma = Sigma_1, log=F))/
    (pi_2*apply(data_1, MARGIN = 2,  dmvnorm, mean = rep(0, N), sigma = Sigma_2, log=F) +
       pi_1*apply(data_1, MARGIN = 2,  dmvnorm, mean = rep(0, N), sigma = Sigma_1, log=F))
  
  # cluster 2 responsibilities
  r_j2 <- 
    (pi_2*apply(data_1, MARGIN = 2,  dmvnorm, mean = rep(0, N), sigma = Sigma_2, log=F))/
    (pi_2*apply(data_1, MARGIN = 2,  dmvnorm, mean = rep(0, N), sigma = Sigma_2, log=F)  +
       pi_1*apply(data_1, MARGIN = 2,  dmvnorm, mean = rep(0, N), sigma = Sigma_1, log=F))
  
  # calculate m_g, m
  m_g <- 1:2
  m_g[1] <- sum(r_j1)
  m_g[2] <- sum(r_j2)
  
  # calculate pi_g
  pi_1 <- m_g[1]/sum(m_g)
  pi_2 <- m_g[2]/sum(m_g)
  
  # update parameters
  # length parameter
  l_1 <- sum(r_j1*ML_solutions$.)/m_g[1]
  l_2 <- sum(r_j2*ML_solutions$.)/m_g[2]
  
  # update sigma_f
  s_1 <- 1
  s_2 <- 1
  
  # new likelihood:
  updated_log_like <- log_like_mixture(data_1, N, K, t, l_1, l_2, pi_1, pi_2, 
                                       s_1, s_2)
  # Print output:
  print(paste(c("Iteration #", counter_,
                " log-like:",
                round(updated_log_like, digits=5), " pi 1:", round(pi_1, digits=3)), collapse = "") )
  counter_ = counter_ + 1
  # continue:
  continue <- ifelse( (abs(old - updated_log_like) > .00001) & 
updated_log_like > old, T, F)
}


round(r_j1, digits=2)
l_1
l_2
pi_1
pi_2
ML_solutions$.




############# Run full search:
tofino_l_opt 
van_l_opt 
port_l_opt 
vic_l_opt

# data_1 <- van_seasoned
data_1 <- tofino_seasoned
# set dimensions:
N <- nrow(data_1)
t <- 1:N
# optimize GP hyperparameters:
ML_solutions <- tofino_l_opt   %>% as.data.frame()

# create one hot encoding matrix:
resp_matrix <- permutations(2,10,0:1,repeats.allowed=TRUE) 
# intialize vector:
dat_log_like <- 1:nrow(resp_matrix)
# loop through values:
for (i in 1:nrow(resp_matrix) ) {
  # cluster 1 responsibilities
  r_j1 <- resp_matrix[i, ]
  # cluster 2 responsibilities
  r_j2 <- 1-r_j1
  # calculate m_g, m
  m_g <- 1:2
  m_g[1] <- sum(r_j1)
  m_g[2] <- sum(r_j2)
  
  # calculate pi_g
  pi_1 <- m_g[1]/sum(m_g)
  pi_2 <- m_g[2]/sum(m_g)
  # update parameters
  # length parameter
  l_1 <- sum(r_j1*ML_solutions$.)/m_g[1]
  l_2 <- sum(r_j2*ML_solutions$.)/m_g[2]
  # update sigma_f
  s_1 <- 1 #sum(r_j1*ML_solutions$V2)/m_g[1]
  s_2 <- 1 #sum(r_j2*ML_solutions$V2)/m_g[2]
  # new likelihood:
  dat_log_like[i] <- log_like_mixture(data_1, N, K, t, l_1, l_2, pi_1, pi_2, 
                                      s_1, s_2)
}
# store values:
final_data_tofino <- data.frame(cbind(resp_matrix, dat_log_like))

Sigma_1 <- get_sigmad(K, t, s_1, l_1)
Sigma_2 <- get_sigmad(K, t, s_2, l_2)

xx <- (pi_1*apply(data_1, MARGIN = 2,  dmvnorm, mean = rep(0, N), sigma = Sigma_1, log=F))/
  (pi_2*apply(data_1, MARGIN = 2,  dmvnorm, mean = rep(0, N), sigma = Sigma_2, log=F) +
     pi_1*apply(data_1, MARGIN = 2,  dmvnorm, mean = rep(0, N), sigma = Sigma_1, log=F))

r_j1 <- xx
# cluster 2 responsibilities
r_j2 <- 1-xx
# calculate m_g, m
m_g <- 1:2
m_g[1] <- sum(r_j1)
m_g[2] <- sum(r_j2)

# calculate pi_g
pi_1 <- m_g[1]/sum(m_g)
pi_2 <- m_g[2]/sum(m_g)
# update parameters
# length parameter
l_1 <- sum(r_j1*ML_solutions$.)/m_g[1]
l_2 <- sum(r_j2*ML_solutions$.)/m_g[2]

log_like_mixture(data_1, N, K, t, l_1, l_2, pi_1, pi_2, s_1, s_2)

round(xx, digits=3)


## COMPLETE DATA FULL RUN:

data_1 <- tofino_seasoned
# set dimensions:
N <- nrow(data_1)
t <- 1:N
# optimize GP hyperparameters:
ML_solutions <- tofino_l_opt %>% as.data.frame()

# create one hot encoding matrix:
resp_matrix <- permutations(2,10,0:1,repeats.allowed=TRUE) 

# intialize vector:
complete_likelihood <- 1:nrow(resp_matrix)

# first compute L1 and L2
for (i in 1:nrow(resp_matrix) ) {
  # Set complete assignment:
  r_j1 <- resp_matrix[i, ]
  r_j2 <- 1-r_j1
  # estimate params:
  
  # calculate m_g, m
  m_g <- 1:2
  m_g[1] <- sum(r_j1)
  m_g[2] <- sum(r_j2)
  
  l_1 <- sum(r_j1*ML_solutions$.)/m_g[1] 
  l_2 <- sum(r_j2*ML_solutions$.)/m_g[2] 
  pi_1 <- mean(r_j1)
  pi_2 <- mean(r_j2)
  s_1 <- 1
  s_2 <- 1
# compute complete data likelihood:
  # two likelihoods:
  like_1 <- apply(data_1, MARGIN = 2, dmvnorm, mean = rep(0, N), 
                  sigma = get_sigmad(K, t, s_1, l_1)) 
  
  like_2 <- apply(data_1, MARGIN = 2, dmvnorm, mean = rep(0, N), 
                  sigma = get_sigmad(K, t, s_2, l_2)) 
  
  # compute complete data like:
  complete_likelihood[i] <- prod(like_1*pi_1*r_j1 + like_2*pi_2*r_j2)
  
}

tofino_complete <- data.frame(cbind(resp_matrix, complete_likelihood))


tofino_complete[which.max(tofino_complete$complete_likelihood[-1]),]
# Victoria






