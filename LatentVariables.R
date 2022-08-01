
#  ======================  #
#  Latent variable models  #
#  ======================  #
library(data.table)
library(dplyr)
library(bit64)
library(CVXR)
library(tidyr)
library(igraph)
library(parallel)

autPair <- fread('autPair.txt') %>% 
  .[Year <= 2019]

#  bootstrap estimation of latent variable models
bootstrap <- function(autPair1){
  
  #  sample papers with replacement
  autPair1 <- autPair1 %>% 
    .[, idx := c(1:.N)] #  row index to match sampled papers
  rowIdx <- data.table(idx = sample(c(1:nrow(autPair1)), 
                                    nrow(autPair1), replace = TRUE))
  autPair1 <- rowIdx %>%  #  match the paper index
    .[autPair1, on = .(idx), nomatch = 0, allow.cartesian = TRUE] %>% 
    .[, idx := NULL]
  
  autPair1 <- autPair1[, `:=`(maxYr = max(Year), minYr = min(Year)), 
                       .(AuthorId, i.AuthorId)] %>% 
    .[, duration := maxYr - minYr + 2]
  
  A <- data.table(i0 = autPair1$AuthorId,
                  j0 = autPair1$i.AuthorId,
                  elitePap = autPair1$elitePap,
                  duration = autPair1$duration) %>% 
    .[, .(pub = .N, elite = sum(elitePap)), .(duration, i0, j0)]
  
  #  reorder author id; starting from 1
  AuthorId <- unique(c(autPair1$AuthorId, autPair1$i.AuthorId))
  id0 <- as.integer(as.factor(AuthorId))
  id <- data.table(id = id0, AuthorId = AuthorId)
  
  #  build matrix ready for model
  A1 <- A[id, on = c('i0==AuthorId'), nomatch = 0] %>% 
    .[id, on = c('j0==AuthorId'), nomatch = 0] %>% 
    setnames(c('t', 'i0', 'j0', 'pub', 'elite', 'i', 'j')) %>% 
    .[, .(i, j, pub, elite, t, i0, j0)]
  
  ###  remove tree-like structures  ###
  m <- as.matrix(A1[, .(i, j)])
  g <- graph.edgelist(m, directed = FALSE)
  membership0 <- components(g)
  memb.tmp <- data.table(membership = 1:length(membership0$csize),
                         csize = membership0$csize) #  component size
  membership <- data.table(id = c(1:length(membership0$membership)),
                           membership = membership0$membership) %>% 
    .[memb.tmp, on = 'membership']
  
  #  components, number of students + number of advisors
  #  should be <= number of links; for valid MCMC estimation
  A2 <- A1 %>% 
    .[membership, on = 'i==id', nomatch = 0] %>% 
    .[, `:=`(N_i = n_distinct(i), N_j = n_distinct(j),
             N_links = .N), .(membership)] %>% 
    .[csize <= N_links & N_links > 2] %>% 
    .[, .(i0, j0, pub, elite, t)]
  
  if(nrow(A2) == 0){return(NULL)}  # return NULL DT
  
  #  reorder author id AGAIN; starting from 1
  lvl.1 <- unique(c(A2$i0, A2$j0))
  id0.1 <- as.integer(as.factor(lvl.1))
  Id <- data.table(id.1 = id0.1, AuthorId = lvl.1)
  
  A3_0 <- A2[Id, on = c('i0==AuthorId'), nomatch = 0] %>% 
    .[Id, on = c('j0==AuthorId'), nomatch = 0] %>% 
    setnames(c('i0', 'j0', 'pub', 'elite', 't', 'i1', 'j1'))
  return(A3_0)
}

#  AuthorId for all data
ID <- data.table(AuthorId = unique(autPair$AuthorId)) %>% 
  .[order(AuthorId)]

binomModel <- function(MODEL, YEAR){
  #  YEAR: include papers published before it
  
  autPair2 <- autPair[Year <= YEAR]
  
  #  bootstrap publications
  A3 <- bootstrap(autPair2)
  if(is.null(A3)){return(NULL)} # return NULL DT
  #  AuthorId for all data
  ID <- data.table(AuthorId = unique(autPair2$AuthorId)) %>% 
    .[order(AuthorId)]
  
  #  all author ids
  A3.id <- data.table(AuthorId = c(A3$i0, A3$j0), 
                      id = c(A3$i1, A3$j1)) %>% 
    .[, .I[1], .(AuthorId, id)] %>% 
    .[, V1 := NULL] %>% 
    .[order(id)]
  
  ###  Productivity model  ###
  if(MODEL == 'productivity'){
    m2 <- A3[, c('i1', 'j1')] %>% as.matrix
    pub <- Constant(A3$pub)
    t <- Constant(A3$t)
    
    N <- c(m2[, 1:2]) %>% n_distinct
    lambda <- Variable(N)
    obj <- sum(-(lambda[m2[, 1]] + lambda[m2[, 2]]) * t + 
                 pub * log((lambda[m2[, 1]] + lambda[m2[, 2]]) * t))
    constraints <- list(lambda >= 0.01)
    
    problem <- Problem(Maximize(obj), constraints)
    result <- psolve(problem, solver = "SCS", max_iters = 500) # installed_solvers()
    print(result$value)
    lambda_res <- result$getValue(lambda)
  }
  
  
  ###  Prominence model  ###
  if(MODEL == 'prominence'){
    elite <- A3$elite %>% as.integer
    pub <- A3$pub %>% as.integer
    student <- A3$i1 %>% as.integer
    supervisor <- A3$j1 %>% as.integer
    
    N <- max(c(student, supervisor)) %>% as.integer
    V <- length(student) %>% as.integer
    lambda <- Variable(N)
    
    obj <- sum(elite * log(lambda[student] + lambda[supervisor]) +
                 (pub - elite) * log(1 - (lambda[student] + lambda[supervisor])))
    constraints <- list(lambda >= 1e-6, lambda <= 1 - 1e-4,
                        lambda[student] + lambda[supervisor] >= 1e-6, 
                        lambda[student] + lambda[supervisor] <= 1 - 1e-4)
    
    problem <- Problem(Maximize(obj), constraints)
    result <- psolve(problem, solver = "SCS", max_iters = 1000) 
    print(result$value)
    lambda_res <- result$getValue(lambda)
  }
  
  Res <- data.table(id = c(1:length(lambda_res)), 
                    lambda = lambda_res)
  
  A3.id <- A3.id %>% 
    .[, lambda := lambda_res]
  A3.id <- A3.id[ID, on = .(AuthorId)] %>% 
    .[, id := NULL] # AuthorId, theta
  
  return(A3.id)
}


#  choose which model to run: productivity or prominence
model <- 'productivity'

REP <- 100  #  number of replications
for(yr in c(1975:2017)){
  print(yr)
  
  #  prepare data for the cluster
  t1 <- Sys.time()
  ncores <- 20
  cl <- makeCluster(mc <- getOption("cl.cores", ncores))
  clusterEvalQ(cl, {library(CVXR); library(dplyr); 
    library(data.table); library(igraph)})
  clusterExport(
    cl = cl, 
    varlist = c("autPair", 'binomModel', 'bootstrap', 'model', 'yr')
  )
  
  RES0 <- parLapply(cl, c(1:REP), function(x) binomModel(model, yr))
  stopCluster(cl)
  
  RES <- ID
  
  if(is.null(RES0[[1]])){next}
  for(i in 1:length(RES0)){
    RES <- RES0[[i]][RES, on = .(AuthorId)]
  }
  
  fwrite(RES, paste0('lambda_', model, '_', yr, '.txt'))
  print(Sys.time() - t1)
}












#  ==============================================  #
#  Aggregate all fields for lambda-theta in years  #
#  ==============================================  #
library(data.table)
library(dplyr)
library(bit64)
library(matrixStats)
library(igraph)

lambda1 <- NULL
autPair1 <- NULL

#  function to replace NA with 0
f_dowle2 = function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)), (i):=0]
}

FLD <- 'mathematics'
print(FLD)
autPair0 <- fread('autPair.txt')
Res <- NULL
for(yr in c(1975:2017)){
  lambda0 <- fread(paste0('lambda_poisson_', yr, '.txt'))
  theta0 <- fread(paste0('lambda_binomial_', yr, '.txt'))
  
  #  lambda mean, SD for individual authors
  lambda <- lambda0 %>% 
    .[, `:=`(n_1 = rowCounts(as.matrix(.SD), value = NA),
             nAll = ncol(.) - 1), 
      .SDcols = c(2:ncol(.))] %>% 
    .[, nLam := nAll - n_1] %>% 
    .[nLam >= 1]  #  at least one non-NA value
  f_dowle2(lambda)

  lambda <- lambda %>% 
    .[, `:=`(lambda = rowMeans(.SD, na.rm = TRUE)#, lambdaSD = rowSds(as.matrix(.SD), na.rm = TRUE)
    ), 
    .SDcols = c(2:(ncol(.) - 3))] %>% 
    .[, lambda := ifelse(lambda < 1e-2, 1e-2, lambda)] %>% 
    .[, .(AuthorId, lambda)] 
  
  #  theta mean, SD for individual authors
  theta <- theta0 %>% 
    .[, `:=`(n_1 = rowCounts(as.matrix(.SD), value = NA),
             nAll = ncol(.) - 1), 
      .SDcols = c(2:ncol(.))] %>% 
    .[, nLam := nAll - n_1] %>% 
    .[nLam >= 1]  #  at least one non-NA value
  f_dowle2(theta)

  theta <- theta %>% 
    .[, `:=`(theta = rowMeans(.SD, na.rm = TRUE) #, thetaSD = rowSds(as.matrix(.SD), na.rm = TRUE)
    ), 
    .SDcols = c(2:(ncol(.) - 3))] %>% 
    .[, theta := ifelse(theta < 1e-6, 1e-6, theta)] %>%  # min(theta) = 1e-6
    .[, .(AuthorId, theta)] 
  
  #  combine lambda & theta
  Res_tmp <- lambda %>% 
    .[theta, on = .(AuthorId), nomatch = 0] %>% 
    .[, Year := yr]
  Res <- rbind(Res_tmp, Res)
}

fwrite(Res, paste0(FLD, '/autLamTheta_yr.txt'))







