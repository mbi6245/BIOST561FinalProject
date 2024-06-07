#' Title
#'
#' @param obj a fitted glmer model for which to calculate the robust variance-covariance matrix
#' @param cluster optional expression/vector indicating which observations belong to the same cluster
#' @param type character string specifuing which robust variance form should be used with available options "classic", "DF", "KC", "MD"
#'
#' @return robust variance/covariance matrix for fitted glmer model
#' @export vcovCR_glmerMod_Max
#'
#' @examples
vcovCR_glmerMod_Max = function(obj, cluster = NULL, type="classic") {
  beta=matrix(fixef(obj),ncol=1)
  np=dim(beta)[1]
  ranef_names = names(ranef(obj))
  gamma = matrix(lapply(ranef_names, function(name) ranef(obj)[[name]]))
  # gamma = matrix(c(ranef(obj)$clustime,ranef(obj)$fcluster),ncol=1)
  nq = dim(gamma)[1]
  X = model.matrix(obj,type="fixed")
  Z = model.matrix(obj,type="random")
  Y = obj@resp$y # example of "slots" in R
  if(is.null(cluster)) {
    cluster = clubSandwich:::get_outer_group(obj) # maybe there is a better way to do this
  }
  eta = predict(obj,type="link")
  ginv_eta = predict(obj,type="response")
  link = family(obj)$link
  switch(link,
         "identity" = {
           delta = diag(nobs(obj))
           deltainv = delta
         },
         "logit" = {
           delta = diag(exp(eta) / (1 + exp(eta))^2)
           deltainv = diag((1 + exp(eta))^2 / exp(eta))
         },
         "log" = {
           delta = diag(exp(eta))
           deltainv = diag(1 / exp(eta))
         }
  )

  if (type == "DF") {
    c = dim(X)[1] / (dim(X)[1] - dim(X)[2])
  }

  theta = as.data.frame(VarCorr(obj))
  G = length(levels(cluster))
  # G = ngrps(obj)["fcluster"]
  # sigma = theta$vcov[theta$grp=="Residual"] * diag(nobs(obj))

  # `family(obj)$variance returns function of variance wrt mu`
  sigma2 = sigma(obj)^2
  sigma_mat = family(obj)$variance(ginv_eta) * sigma2
  lambda = getME(obj,"Lambda")
  R = lambda%*%t(lambda)*sigma2
  # R = numeric(0)
  # for (i in 1:(dim(theta)[1] - 1)) {
  #  R = c(R, rep(theta$vcov[i],
  #               ngrps(obj)[ranef_names[i]]))
  # }
  # R = diag(R)
  # R = diag(c(rep(theta$vcov[1], ngrps(obj)["clustime"]),rep(theta$vcov[2],ngrps(obj)["fcluster"])))
  mbv = vcov(obj)

  # create gamma vector to calculate P and consequently e
  # gamma_vec = numeric(0)
  # first component of original gamma_vec is gamma[1, 1][[1]][, 1]
  # for(i in 1:nq) {
  #   gamma_temp = gamma[i, 1][[1]]
  #   for (j in 1:dim(gamma_temp)[2]) {
  #     gamma_vec = c(gamma_vec, gamma_temp[, j])
  #   }
  # }
  # gamma_vec = c(gamma[1, 1][[1]], gamma[2, 1][[1]])

  P = deltainv %*% (Y - ginv_eta) + eta # X %*% beta + Z %*% gamma_vec

  e = P - X %*% beta

  # sum product of matrices over over clusters
  running_sum = matrix(0, nrow = dim(X)[2], ncol = dim(X)[2])
  WB_C = solve(R)
  for (clus in unique(cluster)) {
    clus_idx = which(cluster == clus)

    WB_A = diag(1 / diag(deltainv[clus_idx, clus_idx] %*% diag(sigma_mat[clus_idx]) %*% deltainv[clus_idx, clus_idx]))
    WB_U = Z[clus_idx, ]
    WB_V = t(Z[clus_idx, ])
    W = WoodburyMatrix(A = WB_A, B = WB_C, U = WB_U, V = WB_V)

    Vinv_g = solve(W)
    # mbv_g = mbv[clus_idx, ]
    X_g = X[clus_idx, ]
    e_g = e[clus_idx, ]

    switch(type,
           "KC" = {
             H_g = X_g %*% mbv %*% t(X_g) %*% Vinv_g
             pre_F_g = diag(dim(H_g)[1]) - t(H_g)
             F_g = sqrtm(pre_F_g)
             running_sum = running_sum + t(X_g) %*% Vinv_g %*% t(F_g) %*% e_g %*% t(e_g) %*% F_g %*% Vinv_g %*% X_g
           },
           "MD" = {
             H_g = X_g %*% mbv %*% t(X_g) %*% Vinv_g
             pre_F_g = diag(dim(H_g)[1]) - t(H_g)
             F_g = solve(pre_F_g)
             running_sum = running_sum + t(X_g) %*% Vinv_g %*% t(F_g) %*% e_g %*% t(e_g) %*% F_g %*% Vinv_g %*% X_g
           },
           {
             running_sum = running_sum + t(X_g) %*% Vinv_g %*% e_g %*% t(e_g) %*% Vinv_g %*% X_g
           }
    )
  }
  c = 1

  # final robust variance formula
  robustVar = c * mbv %*% running_sum %*% mbv
  return(diag(robustVar))
}
