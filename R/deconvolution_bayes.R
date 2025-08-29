


#' Bayesian Inference with two possible precisions for all genes
#'
#' @param X training data(nsamples*ngenes)
#' @param Y test data (vector with ngenes)
#' @param tau_eta precision for normal distribution of etas (logit link of probabilities)
#' @param alpha beta distribution parameter for the proportion of marker genes
#' @param beta beta distribution parameter for the proportion of marker genes
#' @param kappa1 gamma distribution parameter for the precision of useful genes
#' @param lambda1 gamma distribution parameter for the precision of useful genes
#' @param kappa2 gamma distribution parameter for the precision of nuisance genes
#' @param lambda2 gamma distribution parameter for the precision of nuisance genes
#' @param n.chains number of chains for MCMC
#' @param n.iter number of iterations in each chain
#' @param n.burnin number of burnins
#' @param n.thin thinning for posterior sampling
#'
#' @returns posterior mean of each parameters and mixing check
#' @importFrom R2jags jags
#' @importFrom coda mcmc gelman.diag
#' @importFrom stats rnorm rbinom
#' @export
#'
deconv_bayes_2 <- function(X, Y, tau_eta=0.1, alpha=1, beta=5,
                           kappa1=10, lambda1=1, kappa2=1, lambda2=3,
                           n.chains=5, n.iter=5e4, n.burnin=2e4, n.thin=4){

    # specify JAGS
    model_string <- "
    model {
      for (i in 1:N) {
        y[i] ~ dnorm(mu[i], tau[i])
        mu[i] <- inprod(X[i,], gamma[1:P])
        tau[i] <- Z[i] * tau1 + (1-Z[i]) * tau2
        Z[i] ~ dbern(pi)
      }

      for (j in 1:P){
        expeta[j] <- exp(eta[j])
      }
      sum_expeta <- sum(expeta[])
      for (j in 1:P){
        gamma[j] <- expeta[j]/sum_expeta
      }


      # normal prior for all the betas
      for (j in 1:P){
        eta[j] ~ dnorm(0, tau_eta)
      }

      # Mixture prior for variance
      pi ~ dbeta(alpha, beta)
      tau1 ~ dgamma(kappa1, lambda1)
      tau2 ~ dgamma(kappa2, lambda2)

    }
    "

    N <- ncol(X)
    P <- nrow(X)

    # specify prior parameters
    data_jags <- list(
        N = N,
        P = P,
        X = t(X),
        y = Y,
        tau_eta = tau_eta,
        alpha= alpha, beta=beta,  # prior for pi
        kappa1 = kappa1, lambda1=lambda1, # prior for tau1
        kappa2 = kappa2, lambda2=lambda2 # prior for tau2
    )

    # initialize parameters
    inits <- function() {
        list(eta = rnorm(n=P, mean=0, sd=1/sqrt(tau_eta)),
             Z = rbinom(n=N, size=1, prob=0.5),
             tau1=kappa1/lambda1, tau2=kappa2/lambda2, pi=alpha/(alpha+beta))
    }


    params <- c("eta", "Z", "tau1", "tau2", "pi")

    # run MCMC
    fit <- jags(
        data = data_jags,
        parameters.to.save = params,
        model.file = textConnection(model_string),
        n.chains = n.chains, n.iter = n.iter,
        n.burnin = n.burnin, n.thin = n.thin,
        inits = inits
    )

    # check mixing
    sim_array <- fit$BUGSoutput$sims.array
    deviance_chains <- sim_array[, , "deviance"]
    deviance_chain_list <- split(deviance_chains, col(deviance_chains))
    deviance_chains_mcmc_list <- lapply(deviance_chain_list,
                                        function(onechain) mcmc(onechain))
    gelman_stats <- gelman.diag(deviance_chains_mcmc_list)


    eta_mean <- fit$BUGSoutput$mean$eta
    gamma_mean <- exp(eta_mean) / sum(exp(eta_mean))
    Z_mean <- fit$BUGSoutput$mean$Z
    tau1_mean <- fit$BUGSoutput$mean$tau1
    tau2_mean <- fit$BUGSoutput$mean$tau2
    pi_mean <- fit$BUGSoutput$mean$pi

    output <- list(gamma=gamma_mean, Z=Z_mean, tau1=tau1_mean, tau2=tau2_mean,
                   pi=pi_mean, gelman=gelman_stats)

    return(output)

}



