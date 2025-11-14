# library(devtools)
# load_all()
# document()

library(ProbBreed)
# source("Data/bpsi.R")
getwd()
met_df=read.csv("https://raw.githubusercontent.com/tiagobchagas/BPSI/refs/heads/main/Data/blues_long.csv",header=T)

head(met_df)
mod = bayes_met(data = met_df,
                gen = "gen",
                loc = "env",
                repl = NULL,
                trait = "PH",
                reg = NULL,
                year = NULL,
                res.het = T,
                iter =400, cores = 4, chain = 4) #recommended run at least 4k iterations


mod2 = bayes_met(data = met_df,
                 gen = "gen",
                 loc = "env",
                 repl = NULL,
                 trait = "GY",
                 reg = NULL,
                 year = NULL,
                 res.het = T,
                 iter = 400, cores = 4, chain = 4) #recommended run at least 4k iterations

mod3 = bayes_met(data = met_df,
                 gen = "gen",
                 loc = "env",
                 repl =  NULL,
                 trait = "NDM",
                 reg = NULL,
                 year = NULL,
                 res.het = T,
                 iter = 400, cores = 4, chain = 4) #recommended run at least 4k iterations



models=list(mod,mod2,mod3)
names(models) <- c("PH","GY","NDM")
inc=c(FALSE,TRUE,FALSE)
names(inc) <- names(models)

results <- lapply(names(models), function(model_name) {
  x <- models[[model_name]]  # actual model object

  outs <- extr_outs(model = x,
                    probs = c(0.05, 0.95),
                    verbose = TRUE)

  a <- prob_sup(extr = outs,
                int = .2,
                increase = inc[[model_name]],  # ← now model_name is a character!
                save.df = FALSE,
                verbose = TRUE)

  return(a)
})
names(results) <- names(models)




problist=results;
increase = c(FALSE,TRUE,FALSE);
int = 0.1;
lambda =c(1,2,1);
save.df = F

bpsi=BPSI(problist=results,
          increase = c(FALSE,TRUE,FALSE),
          int = 0.1,
          lambda =c(1,2,1),
          save.df = F)


# BPSI_result=bpsi
# head(bpsi)
plot(bpsi,category = "Ranks")

plot(bpsi,category = "BPSI")

df=print(bpsi)
