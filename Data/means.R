
library(ProbBreed)
library(tidyverse)
library(asreml)
library(openxlsx)
# pheno <- read.csv("https://raw.githubusercontent.com/mauricioaraujj/Pan_African_Trials_Network/refs/heads/main/data/Malawi_data.csv",sep=";")
pheno <- read.csv("Data/pat_data.csv",sep=";")
# pheno <- pheno[,c(3,4,5,7,14,15,20,11,33,34)]
head(pheno)
colnames(pheno)
traits <- c(colnames(pheno)[c(4,5,6)])
head(pheno)

# hist(pheno$LOD)
# colnames(pheno)
# str(traits)
factors <- c(colnames(pheno)[c(1,2,3)])
factors
pheno[, factors] <- lapply(pheno[, factors], as.factor)
str(pheno)

# data.list = split(pheno, f = pheno$env)

pheno[, traits] <- lapply(pheno[, traits], as.numeric)
str(pheno)
stmt = list()
j="E01"
i="PH"
unique(pheno$env)
unique(pheno$gen)
str(pheno)
data.list <- split(pheno,f=pheno$env)
Blue=list()


for (j in names(data.list)) {

  summar.mat = matrix(nrow = length(traits), ncol = 4,
                      dimnames = list(traits,
                                      c("geno","error","H2", "CV")))

  vc.env = list()
  Blue.env = list()
  head(pheno)


  for (i in traits) {

  message("====> Trait:", i)




  x = data.list[[j]]
  message("====> Environment:", j)




    if(all(is.na(x[,i]))) next

  str.mod <- as.character(NA)
  str.mod[1] <- paste0("asreml::asreml(fixed=", i,  "~ rep,
                                                          random = ~ gen,maxit = 50,
                                                          data = x,na.action = na.method(x = 'include', y = 'include'))")
  str.mod[2] <- paste0("asreml::asreml(fixed=", i,  "~ rep+gen,
                     maxit = 50,
                     residual = ~units,
                     data = x, na.action = na.method(x = 'include', y = 'include'))")


  # Combine into a complete model formula
  code <- paste(str.mod[1])

  # Fit the model for gmoda
  model.r <- tryCatch(expr = {
    eval(parse(text = code))
  }, error = function(msg) {
    message("Error in gmoda fitting: ", msg)
    return(NULL)
  })

  if(is.null(model.r)) next

  code <- paste(str.mod[2])

  # Fit the model for gmoda
  model.f <- tryCatch(expr = {
    eval(parse(text = code))
  }, error = function(msg) {
    message("Error in gmoda fitting: ", msg)
    return(NULL)
  })


    # colocar nugget and spline
    # model.f = asreml(fixed = x[,i] ~ rep+gen,
    #                  maxit = 50,
    #                  residual = ~units,
    #                  data = x, na.action = na.method(x = "include", y = "include"))


    # model.r = asreml(fixed = x[,i] ~ rep,
    #                  random = ~ gen,maxit = 50,
    #                  data = x)

    if(!model.r$converge) next

    if(any(na.exclude(model.r$vparameters.pc > 1))) model.r = update.asreml(model.r)
    # varcomp = summary(modelhvan)$varcomp
    varcomp = summary(model.r)$varcomp
    varcomp.df = data.frame(
      effect = c('genotypic', "residual"),
      component = c(round(
        varcomp[grep('geno', rownames(varcomp)),1],4),
        round(varcomp[grep('!R', rownames(varcomp)),1],4))
    )
    # if(!modelad$converge) next
    #
    # if(any(na.exclude(modelad$vparameters.pc > 1))) modelad = up.mod(modelad)
    #
    if(!model.f$converge) next

    if(any(na.exclude(model.f$vparameters.pc > 1))) model = update.asreml(model.f)

    lrtest = as.data.frame(lrt(asreml(fixed = x[,i] ~ rep,
                                      residual = ~units, maxit = 50,
                                      data = x,na.action = na.method(x = 'exclude', y = 'exclude')),model.r))
    # lrtest=update.asreml(lrtest)

    lrtest$effect = c("genotypic")

    varcomp= full_join(varcomp.df, lrtest[,-1], by = "effect")
    varcomp$signi = ifelse(varcomp$`Pr(Chisq)` <= 0.06, "*", "ns")
    varcomp$signi2 = ifelse(varcomp$`Pr(Chisq)` <= 0.06,
                            paste(round(varcomp$component, 4), "*"),
                            paste(round(varcomp$component, 4), "ns"))
    rownames(varcomp)=varcomp$effect

    pred = predict(model.r, classify = "gen", sed = TRUE)
    pred$sed= pred$sed[-which(pred$pvals$status=="Aliased"),-which(pred$pvals$status=="Aliased")]
    MVdelta = mean((pred$sed^2)[upper.tri(pred$sed^2,diag = F)])
    H2 = 1-(MVdelta/(2*varcomp[varcomp$effect == "genotypic","component"]))
    h2 = vpredict(model.r, h2~(V1)/(V1+V2))



    #h2ad = vpredict(modelad, h2~(V3)/(V3+V4))

    CV = sqrt(varcomp[grep("r", rownames(varcomp)),"component"])/
      mean(x[,i], na.rm = TRUE)*100

    summar.mat[i,] = c(varcomp$signi2[1], round(varcomp$component[2],4),
                       round(h2$Estimate,4), round(CV,4))



    vc.env[[i]] = varcomp.df

    Blue[[i]][[j]] = pred$pvals


  stmt[[j]] = list(summary = summar.mat,vc=vc.env)

  }

  rm(summar.mat, varcomp ,
     h2, CV, x)

}



  diagnostic = do.call(rbind, lapply(stmt, function(x){
  aa = as.data.frame(x$summary) |> rownames_to_column("trait")
  aa
})) |> rownames_to_column("env") |>
  mutate_at("env", str_replace, "\\..*", "") |>
  dplyr::select(env, trait, geno) |>
  separate(geno, into = c("genovar", "sig"), sep = " ")
print(diagnostic)
i="LOD"

pheno_filt <- droplevels(pheno[
  which(pheno$env %in% diagnostic[which(diagnostic$trait %in% i &
                                          diagnostic$sig == "*"), "env"]), ])


a=do.call(rbind,lapply(lapply(Blue, function(x) {
  do.call(rbind,(x))}),function(x){
   aa = as.data.frame(x) |> rownames_to_column("env")
  aa  }))
a

Blues=a |> rownames_to_column("trait") |> select(-status,-std.error)
Blues$trait=sub("\\..*","",Blues$trait)
Blues$env=sub("\\..*","",Blues$env)
head(Blues)
Blues_filt <- droplevels(Blues[
  which(c(Blues$trait,Blues$env)  %in% diagnostic[which(diagnostic$trait %in% i &
                                          diagnostic$sig == "*"), c("trait","env")]), ])

Blues_filt <- droplevels(Blues[
  which(paste(Blues$trait, Blues$env) %in%
          paste(diagnostic$trait[diagnostic$trait %in% traits & diagnostic$sig == "*"],
                diagnostic$env[diagnostic$trait %in% traits & diagnostic$sig == "*"])),
])
head(Blues_filt)

traits
Blues_long=Blues_filt |> pivot_wider(names_from = c("trait"),values_from = "predicted.value")
Blues_long


hist(Blues_long$PH)
tail(Blues_long)
write.csv(Blues_long,file="Data/blues_long.csv")
