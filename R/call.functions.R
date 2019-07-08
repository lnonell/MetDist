#calls necessary functions before package creation
codeDir <- "D:/Doctorat/Simplex/MetDist/R"

source(file=file.path(codeDir,"best.dist.R"))
source(file=file.path(codeDir,"best.aic.ks.R"))
source(file=file.path(codeDir,"best.aic.ks.with.betabin.R"))

source(file=file.path(codeDir,"est.all.params.R"))
source(file=file.path(codeDir,"est.betabin.params.R"))

source(file=file.path(codeDir,"fn.models.parallel.R"))
#source(file=file.path(codeDir,"fn.models.parallel.onlyZOIP.R")) #només per a completar els models. dESPRES BORRAR
source(file=file.path(codeDir,"fn.models.betabin.parallel.R"))
source(file=file.path(codeDir,"fn.simu.dml.R"))
source(file=file.path(codeDir,"fn.simu.non.dml.R"))
source(file=file.path(codeDir,"fn.simulations.R"))
source(file=file.path(codeDir,"apply.limma.R"))
source(file=file.path(codeDir,"fn.simulations.betabin.R"))

source(file=file.path(codeDir,"models.eval.R"))
source(file=file.path(codeDir,"res.chr.R"))
source(file=file.path(codeDir,"find.dml.n.R"))
