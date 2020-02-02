devtools::document()
devtools::check(vignettes=FALSE)
devtools::build()

R CMD REMOVE fusionImage
R CMD INSTALL  /media/disk/Proyectos_2020/GITS/fusionImage/linux/fusionImage_0.0.1.tar.gz 

detach(package:fusionImage, unload=TRUE)
library(fusionImage)

devtools::build_vignettes()
devtools::check()
devtools::build()


git add .
git status
git commit -m "Definitivo"
git push

