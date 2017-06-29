all: slides

slides: *.Rtex
	Rscript -e "library(knitr); library(methods); knit2pdf('master.Rtex')"
