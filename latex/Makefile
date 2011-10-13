DATE=`date +%F`
NAME=ab_paper_draft
TEXFILE=whole_paper.tex

.PHONY: all clean html pdf epub embed

pre:
	mkdir -p build
post:
	@echo POST

clean:
	rm -rf build

pdf: pre
	# pdflatex --jobname=$(NAME)_$(DATE) --output-directory=build $(TEXFILE)
	pdflatex --jobname=$(NAME)_$(DATE) --output-directory=build whole_paper.tex

default: pdf
