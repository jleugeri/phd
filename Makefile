SRC_FILES := $(wildcard fig/*.tex)
PDF_FIGURES := $(patsubst fig/%.tex,fig/%.pdf,$(SRC_FILES))
SVG_FIGURES := $(patsubst fig/%.tex,fig/%-1.svg,$(SRC_FILES))

all: pdfs svgs
pdfs: $(PDF_FIGURES)
svgs: $(SVG_FIGURES) pdfs

list:
	@echo $(PDF_FIGURES)

fig/%.pdf: fig/%.tex
	lualatex --output-directory _build --jobname $* --shell-escape "\def\inputname{$<} \input{fig_template.tex}"
	cp _build/$*.pdf fig/

fig/%-1.svg: fig/%.pdf
	pdf2svg $< fig/$*-%d.svg all

