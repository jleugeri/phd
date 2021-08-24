SRC_FILES := $(wildcard fig/*.tex)
FIGURES := $(patsubst fig/%.tex,fig/%.pdf,$(SRC_FILES))


all: $(FIGURES)

list:
	@echo $(FIGURES)

fig/%.pdf: fig/%.tex
	lualatex --output-directory _build --jobname $* --shell-escape "\def\inputname{$<} \input{fig_template.tex}"
	cp _build/$*.pdf fig/

