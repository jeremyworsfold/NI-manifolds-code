FIGURES := ants ants-errors tk2-leading tk2-correction tk2-volume

.PHONY: all figures clean

all: figures

figures:
	@mkdir -p figures/out
	@cd figures && for figure in $(FIGURES); do \
		latexmk -silent -pdf -interaction=nonstopmode -outdir=out $$figure.tex && \
		cp out/$$figure.pdf $$figure.pdf || exit 1; \
	done

clean:
	rm -rf figures/build figures/out
