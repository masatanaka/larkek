#
draft.pdf: draft.dvi
	dvipdfmx draft

draft.dvi: draft.tex tex/Data.tex tex/Introduction.tex tex/Kaon.tex tex/Pion.tex tex/Proton.tex tex/Recombination.tex tex/Setup.tex tex/Simulation.tex tex/Summary.tex
	latex draft;latex draft;latex draft

clean: 
	-rm draft.dvi draft.pdf

distclean:
	-rm draft.aux tex/*.aux  