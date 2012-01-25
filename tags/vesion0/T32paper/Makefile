#
draft.pdf: draft.dvi
	dvipdfmx draft

draft.dvi: draft.tex tex/beamenergy.tex tex/BeamQuality.tex tex/CrossTalk.tex tex/Data-MCComparison.tex tex/DataQuality.tex tex/EventReconstruction.tex tex/EventSimulation.tex tex/hitclustering.tex tex/Introduction.tex tex/LiquidArgonPurity.tex tex/NoiseSim.tex tex/RecombinationFactor.tex tex/SoftwareFramework.tex tex/StoppedKaon2.tex tex/StoppedKaon.tex tex/StoppedProton2.tex tex/Summary.tex
	latex draft;latex draft;latex draft

clean: 
	-rm draft.dvi draft.pdf

distclean:
	-rm draft.aux tex/*.aux  