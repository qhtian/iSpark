iSpark algorithm

The identification of spatiotemporally restricted Ca2+ signals, Ca2+ sparks, was instrumental for our understanding of cardiac Ca2+ homeostasis. High-speed 2D confocal imaging enables acquisition of such Ca2+ sparks with high-content information but their full appreciation is constrained by the lack of unbiased and easy-to-use analysis tools. We developed a software toolset for unbiased and automatic Ca2+ spark analysis for huge data sets of subcellular Ca2+ signals. iSpark was developed to be scanner and detector independent. In myocytes from hearts subjected to various degrees of hypertrophy we acquired more than 5.000.000 Ca2+ sparks from 14 mice. The iSpark-enabled analysis of this large Ca2+ spark data set showed that the highly organized distribution of Ca2+ sparks present in healthy cells disarrayed concomitant with the development of aberrant transverse tubules and disease severity. Thus, iSpark represents a versatile and universal tool for analyzing local Ca2+ signaling in healthy as well as diseased, aberrant local Ca2+ signal transduction. The results from the unbiased analysis of large data sets provide a deeper insight into possible mechanisms contributing to the onset and progression of cardiac diseases such as hypertrophy.

Please refer to the formal paper published in Journal of Molecular and Cellular Cardiology (in press):
Large scale, unbiased analysis of elementary calcium signaling events in cardiac myocytes
by Qinghai Tian (1), Laura Schröder (1), Yvonne Schwarz (2), Aline Flockerzi (1), Lars Kaestner (1,4), Andre Zeug (3), Dieter Bruns (2), and Peter Lipp (1)
1 Center for Molecular Signaling (PZMS), Institute for Molecular Cell Biology, Research Center for Molecular Imaging and Screening; Medical Faculty, Saarland University, Homburg/Saar, Germany; 
2 Department of Physiology, Medical Faculty, Saarland University, Homburg, Germany; 
3 Cellular Neurophysiology, Hannover Medical School, Hannover, Germany.
4 Experimental Physics, Faculty NT, Saarland University, Saarbrücken, Germany


Notes to use the codes:
  1. Please install a copy of ImageJ, and locate the pacage "ij.jar" or similar "ij-xxx.jar";
  2. Please get a copy of PureDenoise_.jar from http://bigwww.epfl.ch/algorithms/denoise/.
  3. Please get a copy of bioformats_package.jar from https://www.openmicroscopy.org/bio-formats/downloads/.
  4. Open MatLab, run the command:
    >> cd(prefdir)
    >> edit javaclasspath.txt
  5. Write down the path to the pacages listed above in the editor like (e.g. in my Mac OS):
      /Applications/MATLAB_R2014a.app/toolbox/imagejlib/ij.jar
      /Applications/MATLAB_R2014a.app/toolbox/imagejlib/PureDenoise_.jar
      /Applications/MATLAB_R2014a.app/toolbox/imagejlib/bioformats_package.jar
  6. In MatLab's preferences, navigate to MATLAB / General / Java Heap Memory, adjust it to a maximal value if possible.
  7. Restart MatLab and you are ready to go.

  If you do not have a MatLab license and thus need a copy of compiled version, please contact me directly.


Qinghai Tian (tian_qhcn@icloud.com)

	Center for Molecular Signaling (PZMS)
	Institute for Molecular Cell Biology
	Research Center for Molecular Imaging and Screening
	Building 61, Medical Faculty
	Saarland University Hospital
	Saarland University
	66421 Homburg/Saar, Germany
	www.pzms.uni-saarland.de   -or-   www.lipplab.de