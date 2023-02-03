August 2021
vadr-models-corona-1.3-1

https://github.com/ncbi/vadr

Instructions for SARS-CoV-2 annotation using VADR:
https://github.com/ncbi/vadr/wiki/Coronavirus-annotation

VADR documentation:
https://github.com/ncbi/vadr/blob/master/README.md

See RELEASE-NOTES.txt for details on changes between model versions. 

------------

Recommended command for SARS-CoV-2 annotation 
(used by GenBank as of writing with vadr 1.3):

v-annotate.pl \
--split --cpu 8 \
--glsearch \
-s -r --nomisc \ 
--mkey sarscov2 \
--lowsim5seq 6 --lowsim3seq 6 \
--alt_fail lowscore,insertnn,deletinn \
--mdir <PATH-TO-THIS-MODEL-DIR> \
<fastafile> <outputdir>

The '--split --cpu 8' options will run v-annotate.pl multi-threaded on
8 threads. To change to '<n>' threads use '--split --cpu <n>', but
make sure you have <n> * 2G total RAM available. 

To run single threaded remove the '--split --cpu 8' options.

------------

Example of mapping model coordinates in VADR output files to the
NC_045512 model:

perl $VADRSCRIPTSDIR/miniscripts/vadr-map-model-coords.pl \
<PATH-TO-VADR-.alt.list-FILE-OR-GENBANK-SUBMISSION-PORTAL-DETAILED-ERROR-REPORT-TSV-FILE>
<PATH-TO-THIS-MODEL-DIR>/sarscov2.mmap \
NC_045512

This will add an additional tab-delimited field at the end of each
line with releavnt model oordinates in the NC_045512 model.

------------

Example of outputting a map of all model positions between two models
(from model NC_045512-MW422255 to model NC_045512):

perl $VADRSCRIPTSDIR/miniscripts/vadr-map-two-complete-models \
<PATH-TO-THIS-MODEL-DIR>/sarscov2.mmap \
NC_045512-MW422255 \
NC_045512 

-------------

Contact eric.nawrocki@nih.gov for help.
