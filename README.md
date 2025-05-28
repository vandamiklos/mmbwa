# mmbwa

Developed to improve alignment accuracy of highly fragmented variants while minimizing runtime 
using minimap2 and bwa-mem. 
Reads are first aligned with minimap2, the output is piped into a filter script that selects 
reads with long soft-clips and directs them to bwa-mem. 
The minimap2 and bwa-mem alignments are merged at the end.