# KR_Python
*Examination work for Python course*

**Introduction into name-giver tool**

I am glad to present you the name-giver tool for contigs written in .fasta format.

This simple script makes requests to NCBI blastn service, and names your contigs by given alignment information.

Name for contig is choosing by expect value: if in 5 main alignments expect values are splited too much - tool will choose name from alignment with minimal except value. Otherwise, if some except values are close enough - will be choosen name, which represented oftenly then others.

---Argparse values---

Arguments for launching:
  * -i --input  (path to input file)
  * -fr --fresh (fresh start. Tool automatically saves progress to avoid accidents. By using this parameter you can launch tool without dump of last session)
  
  
---Output---

As output you will get an new .fasta file, named such as your input file, but with suffix "_ out"
Contigs in output file will be grouped by given names.
