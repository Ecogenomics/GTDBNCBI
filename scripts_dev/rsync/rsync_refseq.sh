#!/bin/bash
rsync  -rtlvL --chmod=Dug=rwx,Do=rx,Fug=rw,Fo=r  --exclude=*/all_assembly_versions --exclude=*/representative --exclude=README.txt ftp.ncbi.nlm.nih.gov::genomes/refseq/archaea/ /srv/db/ncbi/from_ftp/genomes/refseq/archaea
rsync  -rtlvL --chmod=Dug=rwx,Do=rx,Fug=rw,Fo=r  --exclude=*/all_assembly_versions --exclude=*/representative --exclude=README.txt ftp.ncbi.nlm.nih.gov::genomes/refseq/bacteria/ /srv/db/ncbi/from_ftp/genomes/refseq/bacteria
