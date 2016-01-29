#!/bin/bash
rsync  -avL  --exclude=*/all_assembly_versions --exclude=*/representative --exclude=README.txt ftp.ncbi.nlm.nih.gov::genomes/genbank/archaea/ /srv/db/ncbi/from_ftp/genomes/genbank/archaea
rsync  -avL  --exclude=*/all_assembly_versions --exclude=*/representative --exclude=README.txt ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/ /srv/db/ncbi/from_ftp/genomes/genbank/bacteria
