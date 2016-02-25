#!/bin/bash
#Permission is always an issue during the rsync because it is rarely the same user running
# we use -rtl instead of -a to avoid resetting owner,group and permission
# we use the --chmod flag to allow all dataadmin to edit the records
# -r recurse into directories
#-t preserve modification times
#-l copy symlinks as symlinks
#-v increase verbosity
#-L transform symlink into referent file/dir
#--chmod=Dug=rwx,Do=rx,Fug=rw,Fo=r give 775 for directories and 664 for files 
rsync -rtlvL --chmod=Dug=rwx,Do=rx,Fug=rw,Fo=r --exclude=*/all_assembly_versions --exclude=*/representative --exclude=README.txt ftp.ncbi.nlm.nih.gov::genomes/genbank/archaea/ /srv/db/ncbi/from_ftp/genomes/genbank/archaea
rsync -rtlvL --chmod=Dug=rwx,Do=rx,Fug=rw,Fo=r --exclude=*/all_assembly_versions --exclude=*/representative --exclude=README.txt ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/ /srv/db/ncbi/from_ftp/genomes/genbank/bacteria
