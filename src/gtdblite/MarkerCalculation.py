import tempfile
import os
import shutil

from gtdblite.Exceptions import GenomeDatabaseError

from checkm_static import prodigal

def RunProdigalOnGenomeFasta(fasta_path):
    try:
        prodigal_dir = tempfile.mkdtemp()
    
        runner = prodigal.ProdigalRunner(prodigal_dir)
        if runner.run(fasta_path) is None:
            raise GenomeDatabaseError("Failed to run prodigal.")
        
        return prodigal_dir
    
    except:
        if prodigal_dir:
            shutil.rmtree(prodigal_dir)
        raise

    
def CalculateBestMarkerOnProdigalDir(marker_hmm_file, prodigal_dir):
    
    os.system("hmmalign --outformat Pfam -o %s %s %s" % (
        os.path.join(prodigal_dir, "%i.aligned" % (marker_id,)),
        marker_path,
        os.path.join(prodigal_dir, 'genes.faa')
    ))
    
    fh = open(os.path.join(prodigal_dir, "%i.aligned" % (marker_id,)))
    fh.readline()
    fh.readline()
    seqline = fh.readline()
    seq_start_pos = seqline.rfind(' ')
    fh.readline()
    fh.readline()
    mask = fh.readline()
    seqline = seqline[seq_start_pos:]
    mask = mask[seq_start_pos:]
    seqline = ''.join([seqline[x] for x in range(0, len(seqline)) if mask[x] == 'x'])
    
    if (seqline.count('-') / float(len(seqline))) > 0.5: # Limit to less than half gaps
        return None
    
    return seqline
