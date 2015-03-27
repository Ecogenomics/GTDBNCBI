import tempfile
import os
import shutil
import re

from gtdblite.Exceptions import GenomeDatabaseError

from checkm_static import prodigal, hmmer, resultsParser

from simplehmmer import hmmmodelparser



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


def CalculateBestMarkerOnProdigalDir(output_prefix, marker_hmm_file, prodigal_dir):
    
    hmmsearch_result_filepath = os.path.join(prodigal_dir, "%s.search" % output_prefix)
    
    # Do a hmmsearch to find the best hit
    os.system("hmmsearch --tblout %s %s %s > /dev/null" % (
        hmmsearch_result_filepath,
        marker_hmm_file,
        os.path.join(prodigal_dir, 'genes.faa')
    ))
    
    model = hmmmodelparser.HmmModelParser(marker_hmm_file).parse().next()
    
    blank_result = model.leng * '-'
    
    best_hit = FindBestMarkerInHMMSearchResult(hmmsearch_result_filepath)
    if best_hit is None:
        return blank_result

    if not resultsParser.vetHit(model, best_hit):
        return blank_result
    
    hmmalign_result_filepath = os.path.join(prodigal_dir, "%s.aligned" % (output_prefix,))
    
    # Do a hmmalign and retrieve the aligned best hit given by hmmsearch
    os.system("hmmalign --outformat Pfam -o %s %s %s" % (
        hmmalign_result_filepath,
        marker_hmm_file,
        os.path.join(prodigal_dir, 'genes.faa')
    ))
    
    result = GetAlignedMarker(best_hit.target_name, hmmalign_result_filepath)
    if result is None:
        return blank_result
    
    if len(result) != model.leng:
        raise Exception("Result length doesn't equal the model length")
    
    return result
    
    
def FindBestMarkerInHMMSearchResult(hmmsearch_result_file):
    
    fh = open(hmmsearch_result_file)
    parser = hmmer.HMMERParser(fh)
    best_score = 0
    best_hit = None
    while True:
        hit = parser.next()
        if not hit:
            break
        if hit.full_score > best_score:
            best_score = hit.full_score
            best_hit = hit    

    return best_hit

def GetAlignedMarker(hit_name, hmmalign_result_filepath):

    fh = open(hmmalign_result_filepath)
    for line in fh:
        if line[0:len(hit_name)] == hit_name:
            splitline = line.rsplit(" ", 1)
            return re.sub('[^A-Z\-]','', splitline[-1])

    return None
