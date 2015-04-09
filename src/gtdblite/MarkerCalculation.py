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
        
        # Rename the genes.faa fasta entries to remove the chance id parsing mixups, etc....
        modified_genes_filepath = os.path.join(prodigal_dir, 'genes_id_modified.faa')
        
        out_fh = open(modified_genes_filepath, 'wb')
        genes_fh = open(os.path.join(prodigal_dir, 'genes.faa'), 'rb')
        
        count = 1
        
        for line in genes_fh:
            
            if line[0] == ">":
                out_fh.write(">%i\n" % count)
                count += 1
                continue
            
            out_fh.write(line)
        
        out_fh.close()
        genes_fh.close()
        
        return prodigal_dir
    
    except:
        if prodigal_dir:
            shutil.rmtree(prodigal_dir)
        raise


def CalculateBestMarkerOnProdigalDir(output_prefix, marker_hmm_file, prodigal_dir):
    
    hmmsearch_result_filepath = os.path.join(prodigal_dir, "%s.search" % output_prefix)
    
    #TODO: Instead of using the prodigal dir, use a temp file.
    modified_genes_filepath = os.path.join(prodigal_dir, 'genes_id_modified.faa')
    
    # Do a hmmsearch to find the best hit
    os.system("hmmsearch --tblout %s %s %s > /dev/null" % (
        hmmsearch_result_filepath,
        marker_hmm_file,
        modified_genes_filepath
    ))
    
    model = hmmmodelparser.HmmModelParser(marker_hmm_file).parse().next()
    
    blank_result = model.leng * '-'
    
    best_hit = FindBestMarkerInHMMSearchResult(hmmsearch_result_filepath)
    if best_hit is None:
        return blank_result

    if not resultsParser.vetHit(model, best_hit):
        return blank_result
    
    # Create a file containing only the best hit
    writing_to_best_hit = False
    
    genes_fh = open(modified_genes_filepath, 'rb')
    
    best_hit_filepath = os.path.join(prodigal_dir, 'best_hit_%s_%s.faa' % (output_prefix, best_hit.target_name))
    best_hit_fh = open(best_hit_filepath, 'wb')
    
    for line in genes_fh:
        
        if line[0] == ">":
            writing_to_best_hit = False
        
        if line.startswith('>' + best_hit.target_name):
            writing_to_best_hit = True
            
        if writing_to_best_hit:
            best_hit_fh.write(line)
    
    best_hit_fh.close()

    # Do a hmmalign and retrieve the aligned best hit given by hmmsearch        
    
    hmmalign_result_filepath = os.path.join(prodigal_dir, "%s.aligned" % (output_prefix,))
    
    os.system("hmmalign --outformat Pfam -o %s %s %s" % (
        hmmalign_result_filepath,
        marker_hmm_file,
        best_hit_filepath
    ))
    
    result = GetAlignedMarker(best_hit.target_name, hmmalign_result_filepath)
    if result is None:
        return blank_result
    
    if len(result) != model.leng:
        raise Exception(
            "Result length doesn't equal the model length. Prodigal dir: %s. Output Prefix: %s. Result: %s. Model Length: %i" %
            (prodigal_dir, output_prefix, result, model.leng)
        )
    
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
    
    hit_seq = None
    mask_seq = None
    
    for line in fh:
        splitline = line.split(" ", 1)
        if splitline[0] == hit_name:
            rsplitline = line.rsplit(" ", 1)
            hit_seq = rsplitline[-1]
            continue
        if line[0:len("#=GC RF")] == "#=GC RF":
            rsplitline = line.rsplit(" ", 1)
            mask_seq = rsplitline[-1]
        
    if mask_seq is None:
        raise Exception("Unable to get mask from hmm align result file: %s" % hmmalign_result_filepath)
    
    if hit_seq is None:
        return None
    
    aligned_marker = ""
    for pos in xrange(0, len(mask_seq)):
        if mask_seq[pos] != 'x':
            continue
        aligned_marker += hit_seq[pos]

    return aligned_marker
