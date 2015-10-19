import os
import time
import multiprocessing


def runProdigalMulti(nprocs=None, nogene_list=None):
    def worker(nogene_sublist, out_dict):
        """ The worker function, invoked in a process. 'genome_dict' is a
            dictionnary of genes to align. The results are placed in
            a dictionary that's pushed to a queue.
        """

        for key, value in nogene_sublist:
            if value.get("gene_path") is None:
                prodigal_tmp_dir = MarkerCalculation.RunProdigalOnGenomeFasta(value.get("fasta_path"))
                value["gene_path"] = os.path.join(prodigal_tmp_dir, "genes.faa")
                out_dict[key] = value

        return True

    # Each process will get 'chunksize' nums and a queue to put his out
    # dict into
    manager = multiprocessing.Manager()
    out_dict = manager.dict()
#    chunksize = int(math.ceil(len(genome_dict) / float(nprocs)))
    procs = []

    for item in splitchunks(nogene_list, nprocs):
        p = multiprocessing.Process(
                target=worker,
                args=(item, out_dict))
        procs.append(p)
	p.start()

    # Collect all results into a single result dict. We know how many dicts
    # with results to expect.
    resultlist = []
    while out_q.empty():
	   time.sleep(1)

    for item in splitchunks(nogene_list, nprocs):
        out_dict_get = out_dict.get()
        print "outqget = {0}".format(out_dict_get)
        resultlist.append(out_dict_get)

    # Wait for all worker processes to finish
    for p in procs:
        p.join()

    return resultlist
