#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
import argparse
import urllib
import xml.etree.ElementTree as ET
from datetime import timedelta,datetime
import re
from multiprocessing.pool import ThreadPool as Pool
import multiprocessing
import time

__prog_name__ = 'straininfo_scraper_api.py'
__prog_desc__ = 'Produce metadata files describing type strains according to Straininfo'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2017'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'


class PNUClient(object):
    
    def __init__(self):
        self.threads =32
        self.pool_size = 32  # your "parallelness"
        self.pool = Pool(self.pool_size)

    def worker(self, intpage,out_q):
        strain_dict ={}       
        list_strains = []
        textname = ''
        url = 'http://www.straininfo.net/taxa/{0}'.format(intpage)
        print url
        
        try:
            testfile = urllib.URLopener()
            testfile.retrieve(url, "toparse{0}.html".format(intpage))
            
            with open("toparse{0}.html".format(intpage),'r') as htmlpagespe:
                spename = False
                species = ''
                for line in htmlpagespe:
                    if spename and len(species) == 0 :
                        matchObj = re.match('\s+<td class="value"><span class="speciesname"><em>([^<]+)<\/em><\/span><\/td>',line)
                        if matchObj:
                            print matchObj.group(1)
                            species = matchObj.group(1)
                    elif not spename and len(species) == 0 :
                        matchObj = re.match('\s+<tr><td class="option">species<\/td>',line)
                        if matchObj:
                            spename =True
                    else:
                        matchObj = re.match("\s*<div class='popup'>([^<]+)<strong>type strain<\/strong> of:<br\/>",line)
                        if matchObj:
                            strain =  matchObj.group(1).replace(" is ","")
                            list_strains.append(strain)
                            list_strains.append(strain.replace(" ",""))
                
            os.remove("toparse{0}.html".format(intpage))
            out_q.put((species,"=".join(set(list_strains))))


        except IOError:
            print 'url does not exist.'
            

        return True
        
    def run(self,outfile):
        outf = open(outfile,'w')
        outf.write("straininfo_strains\n")
        full_dict = {}
        manager = multiprocessing.Manager()
        out_q = manager.Queue()
        workers = [self.pool.apply_async(self.worker, (i,out_q)) for i in range(400000)]

        
        # Collect all results into a single result dict. We know how many dicts
        # with results to expect.
        while out_q.empty():
            time.sleep(1)

        self.pool.close()
        self.pool.join()
        
        while not out_q.empty():
            info = out_q.get()
            if info[0] != '' and info[1] != '':
                outf.write("{0} {1}\n".format(info[0],info[1]))
        outf.close()
            

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_file', help='output file')

    args = parser.parse_args()

    try:
        p = PNUClient()
        p.run(args.output_file)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
