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

__prog_name__ = 'plot_genome_msa_stats.py'
__prog_desc__ = 'Create plots of genome completeness, genome contamination, and percent MSA.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2019'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import tempfile
import ntpath
import shutil

from biolib.plots.abstract_plot import AbstractPlot

from matplotlib.ticker import FuncFormatter

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(5 * y * 100) # 5 is the bin width
    return s

class Histogram(AbstractPlot):
    """Histogram plot."""

    def __init__(self, options):
        AbstractPlot.__init__(self, options)
        
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

    def plot(self, plot_num, x, xlabel, ylabel, bins, color):
        """Plot histogram."""

        axis = self.fig.add_subplot(3, 1, plot_num)

        pdf, bins, patches = axis.hist(x, bins=bins, 
                                            normed=True,
                                            cumulative=False,
                                            color=color, 
                                            lw=0.5, 
                                            histtype='bar', 
                                            stacked=False)
        
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        
        axis.set_xticks(list(range(0, 101, 5)))
        
        formatter = FuncFormatter(to_percent)
        self.fig.gca().yaxis.set_major_formatter(formatter)
        
        self.prettify(axis)

        self.fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
        self.draw()
        

class PlotStats(object):
  """Create plots of genome completeness, genome contamination, and percent MSA."""

  def __init__(self):
    """Initialization."""
    pass

  def run(self, metadata_file, msa_info_file, genome_info_file, output_prefix):
    """Create plots of genome completeness, genome contamination, and percent MSA."""
    
    # read metadata
    comp = {}
    cont = {}
    type = {}
    with open(metadata_file) as f:
        headers = f.readline().strip().split('\t')
        
        comp_index = headers.index('checkm_completeness')
        cont_index = headers.index('checkm_contamination')
        type_index = headers.index('gtdb_type_designation')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            comp[gid] = float(line_split[comp_index])
            cont[gid] = float(line_split[cont_index])
            type[gid] = line_split[type_index]
            
    # read MSA info
    msa_perc = {}
    with open(msa_info_file) as f:
        headers = f.readline().strip().split('\t')
        
        msa_perc_index = headers.index('Amino acids (%)')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            msa_perc[gid] = float(line_split[msa_perc_index])
            
    # read species information
    sp = {}
    with open(genome_info_file) as f:
        headers = f.readline().strip().split('\t')
        
        sp_index = headers.index('Species')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            sp[gid] = (line_split[sp_index])
            
    # write out statistics to file
    fout = open(output_prefix + '_table.tsv', 'w')
    fout.write('Genome ID\tCompletenss (%)\tContamination (%)\tMSA completenss (%)\tSpecies\tGTDB type designation\n')
    for gid in msa_perc:
        fout.write('%s\t%.2f\t%.2f\t%.2f\t%s\t%s\n' % (gid, comp[gid], cont[gid], msa_perc[gid], sp[gid], type[gid]))
        
    # plot stats
    options = AbstractPlot.Options(6, 7.5, 10, 8, 300)
    
    hist = Histogram(options)
    hist.plot(1, list(msa_perc.values()), 'MSA completeness (%)', 'Genomes (%)', list(range(0, 101, 5)), 'blue')
    hist.plot(2, list(comp.values()), 'Completeness (%)', 'Genomes (%)', list(range(0, 101, 5)), 'blue')
    hist.plot(3, list(cont.values()), 'Contamination (%)', 'Genomes (%)', list(range(0, 101, 5)), 'blue')
    hist.save_plot(output_prefix + '.png')

if __name__ == '__main__':
  print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
  print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('metadata_file', help='GTDB metadata file')
  parser.add_argument('msa_info_file', help='MSA statistics (genome_msa_stats.tsv)')
  parser.add_argument('genome_info_file', help='Genome species information (gids_bac_canonical.lst)')
  parser.add_argument('output_prefix', help='output prefix')

  args = parser.parse_args()

  try:
    p = PlotStats()
    p.run(args.metadata_file, args.msa_info_file, args.genome_info_file, args.output_prefix)
  except SystemExit:
    print("\nControlled exit resulting from an unrecoverable error or warning.")
  except:
    print("\nUnexpected error:", sys.exc_info()[0])
    raise
