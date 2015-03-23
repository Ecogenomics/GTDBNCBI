###############################################################################
#
# hmmer.py - runs HMMER and provides functions for parsing output
#
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

from re import split as re_split

class FormatError(BaseException):
    pass

class HMMERParser():
    def __init__(self, fileHandle):
        self.handle = fileHandle

    def next(self):
        return self.readHitsTBL()

    def readHitsTBL(self):
        """Process single hit in tblout format."""
        """
We expect line to look like:
NODE_110054_length_1926_cov_24.692627_41_3 -          Ribosomal_S9         PF00380.14   5.9e-48  158.7   0.0   6.7e-48  158.5   0.0   1.0   1   0   0   1   1   1   1 # 1370 # 1756 # 1 # ID=41_3;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None
        """
        while True:
            line = self.handle.readline().rstrip()
            try:
                if line[0] != '#' and len(line) != 0:
                    dMatch = re_split(r'\s+', line.rstrip())
                    if len(dMatch) < 19:
                        raise FormatError("Error processing line:\n%s" % (line))
                    refined_match = dMatch[0:18] + [" ".join([str(i) for i in dMatch[18:]])]
                    return HmmerHitTBL(refined_match)
            except IndexError:
                return None

class HmmerHitTBL():
    """Encapsulate a HMMER hit given in tblout format."""
    def __init__(self, values):
        if len(values) == 19:
            self.target_name = values[0]
            self.target_accession = values[1]
            self.query_name = values[2]

            self.query_accession = values[3]
            if self.query_accession == '-':
                self.query_accession = self.query_name

            self.full_e_value = float(values[4])
            self.full_score = float(values[5])
            self.full_bias = float(values[6])
            self.best_e_value = float(values[7])
            self.best_score = float(values[8])
            self.best_bias = float(values[9])
            self.exp = float(values[10])
            self.reg = int(values[11])
            self.clu = int(values[12])
            self.ov = int(values[13])
            self.env = int(values[14])
            self.dom = int(values[15])
            self.rep = int(values[16])
            self.inc = int(values[17])
            self.target_description = values[18]

    def __str__(self):
        return "\t".join(
            [self.target_name,
            self.target_accession,
            self.query_name,
            self.query_accession,
            str(self.full_e_value),
            str(self.full_score),
            str(self.full_bias),
            str(self.best_e_value),
            str(self.best_score),
            str(self.best_bias),
            str(self.exp),
            str(self.reg),
            str(self.clu),
            str(self.ov),
            str(self.env),
            str(self.dom),
            str(self.rep),
            str(self.inc),
            self.target_description
            ]
        )
