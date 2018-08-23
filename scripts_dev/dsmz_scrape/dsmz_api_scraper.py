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
import requests
import io

from requests.auth import HTTPBasicAuth
from unidecode import unidecode

__prog_name__ = 'dsmz_api_scraper.py'
__prog_desc__ = 'Produce metadata files describing type genera, species, and strains according to DSMZ'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2017'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'


class PNUClient(object):
    # ===============================================================================
    # REST Client for PNU Web Services.
    #
    # The attributes of the class are:
    # * ``headers`` -- sets the content-type of the HTTP request-header to json
    # * ``credentials`` -- attaches the username and the password to HTTPBasicAuth for using with the `requests` library
    # ===============================================================================

    headers = {'Accept': 'application/json'}
    USERNAME = 'To SET'
    PASSWORD = 'TO SET'
    credentials = HTTPBasicAuth(USERNAME, PASSWORD)

    def getGenera(self, outfile, urlreq=None):
        # to get all genera , but only for the first page.
        # to consider the other pages, you have to change the url to
        # 'https://bacdive.dsmz.de/api/pnu/genus/?page=2' etc.

        genus_type_species_dict = {}

        if urlreq is None:
            response = requests.get(
                'https://bacdive.dsmz.de/api/pnu/genus/', headers=self.headers, auth=self.credentials)
        else:
            response = requests.get(
                urlreq, headers=self.headers, auth=self.credentials)

        if response.status_code == 200:
            results = response.json()
            listgenus = results.get("results")
            urlreq = results.get("next")
            for item in listgenus:
                if item.get("label") is not None and item.get('authors') is not None and item.get('taxon') is not None:
                    outfile.write('g__{0}\t\t{1}\n'.format(
                        item.get("label"), item.get('authors') + ',' + item.get('taxon')))
                else:
                    outfile.write('g__{0}\t\t\t\n'.format(item.get("label")))

                if item.get('type_species') is not None:
                    genus_type_species_dict[item.get(
                        'type_species')] = item.get('label')
            if results.get("next") is not None:
                self.getGenera(outfile, urlreq)
#             OUTPUT:
#             object of type 'dict' with the fields 'count', 'previous', 'results', 'next'
# the different genera in field 'results' are separated by ',' e.g.
# {genus1},{genus2},{genus3},
            return genus_type_species_dict

    def getSpecies(self, outfile_species, outfile_strains, dictgenus, urlreq=None):
        # to get a list of all species , but only for the first page.
        # to consider the other pages, you have to change the url to
        # 'https://bacdive.dsmz.de/api/pnu/species/?page=2' etc.

        if urlreq is None:
            response = requests.get(
                'https://bacdive.dsmz.de/api/pnu/species/', headers=self.headers, auth=self.credentials)
        else:
            response = requests.get(
                urlreq, headers=self.headers, auth=self.credentials)

        if response.status_code == 200:
            results = response.json()
            listspe = results.get("results")
            urlreq = results.get("next")
            for item in listspe:
                # if 'subsp.' in item.get("label"):
                #    continue
                if item.get("type_strain") is not None:
                    list_strains = item.get("type_strain")
                    for st in item.get("type_strain"):
                        if " " in st:
                            list_strains.append(st.replace(" ", ''))
                    outfile_strains.write('{0}\t{1}\n'.format(
                        item.get("species"), "=".join(list_strains)))
#                else:
#                    outfile_strains.write('{0} \n'.format(item.get("species")))

                label = ''
                species_authority = ''
                type_spe = ''
                if item.get("label") is not None:
                    label = 's__' + item.get("label")
                if item.get('authors') is not None and item.get('taxon') is not None:
                    species_authority = item.get(
                        'authors') + ',' + item.get('taxon')
                    species_authority = unidecode(species_authority)
                if item.get("label") in dictgenus:
                    type_spe = 'g__' + dictgenus.get(item.get('label'))
                outfile_species.write('{0}\t{1}\t{2}\n'.format(
                    label, type_spe, species_authority))
            if results.get("next") is not None:
                self.getSpecies(outfile_species,
                                outfile_strains, dictgenus, urlreq)

#             OUTPUT:
#             object of type 'dict' with the fields 'count', 'previous', 'results', 'next'
# the species in field 'results' are separated by ',' e.g.
# {species1},{species2},{species3},etc
            return results

    def run(self, outdir, username, password):
        #         Run the client

        outfile_genera = open(os.path.join(outdir, 'dsmz_genera.tsv'), 'w')
        outfile_genera.write(
            "dsmz_genus\tdsmz_type_genus\tdsmz_genus_authority\n")
        dictgenus = self.getGenera(outfile_genera)
        outfile_genera.close()

        outfile_species = io.open(os.path.join(
            outdir, 'dsmz_species.tsv'), 'wb')
        outfile_species.write(
            "dsmz_species\tdsmz_type_species\tdsmz_species_authority\n")
        outfile_strains = open(os.path.join(outdir, 'dsmz_strains.tsv'), 'w')
        outfile_strains.write("dsmz_species\tdsmz_strains\n")
        self.getSpecies(outfile_species, outfile_strains, dictgenus)
        outfile_species.close()
        outfile_strains.close()


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('--username')
    parser.add_argument('--password')

    args = parser.parse_args()

    try:
        p = PNUClient()
        p.run(args.output_dir, args.username, args.password)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
