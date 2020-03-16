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


def version():
    """Read software and NCBI version information from file.

    Returns
    -------
    str
        Software version.
    str
        Software, NCBI database, and GTDB database versions.
    """
    dir_src = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    with open(os.path.join(dir_src, 'VERSION')) as fh:
        software_version = fh.readline().strip()
        software_version = software_version[software_version.find('=') + 1:]

        ncbi_version = fh.readline().strip()
        ncbi_version = ncbi_version[ncbi_version.find('=') + 1:]

        gtdb_version = fh.readline().strip()
        gtdb_version = gtdb_version[gtdb_version.find('=') + 1:]

        taxonomy_version = fh.readline().strip()
        taxonomy_version = taxonomy_version[taxonomy_version.find('=') + 1:]

    return software_version, ncbi_version, gtdb_version, taxonomy_version
