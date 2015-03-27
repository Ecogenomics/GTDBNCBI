###############################################################################
#
# resultsParser.py - Parse and output results.
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

def vetHit(model, hit):
    """Check if hit meets required thresholds."""

    # preferentially use model specific bit score thresholds, before
    # using the user specified e-value and length criteria

    # Give preference to the gathering threshold unless the model
    # is marked as TIGR (i.e., TIGRFAM model)

    if model.nc != None and 'TIGR' in model.acc:
        if model.nc[0] <= hit.full_score:
            return True
    elif model.ga != None:
        if model.ga[0] <= hit.full_score:
            return True
    elif model.tc != None:
        if model.tc[0] <= hit.full_score:
            return True
    elif model.nc != None:
        if model.nc[0] <= hit.full_score:
            return True
    else:
        alignment_length = float(hit.ali_to - hit.ali_from)
        length_perc = alignment_length / float(hit.query_length)
        if length_perc >= 0.75:
            return True

    return False

