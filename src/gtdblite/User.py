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


#  Class: User
#  A class that represents a user of the database
class User(object):

    # Constuctor: User
    # Initialises the object
    def __init__(self):
        pass

    @classmethod
    def createUser(cls, userId, username, rolename, role_id):
        C = cls()
        C.userId = userId
        C.username = username
        C.rolename = rolename
        C.role_id = role_id
        C.rootUser = False
        return C

    @classmethod
    def createRootUser(cls, elevatedFromUsername):
        C = cls()
        C.rootUser = True
        C.elevatedFromUsername = elevatedFromUsername
        return C

    # Function: getUserName
    # Returns the username of the user
    #
    # Returns:
    #   The username of the User as a string
    def getUsername(self):
        return self.username

    # Function: getUserId
    # Returns the id of the user (as stored in PostgreSQL database)
    def getUserId(self):
        return self.userId

    # Function: getRolename
    # Returns the role id of this user
    def getRolename(self):
        return self.rolename

    # Function: getRoleId
    # Returns the role id of this user
    def getRoleId(self):
        return self.roleId

    # Function: isRootUser
    # Check if this is the rootUser
    def isRootUser(self):
        return self.rootUser

    # Function: elevatedFromUsername
    # Get the username of the user before becoming root
    def getElevatedFromUsername(self):
        return self.elevatedFromUsername
