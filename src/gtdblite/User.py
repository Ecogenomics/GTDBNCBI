#  Class: User
#  A class that represents a user of the database
class User(object):
    
    # Constuctor: User
    # Initialises the object
    def __init__(self):
        pass
    
    @classmethod
    def createUser(cls, userId, userName, roleId):
        C = cls()
        C.userId = userId
        C.userName = userName
        C.roleId = roleId
        C.rootUser = False
        return C
    
    @classmethod
    def createRootUser(cls):
        C = cls()
        C.rootUser = True
        return C
    
    # Function: getUserName
    # Returns the username of the user
    #
    # Returns:
    #   The username of the User as a string
    def getUserName(self):
        return self.userName
    
    # Function: getUserId
    # Returns the id of the user (as stored in PostgreSQL database)
    def getUserId(self):
        return self.userId
    
    # Function: getRoleId
    # Returns the role id of this user
    def getRoleId(self):
        return self.roleId
    
    # Function: isRootUser
    # Check if this is the rootUser
    def isRootUser(self):
        return self.rootUser        
        