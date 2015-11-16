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

import logging

import psycopg2

from User import User
from Exceptions import GenomeDatabaseError


class UserManager(object):
    """Manages users in database."""

    def __init__(self, cur):
        """Initialize.

        Parameters
        ----------
        cur : psycopg2.cursor
            Database cursor.
        """

        self.logger = logging.getLogger()

        self.cur = cur
        self.currentUser = None

    # Function: UserLogin
    # Log a user into the database (make the user the current user of the database).
    #
    # Parameters:
    #     username - The username of the user to login
    #
    # Returns:
    # Returns a User calls object on success (and sets the GenomeDatabase
    # current user).
    def userLogin(self, username):
        try:
            self.cur.execute("SELECT users.id, user_roles.id, user_roles.name "
                        "FROM users, user_roles " +
                        "WHERE users.role_id = user_roles.id " +
                        "AND users.username = %s", (username,))

            result = self.cur.fetchone()

            if not result:
                raise GenomeDatabaseError("User not found: %s" % username)

            (user_id, role_id, rolename) = result
            self.currentUser = User.createUser(user_id, username, rolename, role_id)

        except GenomeDatabaseError as e:
            raise e

        return self.currentUser

    # Function: RootLogin
    # Log a user into the database as a root user (make the root user the current user of the database). Check
    # if the current user has permissions to do this.
    #
    # Parameters:
    #     username - The username of the user to login
    #
    # Returns:
    # Returns a User calls object on success (and sets the GenomeDatabase
    # current user).
    def rootLogin(self, username):
        try:
            if not self.conn.IsPostgresConnectionActive():
                raise GenomeDatabaseError(
                    "Unable to establish database connection")

            query = "SELECT id, has_root_login FROM users WHERE username = %s"
            self.cur.execute(query, [username])
            result = self.cur.fetchone()
            self.cur.close()

            if result:
                (_userid, has_root_login) = result
                if not has_root_login:
                    raise GenomeDatabaseError(
                        "You do not have sufficient permissions to logon as the root user.")

                self.currentUser = User.createRootUser(username)
            else:
                raise GenomeDatabaseError("User %s not found." % username)

        except GenomeDatabaseError as e:
            raise e

        return self.currentUser

    # Function: AddUser

    # Add a new user to the database.
    #
    # Parameters:
    #     username - The username of the user to login
    #     usertype - The role of the new user
    #
    # Returns:
    #   True on success, False otherwise.
    def addUser(self, username, rolename=None, has_root=False):
        try:
            if rolename is None:
                rolename = 'user'

            if (not self.currentUser.isRootUser()):
                if has_root:
                    raise GenomeDatabaseError(
                        "Only the root user may grant root access to new users.")

                if rolename == 'admin':
                    raise GenomeDatabaseError(
                        "Only the root user may create admin accounts.")

                if not(self.currentUser.getRolename() == 'admin' and rolename == 'user'):
                    raise GenomeDatabaseError(
                        "Only admins (and root) can create user accounts.")

            self.cur.execute(
                "SELECT username from users where username = %s", (username,))

            if len(self.cur.fetchall()) > 0:
                raise GenomeDatabaseError(
                    "User %s already exists in the database." % username)

            self.cur.execute("INSERT into users (username, role_id, has_root_login) (" +
                        "SELECT %s, id, %s " +
                        "FROM user_roles " +
                        "WHERE name = %s)", (username, has_root, rolename))

            self.conn.commit()

        except GenomeDatabaseError as e:
            self.conn.rollback()
            raise e
        except:
            self.conn.rollback()
            raise

        return True

    def editUser(self, username, rolename=None, has_root=None):
        try:
            if (not self.currentUser.isRootUser()):
                raise GenomeDatabaseError(
                    "Only the root user may edit existing accounts.")

                # The following may be useful in the future if roles change, but at the moment,
                # only the root user can make any meaningful user edits
                """
                if has_root is not None:
                    raise GenomeDatabaseError("Only the root user may edit the root access of users.")

                if rolename == 'admin':
                    raise GenomeDatabaseError("Only the root user may create admin accounts.")

                cur.execute("SELECT users.id, user_roles.id, user_roles.name "
                    "FROM users, user_roles " +
                    "WHERE users.role_id = user_roles.id " +
                    "AND users.username = %s", (username, ))

                result = cur.fetchone()

                if not result:
                    raise GenomeDatabaseError("User not found: %s" % username)

                (user_id, current_role_id, current_rolename) = result
                if current_rolename == 'admin':
                    raise GenomeDatabaseError("Only the root user may edit current admin accounts.")
                """

            conditional_queries = []
            params = []

            if rolename is not None:
                conditional_queries.append(
                    " role_id = (SELECT id from user_roles where name = %s) ")
                params.append(rolename)

            if has_root is not None:
                conditional_queries.append(" has_root_login = %s ")
                params.append(has_root)

            if params:
                self.cur.execute("UPDATE users " +
                            "SET " + ','.join(conditional_queries) + " "
                            "WHERE username = %s", params + [username])

            self.conn.commit()
        except GenomeDatabaseError as e:
            self.conn.rollback()
            raise e
        except:
            self.conn.rollback()
            raise

        return True
