import shutil
import os
import sys

from gtdblite.GenomeDatabaseConnection import GenomeDatabaseConnection
from gtdblite.Exceptions import GenomeDatabaseError


class MetadataManager(object):

    def __init__(self):
        self.conn = GenomeDatabaseConnection()
        self.errorMessages = []
        self.warningMessages = []

    #
    # Group: General Functions
    #
    # Function: ReportError
    # Sets the last error message of the database.
    #
    # Parameters:
    #     msg - The message to set.
    #
    # Returns:
    #   No return value.
    def ReportError(self, msg):
        self.errorMessages.append(str(msg))

    def GetErrors(self):
        return self.errorMessages

    def ClearErrors(self):
        self.errorMessages = []

    def ReportWarning(self, msg):
        self.warningMessages.append(str(msg))

    def GetWarnings(self):
        return self.warningMessages

    def ClearWarnings(self):
        self.warningMessages = []

    # Function: viewMetadata
    #
    # Lists all metadata field available in the database
    #
    # Returns:
    # Print lists
    def viewMetadata(self):
        try:
            self.conn.MakePostgresConnection()
            cur = self.conn.cursor()

            cur.execute("SELECT * FROM view_list_meta_columns")
            print "\t".join(("Table", "Field", "Description"))
            for tup in cur.fetchall():
                print "\t".join(list(tup))
            cur.close()
            self.conn.ClosePostgresConnection()
        except GenomeDatabaseError as e:
            raise self.ReportError(e.message)

    def exportMetadata(self, path):
        try:
            self.conn.MakePostgresConnection()
            cur = self.conn.cursor()
            query = "SELECT exportMeta('{0}')".format(path)
            cur.execute(query)
            print "Export Successful"

            cur.close()
            self.conn.ClosePostgresConnection
        except GenomeDatabaseError as e:
            raise self.ReportError(e.message)
