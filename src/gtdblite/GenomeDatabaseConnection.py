import psycopg2 as pg

from gtdblite import Config


class GenomeDatabaseConnection(object):
    def __init__(self):
        self.conn = None

    # Opens a connection to the PostgreSQL database
    #
    # Returns:
    #   No return value.
    def MakePostgresConnection(self):
        db_name = Config.GTDB_DB_NAME
        conn_string = "dbname=%s user=%s host=%s password=%s" % (
            db_name, Config.GTDB_USERNAME,
            Config.GTDB_HOST, Config.GTDB_PASSWORD
        )
        self.conn = pg.connect(conn_string)

    # Function: ClosePostgresConnection
    # Closes an open connection to the PostgreSQL database.
    #
    # Returns:
    #   No return value.
    def ClosePostgresConnection(self):
        self.conn.close()
        self.conn = None

    # Function: IsPostgresConnectionActive
    # Check if the connection to the PostgreSQL database is active.
    #
    # Returns:
    #   True if connection is active, False otherwise
    def IsPostgresConnectionActive(self):
        if self.conn is not None:
            cur = self.conn.cursor()
            try:
                cur.execute("SELECT count(*) from users")
            except:
                return False
            cur.close()
            return True
        else:
            return False

    # Convenience methods to the pg connection
    def commit(self):
        return self.conn.commit()

    def rollback(self):
        return self.conn.rollback()

    def cursor(self):
        return self.conn.cursor()