import sqlite3


def create_sqlite_database(db):
    """create a database connection to an SQLite database"""
    conn = None
    try:
        conn = sqlite3.connect(db)
        print(sqlite3.sqlite_version)
    except sqlite3.Error as e:
        print(e)
    finally:
        if conn:
            conn.close()

def create_tables(db):
    sql_statements = [
        """CREATE TABLE IF NOT EXISTS genomes (
				taxid INTEGER NOT NULL, 
				species text NOT NULL, 
				fasta_path TEXT,
				status TEXT 
		);""",
        """CREATE TABLE IF NOT EXISTS BUSCOs (
				taxid INTEGER NOT NULL, 
				species text NOT NULL, 
				directory text
		);""",
        """CREATE TABLE IF NOT EXISTS annotations (
				taxid INTEGER NOT NULL,
				species text NOT NULL, 
				repeat_dir text,
				GTF text,
				VCF text
		);""",
        """CREATE TABLE IF NOT EXISTS reseq_reads (
				taxid INTEGER NOT NULL,
				species text NOT NULL, 
				forward_reads text,
				reverse_reads text
		);""",
        """CREATE TABLE IF NOT EXISTS rnaseq_reads (
				taxid INTEGER NOT NULL,
				species text NOT NULL, 
				forward_reads text,
				reverse_reads text
		);""",
    ]
    # return cur.lastrowid from genome table as record in other tables

    # create a database connection
    try:
        with sqlite3.connect(db) as conn:
            cursor = conn.cursor()
            for statement in sql_statements:
                cursor.execute(statement)

            conn.commit()
    except sqlite3.Error as e:
        print(e)


if __name__ == "__main__":
    db_fn = "./data/collembola_data.db"
    create_sqlite_database(db_fn)
    create_tables(db_fn)
