import sqlite3

# how to incorporate taxonomy db?

def get_genome_by_taxid(db, taxid):
    try:
        with sqlite3.connect(db) as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT taxid, species, fasta_path FROM genomes WHERE taxid =?", (taxid,)
            )
            res = cur.fetchall()
            return res
    except sqlite3.Error as e:
        print(e)

def get_all_genomes(db):
    try:
        with sqlite3.connect(db) as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT * FROM genomes "
            )
            res = cur.fetchall()
            return res
    except sqlite3.Error as e:
        print(e)

if __name__ == "__main__":
    db_fn = "./data/collembola_data.db"
    print(get_genome_by_taxid(db_fn, 39272))
    #print(get_all_genomes(db_fn))