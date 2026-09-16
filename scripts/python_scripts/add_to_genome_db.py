import sqlite3
import sys
import pandas as pd


def add_genome(conn, project):
	sql = """ INSERT INTO projects(taxid, species, fasta_path, status)
			  VALUES(?,?,?,?) """
	cur = conn.cursor()
	cur.execute(sql, project)
	conn.commit()
	return cur.lastrowid


def add_BUSCO(conn, project):
	sql = """ INSERT INTO projects(taxid, species, directory)
			  VALUES(?,?,?) """
	cur = conn.cursor()
	cur.execute(sql, project)
	conn.commit()
	return cur.lastrowid


def add_annotation(conn, project):
	sql = """ INSERT INTO projects(taxid, species, repeat_dir, GTF, VCF)
			  VALUES(?,?,?,?,?) """
	cur = conn.cursor()
	cur.execute(sql, project)
	conn.commit()
	return cur.lastrowid


def add_reseq_reads(conn, project):
	sql = """ INSERT INTO projects(taxid, species, forward_reads, reverse_reads)
			  VALUES(?,?,?,?) """
	cur = conn.cursor()
	cur.execute(sql, project)
	conn.commit()
	return cur.lastrowid


def add_rnaseq_reads(conn, project):
	sql = """ INSERT INTO projects(taxid, species, forward_reads, reverse_reads)
			  VALUES(?,?,?,?) """
	cur = conn.cursor()
	cur.execute(sql, project)
	conn.commit()
	return cur.lastrowid


def add_spreadsheet(db, data, table, delimiter="\t"):
	if table == "genomes":
		colnames = ["taxid","species","fasta_path", "status"]

	elif table == "BUSCOs":
		colnames = ["taxid","species","directory"]

	elif table == "annotations":
		colnames = ["taxid","species","repeat_dir","GTF","VCF"]

	elif table == "reseq_reads":
		colnames = ["taxid","species","forward_reads","reverse_reads"]

	elif table == "rnaseq_reads":
		colnames = ["taxid","species","forward_reads","reverse_reads"]

	else:
		sys.exit("Invalid table name")

	try:
		with sqlite3.connect(db) as conn:
			df = pd.read_csv(data, delimiter=delimiter, names=colnames)
			df.to_sql(name=table, con=conn, if_exists="append", index=False)
			conn.commit()
	except sqlite3.Error as e:
		print(e)


def add_entries(db):
	try:
		with sqlite3.connect(db) as conn:
			entries = [
				(
					"X",
					"Y",
					"Z"
				),
				(
					"A",
					"B",
					"C",
				),
			]
			for entry in entries:
				entry_id = add_genome(conn, entry)
				print(f"Added entry with the id {entry_id}")

	except sqlite3.Error as e:
		print(e)


if __name__ == "__main__":
	db_fn = "./data/collembola_data.db"
	data_fn = "./data/db_input_sheets/genomes.tsv"
	# add_entries(fn)
	# add_spreadsheet(db_fn, data_fn, "genomes", "\t")
