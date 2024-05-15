#!/usr/bin/env python3

from Bio import SeqIO
import sqlite3
import numpy as np
import pandas as pd


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def get_sequence(idx, fasta, idlist):
    "Get sequences with ids in idlist in an indexed fasta file"
    records = SeqIO.index_db(idx, fasta, "fasta")
    return [records[str(i + 1)] for i in idlist]  # fasta records start indexing at 1.


def get_procap(procap, idlist):
    "Get procap data with indices in idlist in a procap database"
    conn = sqlite3.connect(procap)
    idstring = "(%s)" % ", ".join([str(id) for id in idlist])
    q = "SELECT * FROM procap where ROWID in %s" % idstring
    q_result = pd.read_sql_query(q, conn)
    return np.array(q_result)[0][1:]  # clean output and slice off index
