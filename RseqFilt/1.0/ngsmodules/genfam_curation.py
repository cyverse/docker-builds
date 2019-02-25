import argparse
import psycopg2


# connect to psql db
try:
    con = psycopg2.connect(database='plant_db', user='renesh', password='927122', host="127.0.0.1", port="5432")
except ValueError:
    print("Not able to connect database")

cur = con.cursor()
cur.execute("SELECT * from ath_gene_fam_group_1")
record = cur.fetchall()
genelist = []

for each_rec in record:
    fam = each_rec[0]
    gene = each_rec[1]
    gene = [x.upper() for x in gene]
    if 'MYB' in fam:
        print(fam)
        genelist.append(gene)
        gene = set(gene)
        print(len(gene))


# print(genelist)
# common = set(genelist[0]) & set(genelist[1]) & set(genelist[2])
# common = set(genelist[0]) & set(genelist[1])
# print("\n")
# print(common,"yes")
# print(len(common))
# common = list(set(genelist[0]).intersection(genelist[1]))
# print("\n")
# print(common,"yes")
# print(len(common))