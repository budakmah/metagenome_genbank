from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys, os

path=os.getcwd()



gb_files=input("input the path for genbank folder :       ")
choice=input("""Please enter your choice for extract
            for all (tRNA,rRNA,CDS): 1
            for CDS : 2
            for CDS and rRNA : 3
            for Translation : 4
            your choice:        """)




if choice not in ["1","2","3","4"]:
    print("Your choice {} is not suitable".format(choice))
    exit()


cds_n={'cytochrome c oxidase subunit 2': 'COX2',
     'cytochrome c oxidase subunit 3': 'COX3',
     'NADH deshydrogenase subunit 4L': 'ND4L',
     'ATP synthase F0 subunit 6': 'ATP6',
     'NADH dehydrogenase subunit 4l': 'ND4L',
     'NADH deshydrogenase subunit 1': 'ND1',
     'NADH deshydrogenase subunit 3': 'ND3',
     'cytochrome c oxidase subunit IIb': 'COX2b',
     'NADH deshydrogenase subunit 2': 'ND2',
     'ATP synthase FO subunit 8': 'ATP8',
     'ATPase subunit 8': 'ATP8',
     'cytochrome c oxidase subunit IIa': 'COX2a',
     'ATP synthase subunit 8': 'ATP8',
     'cytochrome oxidase subunit 2': 'COX2',
     'ATP synthase FO subunit 6': 'ATP6',
     'ATP F0 synthase subunit 6': 'ATP6',
     'ATP F0 synthase subunit 8': 'ATP8',
     'ATP synthase 6': 'ATP6',
     'NADH dehydrogenase subunit 2': 'ND2',
     'cytochrome c oxidase subunit I': 'COX1',
     'cytochrome oxidase subunit III': 'COX3',
     'NADH dehydrogenase subunit 6': 'ND6',
     'cytochrome oxidase subunit 3': 'COX3',
     'ATP8': 'ATP8',
     'NADH dehydrogenase subunit 4': 'ND4',
     'ATP synthase subunit 6': 'ATP6',
     'cytochrome c oxidase subunit 1': 'COX1',
     'NADH dehydrogenase subunit 5': 'ND5',
     'ATPase subunit 6': 'ATP6',
     'cytochrome oxidase subunit 1': 'COX1',
     'NADH dehydrogenase subunit 3': 'ND3',
     'cytochrome oxidase subunit I': 'COX1',
     'cytochrome oxidase subunit II': 'COX2',
     'NADH deshydrogenase subunit 4': 'ND4',
     'NADH dehydrogenase subunit 1': 'ND1',
     'cytochrome c oxidase subunit III': 'COX3',
     'NADH deshydrogenase subunit 5': 'ND5',
     'NADH deshydrogenase subunit 6': 'ND6',
     'cytochrome b': 'CYTB',
     'NADH dehydrogenase subunit 4L': 'ND4L',
     'ATP synthase F0 subunit 8': 'ATP8',
     'ATP synthase 8': 'ATP8',
     'apocytochrome b': 'CYTB',
     'cytochrome c oxidase subunit II': 'COX2'}
rrna_n={'s-rRNA': '12S',
     'large subunit ribosomal RNA': '16S',
     '16S ribosomal RNA': '16S',
     'srRNA': '12S',
     '12S ribosomal RNA': '12S',
     'l6S ribosomal RNA': '16S',
     'small ribosomal RNA': '12S',
     'large ribosomal RNA subunit': '16S',
     'lrRNA': '16S',
     'large ribosomal RNA': '16S',
     'l2S ribosomal RNA': '12S',
     'l-rRNA': '16S',
     'small ribosomal RNA subunit': '12S',
     'small subunit ribosomal RNA': '12S'}



def trna(file):
    gb=SeqIO.read(file,"genbank")
    organism=[]
    if gb.annotations["organism"] not in organism:
        organism.append(gb.annotations["organism"])
        for f in gb.features:
            if f.type == 'tRNA':
                fi=os.path.join(trna_path,f.qualifiers['product'][0]+".fna")
                with open(fi,"a") as out:
                    out.write(">"+gb.annotations["organism"].replace(" ", "_")+"\n"+str(f.location.extract(gb).seq)+"\n")

def cds(file):
    gb=SeqIO.read(file,"genbank")
    organism=[]
    if gb.annotations["organism"] not in organism:
        organism.append(gb.annotations["organism"])
        for f in gb.features:
            if f.type == 'CDS':
                fi=os.path.join(cds_path,cds_n[f.qualifiers["product"][0]]+".fna")
                with open(fi,"a") as out:
                    out.write(">"+gb.annotations["organism"].replace(" ", "_")+"\n"+str(f.location.extract(gb).seq)+"\n")
def rrna(file):
    gb=SeqIO.read(file,"genbank")
    organism=[]
    if gb.annotations["organism"] not in organism:
        organism.append(gb.annotations["organism"])
        for f in gb.features:
            if f.type == 'rRNA':
                fi=os.path.join(rrna_path,rrna_n[f.qualifiers['product'][0]]+".fna")
                with open(fi,"a") as out:
                    out.write(">"+gb.annotations["organism"].replace(" ", "_")+"\n"+str(f.location.extract(gb).seq)+"\n")

def aa(file):
    gb=SeqIO.read(file,"genbank")
    organism=[]
    if gb.annotations["organism"] not in organism:
        organism.append(gb.annotations["organism"])
        for f in gb.features:
            if f.type == 'CDS':
                fi=os.path.join(cds_path,"aa_"+cds_n[f.qualifiers["product"][0]]+".faa")
                with open(fi,"a") as out:
                    out.write(">"+gb.annotations["organism"].replace(" ", "_")+"\n"+f.qualifiers['translation'][0]+"\n")


if choice == "1":
    os.system("mkdir {}".format(os.path.join(path, "tRNAs")))
    trna_path=os.path.join(path, "tRNAs")
    os.system("mkdir {}".format(os.path.join(path, "rRNAs")))
    rrna_path=os.path.join(path, "rRNAs")
    os.system("mkdir {}".format(os.path.join(path, "CDSs")))
    cds_path=os.path.join(path, "CDSs")
    for gb in os.listdir(gb_files):
        if ".gb" in gb:
            trna(os.path.join(gb_files,gb))
            rrna(os.path.join(gb_files,gb))
            cds(os.path.join(gb_files,gb))

if choice == "2":
    os.system("mkdir {}".format(os.path.join(path, "CDSs")))
    cds_path=os.path.join(path, "CDSs")
    for gb in os.listdir(gb_files):
        if ".gb" in gb:
            cds(os.path.join(gb_files,gb))

if choice == "3":
    os.system("mkdir {}".format(os.path.join(path, "CDSs")))
    cds_path=os.path.join(path, "CDSs")
    os.system("mkdir {}".format(os.path.join(path, "rRNAs")))
    rrna_path=os.path.join(path, "rRNAs")
    for gb in os.listdir(gb_files):
        if ".gb" in gb:
            rrna(os.path.join(gb_files,gb))
            cds(os.path.join(gb_files,gb))

if choice == "4":
    os.system("mkdir {}".format(os.path.join(path, "CDSs")))
    cds_path=os.path.join(path, "CDSs")
    for gb in os.listdir(gb_files):
        if ".gb" in gb:
            aa(os.path.join(gb_files,gb))
