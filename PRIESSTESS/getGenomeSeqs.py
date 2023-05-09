from collections import Counter
import glob, os, pandas
from gtfparse import read_gtf # ==1.3.0
os.chdir('/data/realServices/jorge/datos/alessiaGenes')

gtfFile = "/data/realServices/jorge/datos/genome_ref/Homo_sapiens.GRCh38.108.gtf"
gtfDF = read_gtf(gtfFile)

selectedGenesDFfiles = glob.glob("genelists/*")
selectedGenesDFfile = selectedGenesDFfiles[0]
ensembls = set()
for selectedGenesDFfile in selectedGenesDFfiles:
    selectedGenesDF = pandas.read_csv(selectedGenesDFfile,sep="\t")
    ensembls.update(selectedGenesDF.x.dropna())

subgtfDF = gtfDF.loc[gtfDF.gene_id.isin(ensembls)]
subgtfDF = subgtfDF.astype(str)

outBedsDir = 'beds'
os.makedirs(outBedsDir,exist_ok=True)
unmappedDFsDir = 'unmapgenes'
os.makedirs(unmappedDFsDir,exist_ok=True)

basicColumns = ['seqname', 'start', 'end']
infoColums = ['gene_id', 'feature', 'gene_version', 'gene_name','strand', 'frame', 'gene_source',
              'gene_biotype', 'transcript_id', 'transcript_version',
              'transcript_name', 'transcript_source', 'transcript_biotype', 'tag',
              'ccds_id', 'exon_number', 'exon_id', 'exon_version', 'protein_id',
              'protein_version', 'transcript_support_level']
selectedGenesDFfile = selectedGenesDFfiles[0]
importantFeatures = ['gene','utr']
for selectedGenesDFfile in selectedGenesDFfiles:
    selectedGenes = pandas.read_csv(selectedGenesDFfile,sep="\t").x.dropna()
    geneGTFdf = subgtfDF.loc[subgtfDF.gene_id.isin(selectedGenes)]
    selectedGenesDFnotMapped = set(selectedGenes) - set(subgtfDF.gene_id)

    with open(os.path.join(unmappedDFsDir,os.path.basename(selectedGenesDFfile)),'w') as selectedGenesDFnotMappedHNDL:
        selectedGenesDFnotMappedHNDL.write('\n'.join(list(selectedGenesDFnotMapped)))

    for importantFeature in importantFeatures:
        outGTFname = importantFeature+'_'+os.path.basename(selectedGenesDFfile.replace('.csv','.bed'))
        featureGTF = geneGTFdf.loc[geneGTFdf.feature.str.contains(importantFeature)]
        if featureGTF.shape[0] == 0:
            print("EMPTY")
            continue

        subfeatureGTF = featureGTF[basicColumns]
        subfeatureGTF.start = pandas.to_numeric(subfeatureGTF.start) - 200
        subfeatureGTF.end = pandas.to_numeric(subfeatureGTF.end) + 200
        subfeatureGTF['info'] = featureGTF[infoColums].apply(lambda x: '|'.join(x), axis=1)
        subfeatureGTF.to_csv(os.path.join(outBedsDir,outGTFname),sep="\t",index=False,header=False)

genefastasDir = 'genefastas'
os.makedirs(genefastasDir,exist_ok=True)

def bedtoolsGetfasta(bedFile):
    genome = "/data/realServices/jorge/datos/genome_ref/Homo_sapiens.GRCh38.dna.toplevel.fa"
    outFasta = bedFile.replace('beds','genefastas').replace('.bed','.fa')
    cmd = "bedtools getfasta -fi {} -bed {} -name -fo {}".format(genome,bedFile,outFasta)
    os.system(cmd)

bedFiles = glob.glob(outBedsDir+'/utr*')

from multiprocessing import Pool
with Pool(processes=6) as pool:
    pool.starmap(bedtoolsGetfasta, zip(bedFiles))
