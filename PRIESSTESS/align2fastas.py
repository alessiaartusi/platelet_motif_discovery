import glob, subprocess, os
from multiprocessing import Pool
os.chdir('/data/realServices/jorge/datos/alessiaGenes')

def launchjob(cmd,logfile='logfile'):
    cmd = list(map(str,cmd))
    print(' '.join(cmd)+'\n')
    try:
        process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        with open(logfile+'.log','a') as loghndl:
            loghndl.write(process.stdout.decode("utf-8"))
        print('Subprocess finished')
    except Exception as exception:
        with open(logfile+'.log','a') as loghndl:
            loghndl.write(process.stdout.decode("utf-8"))
            loghndl.write(exception)
        print(exception)

def createHISAT2idx(fasta,threads=8):
    outIDXdir = fasta.replace("genefastas","hisat2idx").replace(".fa","")
    cmd = ['/home/eidrian/Software/hisat2-2.2.1/hisat2-build', '-p', threads, fasta, outIDXdir]
    launchjob(cmd)

# 1 STEP CREATE IDXs
fastasDir = "genefastas"
fastas = glob.glob(fastasDir+"/utr*.fa"); fasta = fastas[0]
outIDXdir = "hisat2idx"; os.makedirs(outIDXdir,exist_ok=True)
with Pool(processes=4) as pool:
    pool.starmap(createHISAT2idx, zip(fastas))

def samSort(sam,threads=8):
    bamsorted = sam.replace('.sam','.sorted.bam')
    if os.path.exists(bamsorted) == False:
        cmd = ['/home/eidrian/Software/samtools-1.15/samtools', 'sort', '-@', threads, '-o', bamsorted, sam]
        launchjob(cmd)
        os.system('rm '+sam)
    return(bamsorted)

def customHISAT2align(refGenome,fastq1,threads=12):
    outbasedir = "hisat2align"
    outdir = os.path.join(outbasedir,os.path.basename(refGenome).replace('gene_',''))
    os.makedirs(outdir,exist_ok=True)
    outName = os.path.basename(fastq1).split('.')[0]
    sam = os.path.join(outdir,outName+'.sam')
    bamsorted = sam.replace('.sam','.sorted.bam')
    if os.path.exists(bamsorted):
        print('DONE: '+bamsorted)
        return(bamsorted)
    else:
        if os.path.exists(sam):
            print('DONE: '+sam)
            bamsorted = samSort(sam,threads)
            return(bamsorted)
        cmd = ['/home/eidrian/Software/hisat2-2.2.1/hisat2','-p', threads,'-x',refGenome,fastq1,'-S', sam]
        print('PROCESSING: '+sam+'\n')
        launchjob(cmd,sam)
        bamsorted = samSort(sam)
        return(' '.join(map(str,cmd)))


#fastqsSTR = ' '.join(fastqs)
#big_fastqDir = 'big_fastq'
#fastq1 = '{}/big.single.fastq'.format(big_fastqDir)
#os.makedirs(big_fastqDir,exist_ok=True)
#cmd = 'cat {} > {}'.format(fastqsSTR,fastq1)
#os.system(cmd)

fastqs = glob.glob('/data/realServices/jorge/datos/GSE68086/fastqs/*/*.fastq'); fastq1 = fastqs[0] ; len(fastqs)
refGenomes = glob.glob("hisat2idx/utr*")
refGenomes = list(set(refGenome.split(".")[0] for refGenome in refGenomes)); refGenome = refGenomes[0]
parameters = [(refGenome,fastq1) for refGenome in refGenomes for fastq1 in fastqs]
parameters = tuple(parameters)

with Pool(processes=2) as pool:
    bams = pool.starmap(customHISAT2align, parameters)

def getBAMalignedReads(bamFile):
    outfasta = bamFile.replace('hisat2align','alignedReads').replace('.sorted.bam','.fa')
    os.makedirs(os.path.dirname(outfasta),exist_ok=True)
    if os.path.exists(outfasta):
        print('DONE')
    else:
        cmd = "samtools fasta {} -F 4 > {}".format(bamFile,outfasta)
        print(cmd)
        #os.system(cmd)
        return(bamFile)

bamFiles = glob.glob("hisat2align/utr*/*.bam"); bamFile = bamFiles[0]
len(bamFiles)
with Pool(processes=4) as pool:
    toDobamFiles = pool.starmap(getBAMalignedReads, zip(bamFiles))
