import glob, os, subprocess
from multiprocessing import Pool
################################################################################
############# https://github.com/kaitlin309/PRIESSTESS
################################################################################

# EXAMPLE COMMAND
# PRIESSTESS -fg foreground/three_prime_utr_mega_var40%exp50%.fa -bg background/three_prime_utr_mega_var40%exp50%.BKGR.fa -o prueba/
# /data/realServices/jorge/datos/geneSelection/PRIESSTESS_pipeline/venv/bin/python /data/realServices/jorge/software/PRIESSTESS/bin/PRIESSTESS_logistic_regression.py LR_training_set.tab 10

# 1. Transform fastas to correct format and be sure it is RNA (Ts to Us)
# fastas = glob.glob("datos/geneSelection/fastas/*.fa"); fasta = fastas[0]
# for fasta in fastas:
#     outfasta = fasta.replace("fastas","PRIESSTESS_pipeline/foreground")
#     os.system("awk 'NR%2==0' {} | tr 'T' 'U' > {}".format(fasta,outfasta))
#
# fastas = glob.glob("datos/geneSelection/referenceFastas/*.fa"); fasta = fastas[0]
# for fasta in fastas:
#     outfasta = fasta.replace("referenceFastas","PRIESSTESS_pipeline/background")
#     os.system("awk 'NR%2==0' {} | tr 'T' 'U' > {}".format(fasta,outfasta))

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
    return(' '.join(cmd))

def fuseFastasOfDir(genelistFastasDir):
    outfasta = genelistFastasDir+'.fa'
    cmd = "cat {}/*.fa > {}".format(genelistFastasDir,outfasta)
    print(cmd+'\n')
    os.system(cmd)

genelistFastasDirs = glob.glob("alignedReads/utr*"); genelistFastasDir = genelistFastasDirs[0]
with Pool(processes=4) as pool:
    pool.starmap(fuseFastasOfDir, zip(genelistFastasDirs))

def createRandomDNAFastas(fasta):
    print(fasta+'\n')
    outfasta = fasta.replace('.fa','.random.fa')
    cmd = "/home/eidrian/Software/meme-5.5.1/src/fasta-shuffle-letters -line 200 -kmer 2 -dna -seed 9 -tag random {} {}".format(fasta,outfasta)
    print(cmd+'\n')
    os.system(cmd)
    print('SHUFFLING DONE \n')
    return(outfasta)

def fastas2PriesstessFormat(fasta):
    print(fasta+'\n')
    outfasta = fasta.replace('.fa','.PRIESSTESS.fa')
    cmd = "awk 'NR%2==0' {} | tr 'T' 'U' > {}".format(fasta,outfasta)
    print(cmd+'\n')
    os.system(cmd)
    print('FORMAT DONE \n')
    return(outfasta)

def PRIESSTESScmd(foregroundFasta,backgroundFasta):
    outname = os.path.basename(foregroundFasta).split('.')[0].replace('PRIESSTESS','')
    outdir = 'PRIESSTESS_pipeline/'+outname
    os.makedirs(outdir,exist_ok=True)
    cmd = ['PRIESSTESS','-fg',foregroundFasta,'-bg',backgroundFasta, '-minw', '4', '-maxw', '9', '-o', outdir]
    cmd = launchjob(cmd,outdir)
    print('PRIESSTESS Done: '+outname)
    return(cmd)

def priesstessPipeline(fasta):
    backgroundFasta = createRandomDNAFastas(fasta)
    #foregroundFasta = fastas2PriesstessFormat(fasta)
    foregroundFasta = fasta.replace('.fa','PRIESSTESS.fa')
    backgroundFasta = fastas2PriesstessFormat(backgroundFasta)
    cmd = PRIESSTESScmd(foregroundFasta,backgroundFasta)
    return(cmd)

PRIESSTESS -fg alignedReads/utr_list_logFCthresh/SRR1982704PRIESSTESS.fa -bg alignedReads/utr_list_logFCthresh/SRR1982704.random.PRIESSTESS.fa -minw 4 -maxw 9 -o PRIESSTESS_pipeline/SRR1982704PRIESSTESS

Subprocess finished
PRIESSTESS Done: SRR1982704PRIESSTESS


fastas = glob.glob("alignedReads/utr*/*[0-9].fa"); len(fastas); fasta = fastas[0]
with Pool(processes=4) as pool:
    cmds = pool.starmap(priesstessPipeline, zip(fastas))
