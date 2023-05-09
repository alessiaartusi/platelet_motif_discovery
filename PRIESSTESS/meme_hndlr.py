import glob, subprocess, os
from multiprocessing import Pool

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

############################################################
##### MEME Motifs Finding
############################################################

# def memeMotifs(fasta):
#     #/home/eidrian/Software/meme-5.4.1/src/meme mega50genesSelected.fasta -dna -searchsize 0 -oc /data/realServices/jorge/datos/geneSelection/memeMega50 -p 10
#     basemotifsdir = '/data/realServices/jorge/datos/geneSelection/motifs/'
#     outdir = os.path.join(basemotifsdir,os.path.basename(fasta).replace('.fa',''),'meme')
#     os.makedirs(outdir,exist_ok=True)
#     findMotifs = '/home/eidrian/Software/meme-5.4.1/src/meme'
#     cmd = [findMotifs,fasta,'-dna','-nmotifs','10','-searchsize','0','-oc',outdir]
#     launchjob(cmd,outdir+'/memecmd')
#     return(' '.join(map(str,cmd)))


# PRIESSTESS -fg foreground/gene_var20%plaq_exp50%mega_exp10%.fa -bg background/gene_var20%plaq_exp50%mega_exp10%.BKGR.fa -o prueba/
## PRIESSTESS streme
# streme -verbosity 1 -oc . -alph ${libpath}/alphabets/RNA_${a}_alphabet_MEME -p fg_STREME.fa -n bg_STREME.fa -pvt 0.01 -minw $min_width -maxw $max_width
#/data/realServices/jorge/datos/geneSelection/PRIESSTESS_pipeline/venv/bin/python /data/realServices/jorge/software/PRIESSTESS/bin/PRIESSTESS_logistic_regression.py LR_training_set.tab 10

def memeMotifs(fasta):
    bkgrFasta = fasta.replace("fastas","referenceFastas").replace(".fa",".BKGR.fa")
    basemotifsdir = '/data/realServices/jorge/datos/geneSelection/motifs/meme'
    outdir = os.path.join(basemotifsdir,os.path.basename(fasta).replace('.fa',''))
    if os.path.exists(os.path.join(outdir,'knownResults.html')):
        print('Already Done: '+outdir)
        return
    os.makedirs(outdir,exist_ok=True)
    findMotifs = '/home/eidrian/Software/meme-5.5.1/src/meme'
    cmd = [findMotifs,'--objfun','de',fasta,'-neg',bkgrFasta,'-dna','-searchsize','0','-nmotifs','20','-p','12','-oc',outdir]
    #launchjob(cmd,outdir+'/hommercmd')
    print('New Done: '+outdir)
    return(' '.join(map(str,cmd)))

def stremeMotifs(fasta):
    bkgrFasta = fasta.replace("fastas","referenceFastas").replace(".fa",".BKGR.fa")
    # basemotifsdir = '/data/realServices/jorge/datos/geneSelection/motifs/meme'
    basemotifsdir = '/data/realServices/jorge/datos/geneSelection/motifs/meme'
    outdir = os.path.join(basemotifsdir,os.path.basename(fasta).replace('.fa',''))
    if os.path.exists(os.path.join(outdir,'knownResults.html')):
        print('Already Done: '+outdir)
        return
    os.makedirs(outdir,exist_ok=True)
    findMotifs = '/home/eidrian/Software/meme-5.5.1/src/streme'
    cmd = [findMotifs,'--objfun','de','--p',fasta,'--n',bkgrFasta,'-dna','--nmotifs','20','--oc',outdir]
    # '--p'
    #launchjob(cmd,outdir+'/hommercmd')
    print('New Done: '+outdir)
    return(' '.join(map(str,cmd)))

fastas = glob.glob('/data/realServices/jorge/datos/geneSelection/fastas/*.fa'); fasta = fastas[0]
stremeMotifs(fastas[61])

# /home/eidrian/Software/meme-5.5.1/src/streme --objfun de --p alignedReads/utr_list_FCthresh/SRR2095014.fa -dna --nmotifs 5 --oc streme/SRR2095014
# Negative sequences are shuffled primary sequences (2-order) - training: 663021 hold-out: 73668

/home/eidrian/Software/meme-5.5.1/src/fasta-shuffle-letters -kmer 2 -dna -seed 9 -tag random alignedReads/utr_list_FCthresh/SRR2095014.fa alignedReads/utr_list_FCthresh/SRR2095014.random.fa


# PRUEBA
# /home/eidrian/Software/meme-5.5.1/src/streme --objfun de --p fastas/transcript_mega_var30%exp80%.fa -n referenceFastas/gene_var20%plaq_exp50%mega_exp10%.BKGR.fa -dna --nmotifs 5 --oc /data/realServices/jorge/datos/geneSelection/motifs/meme/gene_var20%plaq_exp50%mega_exp10%

memeMotifs(fastas[10])
# /home/eidrian/Software/meme-5.5.1/src/meme -objfun de /data/realServices/jorge/datos/geneSelection/fastas/transcript_mega_var30%exp80%.fa -neg /data/realServices/jorge/datos/geneSelection/referenceFastas/transcript_mega_var30%exp80%.BKGR.fa -dna -searchsize 0 -nmotifs 20 -p 12 -oc /data/realServices/jorge/datos/geneSelection/motifs/meme/transcript_mega_var30%exp80%


with Pool(processes=4) as pool:
    memecmds = pool.starmap(stremeMotifs, zip(fastas))


/home/eidrian/Software/meme-5.5.1/src/streme --p alignedReads/utr_list_FCthresh/SRR2095014.fa -dna --nmotifs 5 --oc streme/SRR2095014

asta-shuffle-letters [options] <sequence file> [<output file>]

## 1
# meme-5.5.1/bin/meme -objfun de geneSelection/fastas/transcript_mega_var30%exp80%.fa -neg geneSelection/referenceFastas/transcript_mega_var30%exp80%.BKGR.fa -dna -nmotifs 20 -oc geneSelection/motifs/meme/transcript_mega_var30%exp80%

## 2 With Paralell
# meme-5.5.1/bin/meme -objfun de geneSelection/fastas/gene_var50%plaq_exp80%mega_exp40%.fa -neg geneSelection/referenceFastas/gene_var50%plaq_exp80%mega_exp40%.BKGR.fa -dna -nmotifs 20 -oc geneSelection/motifs/meme/transcript_mega_var30%exp80%

meme-5.5.1/src/streme --objfun de --p geneSelection/fastas/transcript_mega_var30%exp80%.fa --n geneSelection/referenceFastas/transcript_mega_var30%exp80%.BKGR.fa -dna --nmotifs 20 --oc /data/realServices/jorge/datos/geneSelection/motifs/meme/transcript_mega_var30%exp80%

meme-5.5.1/src/streme --objfun de --p geneSelection/fastas/gene_var40%plaq_exp60%mega_exp20%.fa --n geneSelection/referenceFastas/gene_var40%plaq_exp60%mega_exp20%.BKGR.fa -dna --nmotifs 20 --oc geneSelection/motifs/meme/gene_var40%plaq_exp60%mega_exp20%
