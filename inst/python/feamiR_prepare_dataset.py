import os
import argparse, gzip
from Bio import SeqIO
from gtfparse import read_gtf
import pandas as pd
import numpy as np
import math
import scipy.stats as ss
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

def complement(seq):
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A','N':'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)
def reverse_complement(s):
    return complement(s[::-1])
def reverse(s):
    return (s[::-1])
def reverse_complement_row(row):
    if row['strand']=='-':
        return reverse_complement(row['sequence'])
    else:
        return str(row['sequence'])

def seedmatch(seed,mR):
	if (seed==mR):
		return True
	else:
		if (seed[1:]==mR[1:]):
			return True
		if (seed[:-1]==mR[:-1]):
			return True
		else:
			return False

def findlastmatch(mirseq,mracomparable):
	currentlast=0
	for i in range(min((len(mracomparable),len(mirseq)))):
		if mirseq[i]==mracomparable[i]:
			currentlast=i
	return currentlast


def findfirstmatch(mirseq,mracomparable):
	for i in range(min(len(mirseq),len(mracomparable))):
		if mirseq[i]==mracomparable[i]:
			return i

def numberateachposition(dataset,i):
    return (len(dataset[dataset['seed_sequence'].str[i] == 'A']),len(dataset[dataset['seed_sequence'].str[i] == 'G']),len(dataset[dataset['seed_sequence'].str[i] == 'C']),len(dataset[dataset['seed_sequence'].str[i] == 'U']))
def numberateachposition_fullmiRNA(dataset,i):
    return (len(dataset[dataset['full_miRNA_sequence'].str[i] == 'A']),len(dataset[dataset['full_miRNA_sequence'].str[i] == 'G']),len(dataset[dataset['full_miRNA_sequence'].str[i] == 'C']),len(dataset[dataset['full_miRNA_sequence'].str[i] == 'U']))
def numberateachposition_upstream_flanking(dataset,i):
    return (len(dataset[dataset['upstream_flanking'].str[i] == 'A']),len(dataset[dataset['upstream_flanking'].str[i] == 'G']),len(dataset[dataset['upstream_flanking'].str[i] == 'C']),len(dataset[dataset['upstream_flanking'].str[i] == 'U']))
def numberateachposition_downstream_flanking(dataset,i):
    return (len(dataset[dataset['downstream_flanking'].str[i] == 'A']),len(dataset[dataset['downstream_flanking'].str[i] == 'G']),len(dataset[dataset['downstream_flanking'].str[i] == 'C']),len(dataset[dataset['downstream_flanking'].str[i] == 'U']))

def outof100_1nt(A,G,C,U):
    total=A+G+C+U
    Aprop = int(100*(A/total)+0.5)
    Gprop = int(100*(G/total)+0.5)
    Cprop = int(100*(C/total)+0.5)
    Uprop = int(100*(U/total)+0.5)
    maximum = max([Aprop,Gprop,Cprop,Uprop])
    if (maximum==Aprop):
        Aprop = 100 - (Gprop + Cprop + Uprop)
    elif (maximum==Gprop):
        Gprop = 100 - (Aprop + Cprop + Uprop)
    elif (maximum==Cprop):
        Cprop = 100 - (Aprop + Gprop + Uprop)
    else:
        Uprop = 100 - (Aprop + Gprop + Cprop)
    return [Aprop,Gprop,Cprop,Uprop]
def outof100posi(A,G,C,U,i):
    return outof100_1nt(A[i],G[i],C[i],U[i])
def outof100_i(function,dataset,i):
    return outof100_1nt(function(dataset,i)[0],function(dataset,i)[1],function(dataset,i)[2],function(dataset,i)[3])
def compareoutof100_i(function,dataset1,dataset2,i):
    return pd.DataFrame({'pos':outof100_i(function,dataset1,i),'neg':outof100_i(function,dataset2,i)},index=['A','G','C','U'])

def numberateachposition2nt_i(dataset,i):
    l = []
    for j in ['A','G','C','U']:
        for k in ['A','G','C','U']:
            l.append(len(dataset[np.logical_and(dataset['miRNA_seq'].str[i] == j,dataset['miRNA_seq'].str[i+1] == k)]))
    return l
def numberateachposition2nt_i_fullmiRNA(dataset,i):
    l = []
    for j in ['A','G','C','U']:
        for k in ['A','G','C','U']:
            l.append(len(dataset[np.logical_and(dataset['full_miRNA_sequence'].str[i] == j,dataset['full_miRNA_sequence'].str[i+1] == k)]))
    return l
def numberateachposition2nt_i_upstream_flanking(dataset,i):
    l = []
    for j in ['A','G','C','U']:
        for k in ['A','G','C','U']:
            l.append(len(dataset[np.logical_and(dataset['upstream_flanking'].str[i] == j,dataset['upstream_flanking'].str[i+1] == k)]))
    return l
def numberateachposition2nt_i_downstream_flanking(dataset,i):
    l = []
    for j in ['A','G','C','U']:
        for k in ['A','G','C','U']:
            l.append(len(dataset[np.logical_and(dataset['downstream_flanking'].str[i] == j,dataset['downstream_flanking'].str[i+1] == k)]))
    return l

def outof100_i_2nt(attribute,dataset,i):
    l = []
    for j in ['A','G','C','U']:
        for k in ['A','G','C','U']:
            app = len(dataset[np.logical_and(dataset[attribute].str[i] == j,dataset[attribute].str[i+1] == k)])
            if app < 0:
                app = 0
            l.append(app)
    l = [round(x * (100/sum(l))) for x in l]
    l[-1]=100-sum(l[:len(l)-1])
    if l[-1]<0:
        l[-2]=l[-1]+l[-2]
        l[-1]=0

    return l
def compare_100_2nt_i(attribute,dataset1,dataset2,i):
    return pd.DataFrame({'pos':outof100_i_2nt(attribute,dataset1,i),'neg':outof100_i_2nt(attribute,dataset2,i)},index=['AA','AG','AC','AU','GA','GG','GC','GU','CA','CG','CC','CU','UA','UG','UC','UU'])
def fishertable(table,col1,col2,letter):
    length1 = sum(table[col1])
    length2 = sum(table[col2])
    return pd.DataFrame([[table[col1][letter],table[col2][letter]],[length1-table[col1][letter],length2-table[col2][letter]]],columns=[col1,col2])
def fishertable_2nt(table,col1,col2,letter1,letter2):
    length1 = sum(table[col1])
    length2 = sum(table[col2])
    return pd.DataFrame([[table[col1][letter1+letter2],table[col2][letter1+letter2]],[length1-table[col1][letter1+letter2],length2-table[col2][letter1+letter2]]],columns=[col1,col2])

def fisher_1nt(positive,negative,numberofpos,function):
    dfforheat = pd.DataFrame([],columns=['chisq','fisherA','fisherG','fisherC','fisherU'])
    for i in range(numberofpos):
        l = []
        comp1 = compareoutof100_i(function,positive,negative,i)
        comp = comp1[(comp1['pos']!=0) | (comp1['neg']!=0)]
        a,p = ss.chisquare(f_obs=comp['pos'],f_exp=comp['neg'])
        l.append(p)
        for j in ['A','G','C','U']:
            if comp1['pos'][j] == 0 and comp1['neg'][j]==0:
                    l.append(1)
            else:
                lettertab = fishertable(comp,'pos','neg',j)
                if lettertab['pos'][0]<0:
                    lettertab['pos'][0]=0
                    lettertab['pos'][1]=100
                b,q = ss.fisher_exact(lettertab)
                l.append(q)
        dfforheat.loc[i] = l
    return dfforheat

def fisher_2nt(positive,negative,numberofpos,attribute):
    dfforheat_2nt = pd.DataFrame([],columns=['chisq','fisherAA','fisherAG','fisherAC','fisherAU','fisherGA','fisherGG','fisherGC','fisherGU','fisherCA','fisherCG','fisherCC','fisherCU','fisherUA','fisherUG','fisherUC','fisherUU'])
    for i in range(numberofpos):
        l = []
        comp1 = compare_100_2nt_i(attribute,positive,negative,i)
        comp = comp1[(comp1['pos']!=0) | (comp1['neg']!=0)]
        a,p = ss.chisquare(f_obs=comp['pos'],f_exp=comp['neg'])
        l.append(p)
        for j in ['A','G','C','U']:
            for k in ['A','G','C','U']:
                if comp1['pos'][j+k] == 0 and comp1['neg'][j+k]==0:
                    l.append(1)
                else:
                    lettertab = fishertable_2nt(comp,'pos','neg',j,k)
                    if lettertab['pos'][0]<0:
                        lettertab['pos'][0]=0
                        lettertab['pos'][1]=100
                    b,q = ss.fisher_exact(lettertab)
                    l.append(q)
        dfforheat_2nt.loc[i] = l
    return dfforheat_2nt

def subsample_neg(positivity,negativity,threshold,j):
    dupsampledneg = negativity.sample(n=len(positivity),random_state=j)
    sampledneg = dupsampledneg.drop_duplicates(subset=['miRNA_id','gene_id'],keep='first').reset_index().drop(labels=['index'],axis=1)
    for pos in range(8):
        sample1A,sample1G,sample1C,sample1U = numberateachposition(sampledneg,pos)
        samplelist = outof100_1nt(sample1A,sample1G,sample1C,sample1U)
        pos1A,pos1G,pos1C,pos1U = numberateachposition(negativity,pos)
        neglist = outof100_1nt(pos1A,pos1G,pos1C,pos1U)
        f1,samplepvalue1 = ss.chisquare(f_obs=samplelist,f_exp=neglist)

        if (samplepvalue1 < threshold):
            representative=False
            j=j+1
        else:
            representative=True
        if representative==False:
            return False   
    if representative==False:
        return (False,'')
    else:
        return (True,sampledneg)

def subsample(positive,negative,output_prefix,num_runs):
    positivity=positive[['miRNA_id','gene_id','seed_sequence']].drop_duplicates()

    j=1
    for k in range(num_runs):
        representative=False
        while representative==False:
            a,subsampledneg = subsample_neg(positive,negative,0.1,j)
            if a==False:
                j=j+1
            else:
                dummies_pos = pd.DataFrame()
                mir_single = positivity['seed_sequence'].apply(lambda x: pd.Series(list(x)))
                for i in mir_single.columns:
                    dummies_pos['seed_position'+str(i+1)]=mir_single[i]
                    if i!=mir_single.columns[-1]:
                        dummies_pos['seed_pair'+str(i+1)]=mir_single[i]+mir_single[i+1]
                dummies_pos['classification']=1
                dummies_neg = pd.DataFrame()
                mir_single = subsampledneg['seed_sequence'].apply(lambda x: pd.Series(list(x)))
                for i in mir_single.columns:
                    dummies_neg['seed_position'+str(i+1)]=mir_single[i]
                    if i!=mir_single.columns[-1]:
                        dummies_neg['seed_pair'+str(i+1)]=mir_single[i]+mir_single[i+1]
                dummies_neg['classification']=0
                collist = [x for x in dummies_pos.columns if (x!='miRNA_id' and x!='gene_id' and x!='classification')]
                dummies_dummies_pos = pd.get_dummies(dummies_pos,columns=collist)
                collist = [x for x in dummies_neg.columns if (x!='miRNA_id' and x!='gene_id' and x!='classification')]
                dummies_dummies_neg = pd.get_dummies(dummies_neg,columns=collist)
                combined_dummies = dummies_dummies_neg.append(dummies_dummies_pos,ignore_index=True)
                combined_dummies.to_csv(output_prefix+str(k)+'.csv')
                j=j+1
                representative=True


#Require 3'UTR file with id as gene name (or whatever is the same as interaction files). Call this file mRNA_3pUTR
#Require miRNA file with id as miRNA identified (same as interaction file). Call this file miRNA_full
#Require interaction dataset with miRNA column, target_gene column and some additional columns if wanted.
def main():

    global args
#First reformat the 3'UTR file and mIRNA file
    parser = argparse.ArgumentParser()
    parser.add_argument('-mRNA_3pUTR', metavar='<mRNA_3pUTR>', dest="mRNA_3pUTR", help="fasta file containing 3'UTR sequences for species in question, with ID matching the mRNA ID in the interaction dataset")
    parser.add_argument('-miRNA_full', metavar='<miRNA_full>', dest="miRNA_full", help="fasta file containing full microRNA sequences, with ID matching the ID matching the microRNA ID in the interaction dataset")
    parser.add_argument('-interactions', metavar='<interactions>', dest="interactions", help="csv file containing miRNA ids (column labelled miRNA_id),mRNA target genes (column labelled target_gene)")
    parser.add_argument('-annotations', metavar='<annotations>', dest="annotations", help="gtf annotations file used to identify 3'UTRs")
    parser.add_argument('-fullchromosomes', metavar='<fullchromosomes>', dest="fullchromosomes", help="fasta file containing full sequence for each chromosome for the species in question (the chromosomes must be labelled consistently with the gtf file)")
    parser.add_argument('-seed', metavar='<seed>', dest="seed", help="1 if seed features should be used in analysis,default:1")
    parser.add_argument('-nonseed_miRNA', metavar='<nonseed_miRNA>', dest="nonseed_miRNA", help="1 if nonseed miRNA features should be used in analysis,default:0")
    parser.add_argument('-flankingmRNA', metavar='<flankingmRNA>', dest="flankingmRNA", help="1 if mRNA flanking features should be used in analysis,default:0")    
    parser.add_argument('-UTR_output', metavar='<UTR_output_file>', dest="UTR_output_file", help="file name for 3'UTR sequences to be saved in")    
    parser.add_argument('-chr', metavar='<num_chr>', dest="num_chr", help="number of chromosomes for the species in question")
    parser.add_argument('-o', metavar='<output_prefix>', dest="output_prefix", help="prefix for any output files")
    parser.add_argument('-positiveset', metavar='<positiveset>', dest="positiveset", help="existing validated set of interactions. This should only be a file output by previous steps in this pipeline.")
    parser.add_argument('-negativeset', metavar='<negativeset>', dest="negativeset", help="existing set of non-validated interactions. This should only be a file output by previous steps in this pipeline.")
    parser.add_argument('-sreformatpath', metavar='<sreformatpath>', dest="sreformatpath", help="Path where sreformat is installed, including sreformat itself. sreformat must be installed to use this pipeline.")
    parser.add_argument('-patmanpath', metavar='<patmanpath>', dest="patmanpath", help="Path where PaTMaN is installed, including patman itself. sreformat must be installed to use this pipeline.")  
    parser.add_argument('-patmanoutput', metavar='<patmanoutput>', dest="patmanoutput", help="Processed csv file from patman alignment. This should only be a file output by previous steps in this pipeline.")  
    parser.add_argument('-minvalidationentries', metavar='<minvalidationentries>', dest="minvalidationentries", help="Minimum number of entries for a validation category to be considered separately (default: 40).")  
    parser.add_argument('-num_runs', metavar='<num_runs>', dest="num_runs", help="Number of subsamples to create, default: 100")  


    args = parser.parse_args()
    if args.output_prefix:
        output_prefix = args.output_prefix
    else:
        output_prefix = ''
        
    if (args.positiveset and args.negativeset):
        status = 'pos-neg'
    elif (args.patmanoutput and args.interactions):
        status = 'patman'
    elif (args.miRNA_full and args.mRNA_3pUTR):
        status = 'utrs extracted'
    elif (args.miRNA_full and args.mRNA_fullchromosomes):
        status = 'extract utr'
    else:
        logging.info('Exit - incorrect inputs given.')
        return 1

    if args.minvalidationentries:
        minvalidationentries = int(args.minvalidationentries)
    else:
        minvalidationentries = 40
    if (status =='extract utr' or status == 'utrs extracted') :
        if not args.interactions:
            logging.info('No interaction dataset provided')
            return 1
        else:
            interaction_dataset = pd.read_csv(args.interactions)
            try:
                interaction_dataset = pd.read_csv(args.interactions)
                interaction_dataset = interaction_dataset[['miRNA','Target Gene','Experiments','Support Type']]
                exp_type='support and type'
                logging.info('Extracted columns of interactions dataset are miRNA, Target Gene, Experiments and Support Type')
            except:
                try:
                    interaction_dataset = pd.read_csv(args.interactions)[['miRNA','Target Gene','Experiments']]
                    exp_type = 'type'
                    logging.info('Extracted columns of interactions dataset are miRNA, Target Gene and Experiments')
                except:
                    try:
                        interaction_dataset = pd.read_csv(args.interactions)[['miRNA','Target Gene','Support Type']]
                        exp_type = 'support'
                        logging.info('Extracted columns of interactions dataset are miRNA, Target Gene and Support Type')
                    except:
                        try:
                            interaction_dataset = pd.read_csv(args.interactions)[['miRNA','Target Gene']]
                            exp_type = 'short'
                            logging.info('Extracted columns of interactions dataset are miRNA and Target Gene')
                        except:
                            logging.info('Interaction dataset of incorrect form. Please provide a csv with columns miRNA, Target Gene and optionally Experiments (containing an experiment type) and/or Support.Type (values Functional MTI, Functional MTI (Weak), Non-Functional MTI, Non-Functional MTI (Weak)')
                            return(1)
        interaction_dataset = interaction_dataset.rename(columns={'miRNA':'miRNA_id','Target Gene':'gene_id'})

        if (status=='extract utrs'):
            
            if args.UTR_output_file:
                prime_output_file = output_prefix+args.UTR_output_file
            else:
                prime_output_file = output_prefix+'UTR_seq.fasta'
            logging.info('Annotation file '+args.annotations+' and full chromosome file '+args.fullchromosomes+' supplied. Extracted three prime UTRs will be saved to '+prime_output_file)
            # combine gtf files and fasta file
            try:
                mRNA_gtf = read_gtf(args.annotations)
                mRNA_gtf = mRNA_gtf[mRNA_gtf['feature']=='three_prime_utr'][['seqname','start','end','gene_name','gene_id','transcript_id','strand']]
            except:
                logging.info('Failed to read gtf annotations file. Please check format.')
                return(1)
            try:
                with open(args.fullchromosomes) as seqfile: 
                    mRNA_fasta = [x for x in SeqIO.parse(seqfile, 'fasta')]
            except:
                logging.info('Failed to read fullchromosomes file')
                return(1)
            chrlist = [str(i) for i in range(1,int(args.num_chr))]
            logging.info('Currently extracting three prime UTRs from fullchromosomes file '+args.fullchromosomes)
            p=0
            for i in chrlist:
                logging.info('Currently on chromosome' + i)
                ongoingseries = pd.Series(dtype='string')
                #df is already filtered to just 3'UTR
                currentchr = mRNA_gtf[mRNA_gtf['seqname']==i]
                currentchr = currentchr.sort_values(by=['transcript_id','start'])
                #freqdf is a frequency table for transcript_ids
                freqdf = currentchr.transcript_id.value_counts()
                #take the right sequence for that chromosome
                for row in range(len(mRNA_fasta)):
                    if mRNA_fasta[row].name == i:
                        sequence = mRNA_fasta[row].seq
                #Take sequences between start and stop for each entry in threeprime and add as column seqpart
                list_of_seqs = list()
                for l in range(len(currentchr)):
                    start = int(currentchr.iloc[l]['start'])-1
                    stop = int(currentchr.iloc[l]['end'])
                    seqpart = sequence[start:stop]
                    list_of_seqs.append(seqpart)
                currentchr['seqpart'] = list_of_seqs
                #This is combining the sequences for multiple entries with the same transcript into a series with transcript_id and full combined sequence
                if len(list_of_seqs) > 0:
                    p=p+1
                    k=0
                    while k<len(currentchr):
                        tid = currentchr.iloc[k]['transcript_id']
                        currentlength = freqdf[tid]
                        newseq = ''
                        for j in range(currentlength):
                            newseq = newseq + currentchr.iloc[k]['seqpart']
                            k=k+1
                        newaddition = pd.Series([newseq],index=[tid])
                        ongoingseries = ongoingseries.append(newaddition)

                    #Take list of unique transcript_ids
                    uniquetids = currentchr[['transcript_id','gene_id','gene_name','strand']].drop_duplicates()
                    ongoingdf = ongoingseries.to_frame()
                    ongoingdf = ongoingdf.reset_index()
                    ongoingdf.columns = ['transcript_id','sequence']
                    #Combine sequences with the gene_name, gene_id and transcript_name and strand
                    mergified = uniquetids.merge(ongoingdf)
                    #Reverse complement the - ones
                    mergified['final_sequence']=mergified.apply(reverse_complement_row,axis=1)
                    mergified = mergified[['gene_name','final_sequence','transcript_id']]
                    mergified = mergified.drop_duplicates()
                    if i == '1':
                        with open(prime_output_file,'w') as primefile:
                            for i in range(len(mergified)):
                                print(">" + mergified.iloc[i]['gene_name'],file=primefile)
                                print(mergified.iloc[i]['final_sequence'], file=primefile)
                    else:
                        with open(prime_output_file,'a') as primefile:
                            for i in range(len(mergified)):
                                print(">" + mergified.iloc[i]['gene_name'],file=primefile)
                                print(mergified.iloc[i]['final_sequence'], file=primefile)

            if p==0:
                logging.info('Failed to extract three prime UTRs. Please check chromosomes are labelled by 1,2,3 etc and that annotations file contains three prime UTRs.')
                return(1)
            chrlist = currentchr = currentlength = freqdf = None
            i = j = k = l = p=None
            list_of_seqs = mRNA_fasta = mRNA_gtf = None
            mergified = newaddition = newseq = ongoingdf = None
            ongoingseries = parser = primefile = None
            row = seqfile = seqpart = sequence = start = stop = tid = uniquetids = None
            mRNAfile_name = prime_output_file
            #Initial processing to get 3'UTR file and save as some string (input?)
        elif args.mRNA_3pUTR:
            mRNAfile_name = args.mRNA_3pUTR
            logging.info('Three prime UTR file '+mRNAfile_name+' supplied so three prime UTR extraction step skipped.')
        else:
            logging.info('Exit - no mRNA file given. Please give either a full chromosomes (e.g. toplevel) file with corresponding anonotation gtf file or a fasta file containing only three prime UTR sequences')
            return 1
        try:
            with open(mRNAfile_name) as mrna_file_open: 
                fullmRNAlist = [x for x in SeqIO.parse(mrna_file_open, 'fasta')]
            idlist = []
            seqlist = []
            for i in range(len(fullmRNAlist)):
                idlist.append(fullmRNAlist[i].name)
                seqlist.append(str(fullmRNAlist[i].seq))
            fullmRNAdf = pd.DataFrame({'gene_id':idlist,'mRNA_sequence':seqlist})
        except:
            logging.info('Failed to open mRNA file. Please check format of input files.')
            return 1
        idlist = None
        seqlist = None
        fullmRNAlist = None



        logging.info('Creating miRNA seed file and saving to '+output_prefix+'seed.fasta')
        try:
            with open(args.miRNA_full) as miRNA_full_open: 
                fulllist = [x for x in SeqIO.parse(miRNA_full_open, 'fasta')]
                with open(output_prefix+'seed.fasta',"w") as output_miRNAseed_open:
                    for i in range(len(fulllist)):
                        print(">"+fulllist[i].id,file=output_miRNAseed_open)
                        print(fulllist[i].seq[:8],file=output_miRNAseed_open)
        except:
            logging.info("Unable to open miRNA file. Please check format.")
            return 1
        seedfile_name = output_prefix+'seed.fasta'
        miRNAfile_name = args.miRNA_full
        if args.mRNA_3pUTR:
            logging.info('Reformating mRNA file and miRNA files to remove hidden characters and check consistent T/Us and capitalisation for PaTMaN alignment. Saving files to '+output_prefix+"mRNA_reformated.fasta"+', '+output_prefix+"reformated_seed.fasta"+' and '+output_prefix+"reformated_fullmiRNA.fasta")
            os.system(args.sreformatpath+' -u -r fasta '+mRNAfile_name +" > "+output_prefix+"mRNA_reformated.fasta")
            mRNAfile_name = output_prefix+"mRNA_reformated.fasta"           
        else:
            logging.info('Reformating mRNA file and miRNA files to remove hidden characters and check consistent T/Us and capitalisation for PaTMaN alignment. Saving files to '+output_prefix+"mRNA_reformated.fasta"+', '+output_prefix+"reformated_seed.fasta"+' and '+output_prefix+"reformated_fullmiRNA.fasta")
            os.system(args.sreformatpath+' -u -r fasta '+mRNAfile_name +" > "+output_prefix+"mRNA_reformated.fasta")
            mRNAfile_name = output_prefix+"mRNA_reformated.fasta"
        os.system(args.sreformatpath+' -u -r fasta '+miRNAfile_name+" > "+output_prefix+"reformated_fullmiRNA.fasta")
        os.system(args.sreformatpath+' -u -r fasta '+seedfile_name+" > "+output_prefix+"reformated_seed.fasta")
        seedfile_name = output_prefix+"reformated_seed.fasta"
        miRNAfile_name = output_prefix+"reformated_fullmiRNA.fasta"

        if (not args.seed or (not args.seed == 0)):
            with open(seedfile_name) as seed_file_open: 
                fullseedlist = [x for x in SeqIO.parse(seed_file_open, 'fasta')]
            idlist = []
            seqlist = []
            for i in range(len(fullseedlist)):
                idlist.append(fullseedlist[i].name)
                seqlist.append(str(fullseedlist[i].seq))
            fullseeddf = pd.DataFrame({'miRNA_id':idlist,'seed_sequence':seqlist})
            fullseedlist = idlist = seqlist = seed_file_open = None
            with open(miRNAfile_name) as miRNA_file_open: 
                fullmiRNAlist = [x for x in SeqIO.parse(miRNA_file_open, 'fasta')]
            idlist = []
            seqlist = []
            for i in range(len(fullmiRNAlist)):
                idlist.append(fullmiRNAlist[i].name)
                seqlist.append(str(fullmiRNAlist[i].seq))
            fullmiRNAdf = pd.DataFrame({'miRNA_id':idlist,'full_miRNA_sequence':seqlist})
            fullmiRNAlist = idlist = seqlist = miRNA_file_open = None

            logging.info('Running PaTMaN seed alignment.')
            os.system(args.patmanpath+' -e 2 -g 0 -D '+mRNAfile_name+' -P '+seedfile_name+' -o '+output_prefix+'patman_seed.txt')
            totalseed = pd.read_csv(output_prefix+'patman_seed.txt',sep='\t',header=None)
            totalseed.columns = ['gene_id','miRNA_id','start','stop','strand','num_mismatches']
            totalseed = totalseed[totalseed['strand']=='-']
            totalseed.miRNA_id = totalseed.miRNA_id.str.rstrip()
            totalseed.gene_id = totalseed.gene_id.str.rstrip()
            logging.info('Preparing to process PaTMaN output.')
            totalseed = totalseed.merge(fullseeddf,on=['miRNA_id'])
            totalseed = totalseed.merge(fullmiRNAdf,on=['miRNA_id'])
            totalseed = totalseed.merge(fullmRNAdf,on=['gene_id'])
            totalseed = totalseed[['miRNA_id','gene_id','start','stop','seed_sequence','full_miRNA_sequence','mRNA_sequence','num_mismatches']]
            totalseed['comparable']=''
            totalseed['matchtype']=''
            totalseed['upstream_flanking']=''
            totalseed['downstream_flanking']=''
            listofcomparable = []
            listofmatchtype = []
            listofupstream_flanking = []
            listofdownstream_flanking = []
            numrows = len(totalseed)
            numrows_over_20 = numrows//20
            logging.info('Processing PaTMaN output, filtering out false seed matches and adding full feature set')
            for i in range(numrows):
                if (i % numrows_over_20 == 0):
                    logging.info(str(5*(i//numrows_over_20))+'%'+' of PaTMaN output processed.')
                start = int(totalseed['start'][i])-1
                stop = int(totalseed['stop'][i])
                seed = totalseed['seed_sequence'][i]
                fullmiRNA = totalseed['full_miRNA_sequence'][i]
                mRNA = totalseed['mRNA_sequence'][i]
                mRNA_rc = reverse_complement(mRNA)
                if totalseed['num_mismatches'][i] == 0 or totalseed['num_mismatches'][i] == '0':
                    listofcomparable.append(seed)
                    listofmatchtype.append('Exact')
                else:
                    if start<0 or stop>=len(totalseed['mRNA_sequence'][i]):
                        listofmatchtype.append('Reject')
                        listofcomparable.append('N')
                    else:
                        comparable = reverse_complement(mRNA[start:stop])
                        listofcomparable.append(comparable)
                        if len(comparable)==8 and len(seed)==8 and (comparable[1:7])==(seed[1:7]):
                            if comparable==seed:
                                listofmatchtype.append('Exact')
                            elif comparable[7]==seed[7] and comparable[0]=='U':
                                listofmatchtype.append('C')
                            elif comparable[0]=='U':
                                listofmatchtype.append('A')
                            elif comparable[7] == seed[7]:
                                listofmatchtype.append('B')
                            else:
                                listofmatchtype.append('Reject')
                        else:
                            listofmatchtype.append('Reject')
                fullcomparable = reverse_complement(mRNA[stop-len(fullmiRNA):stop])
                firstmatch = findfirstmatch(fullmiRNA,fullcomparable)
                lastmatch = findlastmatch(fullmiRNA,fullcomparable)
                if firstmatch is None:
                    listofdownstream_flanking.append('N')
                else:
                    listofdownstream_flanking.append(mRNA[max(0,stop-firstmatch):min(stop-firstmatch+20,len(mRNA)-1)])
                if lastmatch is None:
                    listofupstream_flanking.append('N')
                else:
                    listofupstream_flanking.append(reverse(mRNA[max(0,stop-lastmatch-1-20):min(stop-lastmatch-1,len(mRNA)-1)]))
        totalseed['matchtype']=listofmatchtype
        totalseed['comparable']=listofcomparable
        totalseed['upstream_flanking']=listofupstream_flanking
        totalseed['downstream_flanking']=listofdownstream_flanking
        listofcomparable = listofmatchtype = listofupstream_flanking = listofdownstream_flanking = None
        totalseed = totalseed[totalseed['matchtype']!='Reject']
        totalseed.to_csv(output_prefix+'patmanprocessed.csv')
        interaction_dataset['positive']=1
        interaction_labelled = totalseed.merge(interaction_dataset,how='left')
        positiveset_seed = interaction_labelled[interaction_labelled['positive']==interaction_labelled['positive']].drop(columns=['positive'])
        negativeset_seed = interaction_labelled[interaction_labelled['positive']!=interaction_labelled['positive']].drop(columns=['positive'])
        try:
            negativeset_seed = negativeset_seed.drop(columns=['Experiments'])
        except:
            pass
        try:
            negativeset_seed = negativeset_seed.drop(columns=['Support Type'])
        except:
            pass
        positiveset_seed.to_csv(output_prefix+'seed_positive.csv')
        negativeset_seed.to_csv(output_prefix+'seed_negative.csv')
        logging.info('PaTMaN output processed and split into positive and negative datasets, save at '+output_prefix+'seed_positive.csv and '+output_prefix+'seed_negative.csv')
    elif args.patmanoutput:
        logging.info('PaTMaN output supplied so skipping alignment.')
        if not args.interactions:
            logging.info('no interaction dataset provided')
            return 1
        else:
            interaction_dataset = pd.read_csv(args.interactions)
            try:
                interaction_dataset = pd.read_csv(args.interactions)
                interaction_dataset = interaction_dataset[['miRNA','Target Gene','Experiments','Support Type']]
                exp_type='support and type'
                logging.info('Extracted columns of interactions dataset are miRNA, Target Gene, Experiments and Support Type')
            except:
                try:
                    interaction_dataset = pd.read_csv(args.interactions)[['miRNA','Target Gene','Experiments']]
                    exp_type = 'type'
                    logging.info('Extracted columns of interactions dataset are miRNA, Target Gene and Experiments')
                except:
                    try:
                        interaction_dataset = pd.read_csv(args.interactions)[['miRNA','Target Gene','Support Type']]
                        exp_type = 'support'
                        logging.info('Extracted columns of interactions dataset are miRNA, Target Gene and Support Type')
                    except:
                        try:
                            interaction_dataset = pd.read_csv(args.interactions)[['miRNA','Target Gene']]
                            exp_type = 'short'
                            logging.info('Extracted columns of interactions dataset are miRNA and Target Gene')
                        except:
                            logging.info('Interaction dataset of incorrect form. Please provide a csv with columns miRNA, Target Gene and optionally Experiments (containing an experiment type) and/or Support Type (values Functional MTI, Functional MTI (Weak), Non-Functional MTI, Non-Functional MTI (Weak)')
                            return(1)
        interaction_dataset = interaction_dataset.rename(columns={'miRNA':'miRNA_id','Target Gene':'gene_id'})
        interaction_dataset['positive']=1
        try:
            totalseed = pd.read_csv(args.patmanoutput)
        except:
            logging.info('Failed to open PaTMaN output file.')
            return(1)
        interaction_labelled = totalseed.merge(interaction_dataset,how='left')
        positiveset_seed = interaction_labelled[interaction_labelled['positive']==interaction_labelled['positive']].drop(columns=['positive'])
        negativeset_seed = interaction_labelled[interaction_labelled['positive']!=interaction_labelled['positive']].drop(columns=['positive'])
        try:
            negativeset_seed = negativeset_seed.drop(columns=['Experiments'])
        except:
            pass
        try:
            negativeset_seed = negativeset_seed.drop(columns=['Support Type'])
        except:
            pass
        positiveset_seed.to_csv(output_prefix+'seed_positive.csv')
        negativeset_seed.to_csv(output_prefix+'seed_negative.csv')
        logging.info('PaTMaN output processed and split into positive and negative datasets, save at '+output_prefix+'seed_positive.csv and '+output_prefix+'seed_negative.csv')



    else:
        logging.info('Positive and negative datasets provided so initial processing and PaTMaN alignment is skipped.')
        try:
            positiveset_seed = pd.read_csv(args.positiveset).drop(columns=['Unnamed: 0'])
            negativeset_seed = pd.read_csv(args.negativeset).drop(columns=['Unnamed: 0'])
        except:
            try:
                positiveset_seed = pd.read_csv(args.positiveset)
                negativeset_seed = pd.read_csv(args.negativeset)
            except:
                logging.info('Failed to open positive and/or negative input files. Please check format.')
        if (positiveset_seed.columns).contains('Experiments'):
            if (positiveset_seed.columns).contains('Support Type'):
                exp_type = 'support and type'
            else:
                exp_type = 'type'
        else:
            if (positiveset_seed.columns).contains('Support Type'):
                exp_type = 'support'
            else:
                exp_type = 'short'
    if exp_type == 'support' or exp_type == 'support and type':
        functional = positiveset_seed[positiveset_seed['Support Type']=='Functional MTI']
        if len(functional)>100:
            logging.info('Positive file filtered to only Functional MTI. Number of remaining entries is: '+str(len(functional)))
            positiveset_seed = functional
            functional = None


    positiveset_file = output_prefix+'seed_positive.csv'
    negativeset_file = output_prefix+'seed_negative.csv'

    #dataset prepared

    #statistical results
    os.system('mkdir '+output_prefix+'statistical_analysis')
    #seed
    try:
        seed_dfforheat = fisher_1nt(positiveset_seed,negativeset_seed,8,numberateachposition)
        seed_dfforheat_2nt = fisher_2nt(positiveset_seed,negativeset_seed,7,'seed_sequence')
        seed_dfforheat.to_csv(output_prefix+'statistical_analysis/seed_1nt_statsummary.csv')
        seed_dfforheat_2nt.to_csv(output_prefix+'statistical_analysis/seed_2nt_statsummary.csv')
        heat = sns.heatmap(seed_dfforheat,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
        heat.savefig(output_prefix+'statistical_analysis/seed_1nt_heatmap.jpg', bbox_inches = "tight")
        heat_2nt = sns.heatmap(seed_dfforheat_2nt,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
        heat_2nt.savefig(output_prefix+'statistical_analysis/seed_2nt_heatmap.jpg', bbox_inches = "tight")

        if args.nonseed_miRNA and (args.nonseed_miRNA == 1 or args.nonseed_miRNA == '1'):
            fullmiRNA_dfforheat = fisher_1nt(positiveset_seed,negativeset_seed,22,numberateachposition_fullmiRNA)
            fullmiRNA_dfforheat_2nt = fisher_2nt(positiveset_seed,negativeset_seed,21,'full_miRNA_sequence')
            fullmiRNA_dfforheat.to_csv(output_prefix+'statistical_analysis/fullmiRNA_1nt_statsummary.csv')
            fullmiRNA_dfforheat_2nt.to_csv(output_prefix+'statistical_analysis/fullmiRNA_2nt_statsummary.csv')
            heat = sns.heatmap(fullmiRNA_dfforheat,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
            heat.savefig(output_prefix+'statistical_analysis/fullmiRNA_1nt_heatmap.jpg', bbox_inches = "tight")
            heat_2nt = sns.heatmap(fullmiRNA_dfforheat_2nt,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
            heat_2nt.savefig(output_prefix+'statistical_analysis/fullmiRNA_2nt_heatmap.jpg', bbox_inches = "tight")

        if args.flankingmRNA and (args.flankingmRNA == 1 or args.flankingmRNA == '1'):
            upstreamflanking_dfforheat = fisher_1nt(positiveset_seed,negativeset_seed,20,numberateachposition_upstream_flanking)
            upstreamflanking_dfforheat_2nt = fisher_2nt(positiveset_seed,negativeset_seed,19,'upstream_flanking')
            upstreamflanking_dfforheat.to_csv(output_prefix+'statistical_analysis/upstreamflanking_1nt_statsummary.csv')
            upstreamflanking_dfforheat_2nt.to_csv(output_prefix+'statistical_analysis/upstreamflanking_2nt_statsummary.csv')
            heat = sns.heatmap(upstreamflanking_dfforheat,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
            heat.savefig(output_prefix+'statistical_analysis/upstreamflanking_1nt_heatmap.jpg', bbox_inches = "tight")
            heat_2nt = sns.heatmap(upstreamflanking_dfforheat_2nt,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
            heat_2nt.savefig(output_prefix+'statistical_analysis/upstreamflanking_2nt_heatmap.jpg', bbox_inches = "tight")

            downstreamflanking_dfforheat = fisher_1nt(positiveset_seed,negativeset_seed,20,numberateachposition_downstream_flanking)
            downstreamflanking_dfforheat_2nt = fisher_2nt(positiveset_seed,negativeset_seed,19,'downstream_flanking')
            downstreamflanking_dfforheat.to_csv(output_prefix+'statistical_analysis/downstreamflanking_1nt_statsummary.csv')
            downstreamflanking_dfforheat_2nt.to_csv(output_prefix+'statistical_analysis/downstreamflanking_2nt_statsummary.csv')
            heat = sns.heatmap(downstreamflanking_dfforheat,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
            heat.savefig(output_prefix+'statistical_analysis/downstreamflanking_1nt_heatmap.jpg', bbox_inches = "tight")
            heat_2nt = sns.heatmap(downstreamflanking_dfforheat_2nt,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
            heat_2nt.savefig(output_prefix+'statistical_analysis/downstreamflanking_2nt_heatmap.jpg', bbox_inches = "tight")
    except:
        logging.info('Oops sorry something went wrong with the statistical analysis. Please check the inputs or try again!')
        return(1)
    os.system('mkdir '+output_prefix+'subsamples')
    #subsampling for classifiers

    logging.info('Creating '+str(args.num_runs)+' representative subsamples')
    try: 
        subsample(positiveset_seed,negativeset_seed,output_prefix+'subsamples/full_subsample',int(args.num_runs))
    except:
        logging.info('Problem with creating subsamples. Please check the num_runs parameter is valid and try again')
        return(1)

    if (exp_type == 'type') or (exp_type == 'support and type'):
        totalexperiments = positiveset_seed['Experiments'].drop_duplicates()
        listofexperiments = []
        for exp in totalexperiments:
            if len(positiveset_seed[positiveset_seed['Experiments']==exp])>=minvalidationentries:
                listofexperiments.append(exp)
        logging.info('Considered experiment types are:')
        logging.info(listofexperiments)
        if len(listofexperiments) == 0:
            logging.info('No validation categories large enough to perform analysis. Please relabel and try again')
        else:
            for exp in listofexperiments:
                positive_exp = positiveset_seed[positiveset_seed['Experiments']==exp]
                os.system('mkdir '+output_prefix+'statistical_analysis/'+exp)
                try:
                    seed_dfforheat = fisher_1nt(positive_exp,negativeset_seed,8,numberateachposition)
                    seed_dfforheat_2nt = fisher_2nt(positive_exp,negativeset_seed,7,'seed_sequence')
                    seed_dfforheat.to_csv(output_prefix+'statistical_analysis/’+exp+’/seed_1nt_statsummary.csv')
                    seed_dfforheat_2nt.to_csv(output_prefix+'statistical_analysis/’+exp+’/seed_2nt_statsummary.csv')
                    heat = sns.heatmap(seed_dfforheat,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
                    heat.savefig(output_prefix+'statistical_analysis/’+exp+’/seed_1nt_heatmap.jpg', bbox_inches = "tight")
                    heat_2nt = sns.heatmap(seed_dfforheat_2nt,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
                    heat_2nt.savefig(output_prefix+'statistical_analysis/’+exp+’/seed_2nt_heatmap.jpg', bbox_inches = "tight")

                    if args.nonseed_miRNA and (args.nonseed_miRNA == 1 or args.nonseed_miRNA == '1'):
                        fullmiRNA_dfforheat = fisher_1nt(positive_exp,negativeset_seed,22,numberateachposition_fullmiRNA)
                        fullmiRNA_dfforheat_2nt = fisher_2nt(positive_exp,negativeset_seed,21,'full_miRNA_sequence')
                        fullmiRNA_dfforheat.to_csv(output_prefix+'statistical_analysis/’+exp+’/fullmiRNA_1nt_statsummary.csv')
                        fullmiRNA_dfforheat_2nt.to_csv(output_prefix+'statistical_analysis/’+exp+’/fullmiRNA_2nt_statsummary.csv')
                        heat = sns.heatmap(fullmiRNA_dfforheat,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
                        heat.savefig(output_prefix+'statistical_analysis/’+exp+’/fullmiRNA_1nt_heatmap.jpg', bbox_inches = "tight")
                        heat_2nt = sns.heatmap(fullmiRNA_dfforheat_2nt,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
                        heat_2nt.savefig(output_prefix+'statistical_analysis/’+exp+’/fullmiRNA_2nt_heatmap.jpg', bbox_inches = "tight")

                    if args.flankingmRNA and (args.flankingmRNA == 1 or args.flankingmRNA == '1'):
                        upstreamflanking_dfforheat = fisher_1nt(positive_exp,negativeset_seed,20,numberateachposition_upstream_flanking)
                        upstreamflanking_dfforheat_2nt = fisher_2nt(positive_exp,negativeset_seed,19,'upstream_flanking')
                        upstreamflanking_dfforheat.to_csv(output_prefix+'statistical_analysis/’+exp+’/upstreamflanking_1nt_statsummary.csv')
                        upstreamflanking_dfforheat_2nt.to_csv(output_prefix+'statistical_analysis/’+exp+’/upstreamflanking_2nt_statsummary.csv')
                        heat = sns.heatmap(upstreamflanking_dfforheat,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
                        heat.savefig(output_prefix+'statistical_analysis/’+exp+’/upstreamflanking_1nt_heatmap.jpg', bbox_inches = "tight")
                        heat_2nt = sns.heatmap(upstreamflanking_dfforheat_2nt,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
                        heat_2nt.savefig(output_prefix+'statistical_analysis/’+exp+’/upstreamflanking_2nt_heatmap.jpg', bbox_inches = "tight")

                        downstreamflanking_dfforheat = fisher_1nt(positive_exp,negativeset_seed,20,numberateachposition_downstream_flanking)
                        downstreamflanking_dfforheat_2nt = fisher_2nt(positive_exp,negativeset_seed,19,'downstream_flanking')
                        downstreamflanking_dfforheat.to_csv(output_prefix+'statistical_analysis/’+exp+’/downstreamflanking_1nt_statsummary.csv')
                        downstreamflanking_dfforheat_2nt.to_csv(output_prefix+'statistical_analysis/’+exp+’/downstreamflanking_2nt_statsummary.csv')
                        heat = sns.heatmap(downstreamflanking_dfforheat,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
                        heat.savefig(output_prefix+'statistical_analysis/’+exp+’/downstreamflanking_1nt_heatmap.jpg', bbox_inches = "tight")
                        heat_2nt = sns.heatmap(downstreamflanking_dfforheat_2nt,vmin=0,vmax=0.2,cmap='RdBu',center=0.1).get_figure()
                        heat_2nt.savefig(output_prefix+'statistical_analysis/’+exp+’/downstreamflanking_2nt_heatmap.jpg', bbox_inches = "tight")
                except:
                    logging.info('We had a problem with processing experiment type '+exp+' . Please have a look at that experiment type.')

if __name__ == "__main__":
    main()