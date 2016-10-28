#!usr/bin/python

import sys, os, re, math, pysam, subprocess
import numpy as np
import time

TCRbeta='chr7:141998851-142510972'
TCRalpha='chr14:22090057-23021075'
TCRgamma='chr7:38279625-38407656'

def samread():
    import pysam
	samfile=pysam.AlignmentFile("SRR.bam", "rb")
    for read in samfile.fetch('ch1', 100, 120):
	    print read
		
	pairedreads=pysam.AlignmentFile("allpaired.bam", "wb", template=samfile)	
	for read in samfile.fetch:
        if read.is_paired:
            pairedreads.write(read)
	
    samfile = pysam.AlignmentFile("ex1.bam", "rb" )
    for pileupcolumn in samfile.pileup("chr1", 100, 120):
        print ("\ncoverage at base %s = %s" %
           (pileupcolumn.pos, pileupcolumn.n))
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.
                print ('\tbase in read %s = %s' %
                    (pileupread.alignment.query_name,
                     pileupread.alignment.query_sequence[pileupread.query_position]))

    samfile.close()
	# samtools sort ex1.bam output
	pysam.sort("ex1.bam", "output")
	
	tabixfile = pysam.TabixFile("example.gtf.gz")
    for gtf in tabixfile.fetch("chr1", 1000, 2000):
         print (gtf.contig, gtf.start, gtf.end, gtf.gene_id)
	
	samfile.close()	


def ScreenGenome(fname,Locs=[]):
	if len(Locs)==0:
		return 1
	seqDir={}
	CHRlist=[]
	stlist=[]
	edlist=[]
	for Loc in Locs:
		print Loc
		tmp=Loc.split(':')
		CHR=tmp[0]
		tmp=tmp[1].split('-')
		st=int(tmp[0])
		ed=int(tmp[1])
		CHRlist.append(CHR)
		stlist.append(st)
		edlist.append(ed)
        	handle=pysam.Samfile(fname)
        	print '''Retrieve reads in the %s region''' %(Loc)
		if 'chr7' in handle.references:
			CHR=CHR
		else:
			CHR=re.sub('chr','',CHR)
        	for read in handle.fetch(CHR,st,ed):
                	if read.mapq<=30:
                        	continue
                	if read.qname not in seqDir:
                        	seqDir[read.qname]=[read]
                	else:
                        	seqDir[read.qname].append(read)
        print '''Pair reads in the region'''
        count=0
        count_r=0
	InsertSize=[]
        refs=handle.references
        for rr in handle.fetch(until_eof=True):
                if count % 1000000 ==0:
                        print "--processed %d reads" %(count)
                count+=1
		if count %10000==0:
			insize=np.fabs(rr.pnext-rr.pos)+rr.rlen
			if insize<1000:
				InsertSize.append(str(insize))
                if rr.qname in seqDir:
                        vv=seqDir[rr.qname]
                        flag=0
                        for v in vv:
                                if rr.seq==v.seq:
                                        flag=1
                        if flag==0:
                                seqDir[rr.qname].append(rr)
		else:
			if '/1' in rr.qname or '/2' in rr.qname:
                                if '/1' in rr.qname:
                                        qname_paired=re.sub('/1','/2',rr.qname)
                                else:
                                        qname_paired=re.sub('/2','/1',rr.qname)
			elif '.1' in rr.qname or '.2' in rr.qname:
                                if '.1' in rr.qname:
                                        qname_paired=re.sub('.1','.2',rr.qname)
                                else:
                                        qname_paired=re.sub('.2','.1',rr.qname)
                        else:
                                qname_paired=rr.qname
			if qname_paired in seqDir:	
                        	seqDir[qname_paired].append(rr)
                        	count_r+=1
                        	if count_r % 100 ==0:
                                #print refs[rr.rname], rr.pos
                                	print "---retrived %d unmapped reads" %(count_r)
	nfname=fname+'-Locs'+'.bam'
	HH=handle.header
        ghandle=pysam.Samfile(nfname,mode='wb',header=HH,referencenames=handle.references,referencelengths=handle.nreferences)
        for kk in seqDir:
                for read in seqDir[kk]:
                        ghandle.write(read)
	#ghandle.close()
	g=open(nfname+'.info','w')
	g.write('##'+str(time.time()))
	g.write( "##"+'+'.join(Locs)+'\n')
	g.write( "##"+str(rr.rlen)+'\n')
	g.write( "##total read number=%d, insert size=%s \n" %(count, '\t'.join(InsertSize)))
	g.close()
	return 0

def main():
	fname=sys.argv[1]
	error=ScreenGenome(fname,Locs=[TCRbeta,TCRalpha,TCRgamma])
#	error=ScreenGenome(fname,Locs=[IGH,IGK,IGL])
	if error==1:
		raise "Locus information not given!"

if __name__=='__main__':
	main()
