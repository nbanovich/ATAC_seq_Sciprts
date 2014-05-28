import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages
import pysam
import sys
import pdb

MIN_READ_LEN = 1
MAX_READ_LEN = 1000

MIN_MAP_QUAL = 10
chrom_names = ['chr%d'%i for i in xrange(1,23)]

def plot_footprint(top_profile, bottom_profile, title):

    figure = plot.figure()
    subplot = figure.add_subplot(111)
    L = top_profile.size
    xvals = np.arange(-L/2,L/2)
    subplot.plot(xvals, top_profile, color='r', label='top 20%')
    subplot.plot(xvals, bottom_profile, color='b', label='bottom 20%')

    subplot.axis([xvals.min(), xvals.max(), 0, max([top_profile.max(),bottom_profile.max()])])

    subplot.set_xlabel('position relative to motif')
    subplot.set_ylabel('atacseq footprint signal')

    legend = subplot.legend(loc=1)
    for text in legend.texts:
        text.set_fontsize(6)
    legend.set_frame_on(False)

    subplot.set_title(title)

    return figure

def plot_tss_centered_signal(tss_profile, title):

    figure = plot.figure()
    subplot = figure.add_subplot(111)
    L = tss_profile.size
    xvals = np.arange(-L/2,L/2)
    subplot.plot(xvals, tss_profile, color='r', linewidth=0.5)

    subplot.axvline(0, linestyle='--', linewidth=0.2)
    subplot.axis([xvals.min(), xvals.max(), 0, tss_profile.max()])

    subplot.set_xlabel('position relative to TSS')
    subplot.set_ylabel('atacseq signal')

    subplot.set_title(title)

    return figure

if __name__=="__main__":

    inputbam = sys.argv[1]
    outputpdf = sys.argv[2]
    title = inputbam.split('/')[-1]

    # load TSS from bed file
    tss_handle = open("/mnt/lustre/home/anilraj/immune/dat/refGene_hg19.txt",'r')
    tss_sites = []
    for line in tss_handle:
        row = line.strip().split()
        if row[1] not in chrom_names:
            continue
        if row[2]=='+':
            tss_sites.append((row[1],int(row[3]),row[2]))
        else:
            tss_sites.append((row[1],int(row[4]),row[2]))
    tss_handle.close()
    tss_sites = list(set(tss_sites))
    S = len(tss_sites)

    # load CTCF high chipseq sites
    handle = open("/mnt/lustre/home/anilraj/immune/dat/Ctcf_high_chipseq_locations.bed",'r')
    top_sites = []
    for line in handle:
        row = line.strip().split()
        top_sites.append((row[0],int(row[1]),int(row[2]),row[3]))
    handle.close()
    T = len(top_sites)

    # load CTCF low chipseq sites
    handle = open("/mnt/lustre/home/anilraj/immune/dat/Ctcf_low_chipseq_locations.bed",'r')
    bottom_sites = []
    for line in handle:
        row = line.strip().split()
        bottom_sites.append((row[0],int(row[1]),int(row[2]),row[3]))
    handle.close()
    B = len(bottom_sites)

    #pdb.set_trace()

    N = 600

    pdfhandle = PdfPages(outputpdf)
    top_profile = np.zeros((N,),dtype=int)
    bottom_profile = np.zeros((N,),dtype=int)
    tss_profile = np.zeros((4000,),dtype=int)

    bamfile = pysam.Samfile(inputbam,'rb')

    for site in top_sites:

        profile = np.zeros((N,),dtype=int)
        left = int(site[1])-N/2
        right = int(site[1])+N/2
        reads = bamfile.fetch(site[0], left, right-1)
        for read in reads:

            # reads appear twice, once for each side, only want to consider once
            if not (read.is_read1 or read.is_read2):
                continue

            # require that both sides are uniquely mapped
            if read.is_unmapped or read.mate_is_unmapped:
                continue

            if read.is_reverse == read.mate_is_reverse:
                # reads mapped to same strand...
                continue

            if read.mapq < MIN_MAP_QUAL:
                # read has poor mapping quality
                continue

            # remember pysam pos starts at 0, not 1
            if read.is_reverse:
                isize = -read.isize
                start = read.mpos + 1
            else:
                isize = read.isize
                start = read.pos + 1

            if isize < MIN_READ_LEN or isize > MAX_READ_LEN:
                continue

            end = start + isize - 1

            if site[3]=='+':
                if start+4-left>=0 and start+4-left<N:
                    profile[start+4-left] += 1
                if end-5-left>=0 and end-5-left<N:
                    profile[end-5-left] += 1
            else:
                if right-end+5>=0 and right-end+5<N:
                    profile[right-end+5] += 1
                if right-start-4>=0 and right-start-4<N:
                    profile[right-start-4] += 1

        top_profile += profile

    for site in bottom_sites:

        profile = np.zeros((N,),dtype=float)
        left = int(site[1])-N/2
        right = int(site[1])+N/2
        reads = bamfile.fetch(site[0], left, right-1)
        for read in reads:

            # reads appear twice, once for each side, only want to consider once
            if not (read.is_read1 or read.is_read2):
                continue

            # require that both sides are uniquely mapped
            if read.is_unmapped or read.mate_is_unmapped:
                continue

            if read.is_reverse == read.mate_is_reverse:
                # reads mapped to same strand...
                continue

            if read.mapq < MIN_MAP_QUAL:
                # read has poor mapping quality
                continue

            # remember pysam pos starts at 0, not 1
            if read.is_reverse:
                isize = -read.isize
                start = read.mpos + 1
            else:
                isize = read.isize
                start = read.pos + 1

            if isize < MIN_READ_LEN or isize > MAX_READ_LEN:
                continue

            end = start + isize - 1
    
            if site[3]=='+':
                if start+4-left>=0 and start+4-left<N:
                    profile[start+4-left] += 1
                if end-5-left>=0 and end-5-left<N:
                    profile[end-5-left] += 1
            else:
                if right-end+5>=0 and right-end+5<N:
                    profile[right-end+5] += 1
                if right-start-4>=0 and right-start-4<N:
                    profile[right-start-4] += 1

        bottom_profile += profile

    for site in tss_sites:

        profile = np.zeros((4000,),dtype=int)
        left = int(site[1])-2000
        right = int(site[1])+2000
        reads = bamfile.fetch(site[0], left, right-1)
        for read in reads:

            # reads appear twice, once for each side, only want to consider once
            if not (read.is_read1 or read.is_read2):
                continue

            # require that both sides are uniquely mapped
            if read.is_unmapped or read.mate_is_unmapped:
                continue

            if read.is_reverse == read.mate_is_reverse:
                # reads mapped to same strand...
                continue

            if read.mapq < MIN_MAP_QUAL:
                # read has poor mapping quality
                continue

            # remember pysam pos starts at 0, not 1
            if read.is_reverse:
                isize = -read.isize
                start = read.mpos + 1
            else:
                isize = read.isize
                start = read.pos + 1

            if isize < MIN_READ_LEN or isize > MAX_READ_LEN:
                continue

            end = start + isize - 1

            if site[2]=='+':
                if start+4-left>=0 and start+4-left<4000:
                    profile[start+4-left] += 1
                if end-5-left>=0 and end-5-left<4000:
                    profile[end-5-left] += 1
            else:
                if right-end+5>=0 and right-end+5<4000:
                    profile[right-end+5] += 1
                if right-start-4>=0 and right-start-4<4000:
                    profile[right-start-4] += 1

        tss_profile += profile

    bamfile.close()

    top_profile = top_profile/float(T)
    bottom_profile = bottom_profile/float(B)
    figure = plot_footprint(top_profile, bottom_profile, title)
    pdfhandle.savefig(figure)
    tss_profile = tss_profile/float(S)
    figure = plot_tss_centered_signal(tss_profile, title)
    pdfhandle.savefig(figure)

    pdfhandle.close()
