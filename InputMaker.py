__author__ = 'jorge'

import csv, re, argparse, os, copy, time, pickle
import sqlite3 as lite
from multiprocessing import JoinableQueue, cpu_count, Process, Queue, active_children
from operator import itemgetter

parser = argparse.ArgumentParser()
parser.add_argument('-g', action='store', dest='gene_name',
                    help='The name of the gene to be studied')
parser.add_argument('-G', action='store', dest='gene_info',
                    help='Path of the Gencode (gene annotation) gene_info file')
parser.add_argument('-f', action='store', dest='gtf',
                    help='Path of the Gencode (exons, introns, etc.) gtf file')
parser.add_argument('-p', action='store', dest='prom_ann',
                    help='Path of the FANTOM promoter annotation file')
parser.add_argument('-a', action='store', dest='association',
                    help='Path of the FANTOM enhancer-promoter association and correlation file')
parser.add_argument('-e', action='store', dest='eg_name',
                    help='Path of the EG reference file with the tissues IDs')
parser.add_argument('-x', action='store', dest='expression_dir',
                    help='Path of the directory with the Roadmap Epigenomics expression files')
parser.add_argument('-t', action='store', dest='hist_dir',
                    help='Path of the directory with the Roadmap Epigenomics '
                         'histone tagAlign files in directories named by tissue')
parser.add_argument('-m', action='store', dest='meth_dir',
                    help='Path of the directory with the Roadmap Epigenomics '
                         'methylation bedgraph files in directories named by tissue')
parser.add_argument('-D', default=float('inf'), type=int, action='store', dest='dist',
                    help='Distance value to filter promoters from the '
                         'enhancer-promoter association file. '
                         'Takes the promoters with a SMALLER distance than the one given (default: no filter)')
parser.add_argument('-C', default=0, type=float, action='store', dest='corr',
                    help='Correlation value to filter promoters from the '
                         'enhancer-promoter association file. '
                         'Takes promoters with a HIGHER correlation than the one given (default: no filter)')
parser.add_argument('-P', default=float('inf'), type=float, action='store', dest='pval',
                    help='P-value to filter promoters from the '
                         'enhancer-promoter association file. '
                         'Takes promoters with a SMALLER p-value than the one given (default: no filter)')
parser.add_argument('-F', default=float('inf'), type=float, action='store', dest='fdr',
                    help='FDR value to filter promoters from the '
                         'enhancer-promoter association file. '
                         'Takes promoters with a SMALLER FDR than the one given (default: no filter)')
parser.add_argument('-l', default=None, action='store', dest='permitted_file',
                    help='File with the list of permitted genes, otherwise it will search for it')
pargs = parser.parse_args()


def filteredDatabase(gene_info, prom_ann, association, dist, corr, pval, fdr):
    # Readers for the files
    gene_reader = csv.reader(open(gene_info, 'r'), delimiter='\t')
    prom_reader = csv.reader(open(prom_ann, 'r'), delimiter='\t')
    asso_reader = csv.reader(open(association, 'r'), delimiter='\t')

    # Clean the first lines (no information there)
    for next_line in range(1, 9):
        prom_reader.next()
        if next_line == 7:  # Just one line here
            asso_reader.next()

    # Get a dictionary of all promoters associated to a gene
    promoters = {}  # Keys are genes and values are promoter coordinates
    for pr in prom_reader:
        # In case there are more than 1 associated gene
        split = pr[1].split(',')
        coord = pr[0].split(',')[0]

        for s in split:
            found = re.search(r'p(\d)@(.*)', s)

            if bool(found):
                # Add the entries to the promoters dictionary
                if found.group(2) not in promoters.keys():
                    promoters[found.group(2)] = [coord]

                else:
                    promoters[found.group(2)].append(coord)

    # Get a dictionary of the genes in the Gencode file
    gencode_genes = {} # Keys are genes and values are Encode IDs
    for gr in gene_reader:
        if gr[6] != '' and gr[6] != 'NA':
            gencode_genes[gr[6]] = gr[0]

    # Get a list of the promoters and enhancers in the FANTOM association file with the corresponding filters
    promoter_enhancer = {}  # Promoter coordinates as keys and enhancer coordinates as values
    for ar in asso_reader:
        if float(ar[2]) < dist and float(ar[3]) > corr and float(ar[4]) < pval and float(ar[5]) < fdr:
            prom_coord = ar[1].split(',')[0]

            if prom_coord not in promoter_enhancer.keys():
                promoter_enhancer[prom_coord] = [ar[0]]

            else:
                promoter_enhancer[prom_coord].append(ar[0])

    # Make a deep copy of the promoters dictionary to iterate through it
    copy_promoters = copy.deepcopy(promoters)
    # Purge the promoters dictionary
    for key in copy_promoters.keys():
        # Remove the promoters not found in the associated promoters
        for coord_value in copy_promoters[key]:
            if coord_value not in promoter_enhancer.keys():
                promoters[key].remove(coord_value)

        # If the gene has no associated promoters delete its entry
        # Remove the promoter-associated genes not found in the Gencode file
        if len(promoters[key]) == 0 or key not in gencode_genes.keys():
            promoters.pop(key)

    # Create the txt tab delimited file with the list of genes
    filename = 'IMPermitted-D{}C{}P{}F{}.txt'.format(dist, corr, pval, fdr)
    with open(filename, 'a') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for gene in promoters.keys():
            # reduce() joins the list of lists of associated enhancers for each promoters associated to a gene
            enhancers = reduce(lambda x, y: x + y, [promoter_enhancer[prom] for prom in promoters[gene]])
            # Columns in the file will be Encode ID - gene name - promoters - enhancers
            writer.writerow([gencode_genes[gene]] + [gene] + [','.join(promoters[gene])] + [','.join(enhancers)])


class Producer(object):

    def __init__(self, gene, permitted, hist_queue, meth_queue):
        self.gtf_file = pargs.gtf
        self.gene = gene
        self.permitted = permitted
        self.hist_queue = hist_queue
        self.meth_queue = meth_queue

        # Start the processes
        pr_en_process = Process(target=self.promoters_enhancers)
        pr_en_process.start()
        ex_in_process = Process(target=self.exons_introns)
        ex_in_process.start()

    def promoters_enhancers(self):
        # Search for the gene name in the permitted genes file
        print 'Opening {} file'.format(self.permitted)
        with open(self.permitted, 'r') as perm_file:
            permitted_reader = csv.reader(perm_file, delimiter='\t')
            for line in permitted_reader:
                if line[1] == self.gene:
                    # Split the promoters and enhancers
                    promoters = line[2].split(',')
                    enhancers = line[3].split(',')

                    # Put every promoter in the queues
                    for pr in promoters:
                        print 'Sending promoters to the queues'
                        self.hist_queue.put({'promoter': pr})
                        self.meth_queue.put({'promoter': pr})

                    # Put every enhancer in the queues
                    for en in enhancers:
                        print 'Sending enhancers to the queues'
                        # Change the coordinate format from chr#:#-# to chr#:#..#
                        self.hist_queue.put({'enhancer': en.replace('-', '..')})
                        self.meth_queue.put({'enhancer': en.replace('-', '..')})

                    break

    def exons_introns(self):
        # Search for the gene name in the gtf file
        print 'Opening {} file'.format(self.gtf_file)
        with open(self.gtf_file) as handle:
            gtf_reader = csv.reader(handle, delimiter='\t')
            for line in gtf_reader:
                # The gene name is in the 8th column of each line 4th entry of the semicolon separated string
                gtf_info = line[8].split(';')
                gene_name = re.search(r' gene_name \"(.*)\"', gtf_info[4]).group(1)

                # If there is a match put the gene's elements in the queues with the corresponding label
                if gene_name == self.gene:
                    # For introns
                    if line[2] == 'intron':
                        print 'Sending introns to the queues'
                        coordinates = '{}:{}..{}'.format(line[0], line[3], line[4])
                        self.hist_queue.put({'intron': coordinates})
                        self.meth_queue.put({'intron': coordinates})

                    # For exons
                    elif line[2] == 'exon' and re.search(r' transcript_type \"(.*)\"',
                                                         gtf_info[5]).group(1) != 'retained_intron':
                        print 'Sending exons to the queues'
                        coordinates = '{}:{}..{}'.format(line[0], line[3], line[4])
                        self.hist_queue.put({'exon': coordinates})
                        self.meth_queue.put({'exon': coordinates})

                    # For retained introns
                    elif line[2] == 'exon' and re.search(r'transcript_type \"(.*)\"',
                                                         gtf_info[5]).group(1) == 'retained_intron':
                        print 'Sending retained introns to the queues'
                        coordinates = '{}:{}..{}'.format(line[0], line[3], line[4])
                        self.hist_queue.put({'retained_intron': coordinates})
                        self.meth_queue.put({'retained_intron': coordinates})


class Consumer(object):
    def __init__(self, hist_queue, meth_queue, res_queue):
        self.hist_dir = pargs.hist_dir
        self.meth_dir = pargs.meth_dir
        self.hist_queue = hist_queue
        self.meth_queue = meth_queue
        self.results_queue = res_queue

        # Start the processes
        hist_process = Process(target=self.hist_worker)
        hist_process.start()
        meth_process = Process(target=self.meth_worker)
        meth_process.start()

    def hist_worker(self):
        while True:
            print 'Is hist_queue empty?'.format(self.hist_queue.empty())
            # Get a job and break in case it is a poison pill
            job = self.hist_queue.get()
            print 'Starting with job {}'.format(job)
            if job is None:
                self.hist_queue.task_done()
                break

            # There should only be one key in the dictionary from the queue
            element_type = job.keys()[0]
            coords = job[element_type]
            element_info = re.search(r'(.*):(\d*)\.\.(\d*)', coords)

            # Extract the element information
            chromosome = element_info.group(1)
            element_start = float(element_info.group(2))
            element_end = float(element_info.group(3))

            # The folders in the histones directory should have the tissues names/ids
            hist_walker = os.walk(self.hist_dir)

            # Walk through the directories
            for h_dir_path, h_dirs, h_files in hist_walker:

                print 'Walking through {}'.format(self.hist_dir)
                # Search for tagAlign files in each folder
                for track in h_files:
                    found_tagAlign = re.search(r'(.*)-(.*)\.tagAlign', track)

                    # If found, read the file
                    if bool(found_tagAlign):
                        print 'Found file {}'.format(track)
                        file_path = h_dir_path + '/' + track
                        hist_modification = found_tagAlign.group(2)

                        with open(file_path) as handle:
                            tagAlign_reader = csv.reader(handle, delimiter='\t')

                            # It will count the amount of modifications found in the element
                            modification_counter = 0

                            for tr in tagAlign_reader:
                                # Check if the chromosome is the right one
                                if tr[0] == chromosome:
                                    hist_start = float(tr[1])
                                    hist_end = float(tr[2])

                                    # If the modification is inside or in one of the borders of the element count it
                                    if (element_start <= hist_start < hist_end <= element_end or
                                        hist_start < element_start <= hist_end or
                                            hist_start <= element_end < hist_end):

                                        # Add to the modification counter
                                        print 'Adding to the {} counter'.format(hist_modification)
                                        modification_counter += 1

                            # Make a list with the results
                            print 'Got results for {}'.format(hist_modification)
                            results = [element_type, coords, hist_modification, modification_counter]

                            # Put the results in the queue as a dictionary
                            tissue = h_dir_path.split('/')[-1]
                            self.results_queue.put({tissue: results})
            print 'Task done with job {}'.format(job)
            self.hist_queue.task_done()

    def meth_worker(self):
        while True:
            print 'Is meth queue empty? {}'.format(self.meth_queue.empty())
            job = self.meth_queue.get()
            print 'Starting with job {}'.format(job)
            if job is None:
                self.meth_queue.task_done()
                break

            # There should only be one key in every dictionary from the queue
            element_type = job.keys()[0]
            coords = job[element_type]
            element_info = re.search(r'(.*):(\d*)\.\.(\d*)', coords)

            # Extract the element information
            chromosome = element_info.group(1)
            element_start = float(element_info.group(2))
            element_end = float(element_info.group(3))

            # The folders in the methylation directory should have the tissues names/ids
            meth_walker = os.walk(self.meth_dir)

            # Walk through the directories
            for m_dir_path, m_dirs, m_files in meth_walker:

                print 'Walking through {}'.format(self.meth_dir)
                # Look for a bedgraph file in each folder
                for track in m_files:
                    found_bedgraph = re.search(r'(.*)\.bedgraph', track)

                    if bool(found_bedgraph):
                        file_path = m_dir_path + '/' + track
                        print 'Found file {}'.format(track)

                        # Read the file
                        with open(file_path) as handle:
                            bedgraph_reader = csv.reader(handle, delimiter='\t')

                            # List and counter to get the mean of the methylation scores
                            meth_score_sum = []
                            meth_score_counter = 0.0

                            for br in bedgraph_reader:
                                # Check that the chromosome is the right one
                                if br[0] == chromosome:
                                    meth_start = float(br[1])
                                    meth_end = float(br[2])
                                    meth_score = float(br[3])

                                    # If the methylated region is inside the element or in any border of the element
                                    if (element_start <= meth_start < meth_end <= element_end or
                                        meth_start < element_start <= meth_end or
                                            meth_start <= element_end < meth_end):

                                        # Add the score to the list
                                        print 'Calculating sums for {}'.format(job)
                                        meth_score_sum.append(meth_score)
                                        meth_score_counter += 1

                            # Calculate the average score
                            if meth_score_counter != 0:
                                meth_score_average = sum(meth_score_sum) / meth_score_counter

                            else:
                                meth_score_average = 0

                            # Make a list with the results
                            print 'Got results for {}'.format(job)
                            results = [element_type, coords, 'Methylation', meth_score_average]

                            # Put the results in the queue as a dictionary
                            tissue = m_dir_path.split('/')[-1]  # The containing folder has the name of the tissue
                            self.results_queue.put({tissue: results})

                        # It will only use the first bedgraph file found
                        break
            print 'Ending task with job {}'.format(job)
            self.meth_queue.task_done()

def expression(tissue, gene, permitted, expr_dir):
    expression_score = 0
    # Walk through the expression directory
    walker = os.walk(expr_dir)
    for dir_path, dirs, files in walker:
        for track in files:
            found_pc = re.search(r'(.*)\.pc', track)
            found_nc = re.search(r'(.*)\.nc', track)
            found_rb = re.search(r'(.*)\.rb', track)

            if bool(found_pc):
                file_path = dir_path + '/' + track
                with open(file_path) as handle:
                    pc_reader = csv.reader(handle, delimiter='\t')

            elif bool(found_nc):
                file_path = dir_path + '/' + track
                with open(file_path) as handle:
                    nc_reader = csv.reader(handle, delimiter='\t')

            elif bool(found_rb):
                file_path = dir_path + '/' + track
                with open(file_path) as handle:
                    rb_reader = csv.reader(handle, delimiter='\t')

if __name__ == '__main__':

    start_time = time.time()

    # Create the required queues
    histone_queue = JoinableQueue()
    methylation_queue = JoinableQueue()
    results_queue = Queue()

    # Search for the permitted genes file if it is None
    permitted_file = pargs.permitted_file
    if permitted_file is None:
        print 'Searching for the permitted file'
        permitted_walker = os.walk('./')
        for p_dir_path, p_dirs, p_files in permitted_walker:
            for file_name in p_files:
                found_permitted = re.search(r'IMPermitted-D(.*)C(.*)P(.*)F(.*)\.txt', file_name)

                if bool(found_permitted):
                    permitted_file = file_name
                    break

        # If the permitted file is not found, make one and use it
        if permitted_file is None:
            print 'File not found, creating one'
            filteredDatabase(pargs.gene_info, pargs.prom_ann, pargs.association,
                             pargs.dist, pargs.corr, pargs.pval, pargs.fdr)

            permitted_file = 'IMPermitted-D{}C{}P{}F{}.txt'.format(pargs.dist, pargs.corr, pargs.pval, pargs.fdr)

    # The amount of real processes will be cpu_count() * 2 (see the Consumer class)
    consumer_processes = cpu_count()

    # Start the producers
    Producer(pargs.gene_name, permitted_file, histone_queue, methylation_queue)
    print 'Creating {} Consumers'.format(consumer_processes * 2)
    # Start the consumers
    for i in xrange(consumer_processes):
        Consumer(histone_queue, methylation_queue, results_queue)

    print active_children()

    # Wait for the consumers to stop working
    histone_queue.join()
    methylation_queue.join()
    print 'Queues joined'

    # Insert the poison pills in the queues
    for i in xrange(consumer_processes):
        histone_queue.put(None)
        methylation_queue.put(None)

    # Wait for the poison pills to make effect
    histone_queue.join()
    methylation_queue.join()

    print (time.time() - start_time) / 60.0, 'minutes'

    # # Get the results and organize them
    # results_dictionary = {} # Keys are tissues and values are lists with a dictionary for the element types
    # while True:
    #     element_result = results_queue.get()
    #     tissue = element_result.keys()[0]
    #     element_info = element_result[tissue]
    #     element_type = element_info[0]
    #     entry = (element_info[1], element_info[2], element_info[3])
    #
    #     # If the tissue entry doesn't exist yet in the dictionary
    #     if tissue not in results_dictionary.keys():
    #         results_dictionary[tissue] = {element_type: [entry]}
    #     # If the tissue entry exists in the dictionary but the element type is missing
    #     elif tissue in results_dictionary.keys() and element_type not in results_dictionary[tissue].keys():
    #         results_dictionary[tissue][element_type] = [entry]
    #
    #     # If the tissue entry exists in the dictionary and the element type exists in this entry
    #     elif tissue in results_dictionary.keys() and element_type in results_dictionary[tissue].keys():
    #         # Append the new element entry to the dictionary
    #         results_dictionary[tissue][element_type].append(entry)
    #         # Sort by coordinates and then by feature (methylation and histone modifications)
    #         sorted_list = sorted(results_dictionary[tissue][element_type], key=itemgetter(0, 1))
    #         # Change the list of entries in the dictionary for the sorted one
    #         results_dictionary[tissue][element_type] = sorted_list
    #
    #         ''' Since this will be done for each tissue, all the elements are organized in the same way
    #             for each one of them. This way the columns will be the same when the database is created '''
    #
    #     if results_queue.empty():
    #         break
    #
    # ''' Esto esta muy mal hecho, mejorarlo '''
    # output_file = '{}:IMOutput-D{}C{}P{}F{}.txt'.format(pargs.gene_name, pargs.dist, pargs.corr, pargs.pval, pargs.fdr)
    # with open(output_file, 'a') as handle:
    #     writer = csv.writer(handle, delimiter='\t')
    #
    #     elements = results_dictionary['Spleen'].keys()
    #     for tissue in results_dictionary.keys():
    #         tissue_rows = [tissue]
    #         for element in elements:
    #             for column in results_dictionary[tissue][element]:
    #                 column_entry = '{}@{}-{}={}'.format(element[0:2], column[0], column[1], str(column[2]))
    #                 tissue_rows.append(column_entry)
    #
    #         writer.writerow(tissue_rows)

    # # Create the SQLite database
    # db_file = 'IMOutput-D{}C{}P{}F{}.db'.format(pargs.dist, pargs.corr, pargs.pval, pargs.fdr)
    # connection = lite.connect(db_file)
    #
    # tissue = results_dictionary.keys()[0]
    # elements = results_dictionary[tissue].keys()
    #
    #
    # elements = ' '.join()
    # schema = 'Tissue TEXT Gene_name TEXT Gene_ID TEXT Coord TEXT Expr REAL ' + elements
    #
    # with connection:
    #
    #     db_cursor = connection.cursor()
    #     db_cursor.execute('DROP TABLE IF EXISTS Elements')
    #     db_cursor.execute('CREATE TABLE Elements({})'.format(schema))
    #
    #     for tissue in results_dictionary.keys():
    #         results_dictionary[tissue]
    #
    #
    # with open('Tissue_dict.pickle', 'wb') as handle:
    #     pickle.dump(results_dictionary, handle)





  #   Traceback (most recent call last):
  # File "InputMaker.py", line 424, in <module>
  #   methylation_queue.join()
  # File "/usr/lib/python2.7/multiprocessing/queues.py", line 340, in join
  #   self._cond.wait()
  # File "/usr/lib/python2.7/multiprocessing/synchronize.py", line 246, in wait
  #   self._wait_semaphore.acquire(True, timeout)

